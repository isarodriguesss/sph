#!/usr/bin/env python3
"""Hook PreToolUse (Edit|Write): bloqueia mudanca de parametro fisico sem
leitura fresca de log.csv nesta sessao.

Regra (CLAUDE.md Sec.10, Proibicoes; Sec.2.1): "Nao alterar parametros fisicos
sem antes ler log.csv". Este hook so dispara quando o novo conteudo do arquivo
muda o valor numerico de um parametro fisico conhecido (ver PARAM_NAMES,
espelhando a tabela da Sec.7) em main.py/src/equations.py/src/scheme.py.

Consome (invalida) a marca de leitura apos liberar uma mudanca -- a proxima
mudanca de parametro exige uma nova leitura de log.csv (mesma cadencia do
protocolo Sec.10: ler -> prever -> mudar UMA variavel -> rodar -> comparar).

Escotilha: comentario "# calibracao-sem-log-ok" no trecho novo libera a
mudanca sem leitura previa (uso consciente -- ex.: revert documentado,
mudanca nao-fisica que colide com o regex, primeira calibracao da sessao
com predicao ja discutida na conversa).

Fail-closed (exit 2). Sem dependencias externas.
"""

import json
import re
import sys
from pathlib import Path

STATE_FILE = Path(__file__).resolve().parent.parent / ".state" / "log_read.json"
ESCAPE_MARKER = "# calibracao-sem-log-ok"

WATCHED_PATHS = ("main.py", "src/equations.py", "src/scheme.py")

# Espelha a tabela de parametros calibrados da Sec.7 do CLAUDE.md.
PARAM_NAMES = [
    "beta",
    "sigma",
    "k_consume",
    "D_ext",
    "D_int",
    "D_n",
    "D_n_int",
    "D_o",
    "lambda_",
    "lambda_ext_ratio",
    "lambda_ext",
    "r_growth",
    "f0",
    "tension_ratio",
    "c0",
    "mu",
    "gamma",
    "gamma_mature",
    "alpha_mon",
    "k_n",
    "K_n",
    "cs_max",
    "x_dim",
    "y_dim",
    "h_factor",
    "rho_max",
    "k_o",
    "Q0",
]

NUM_RE = r"(-?\d+\.?\d*(?:[eE][-+]?\d+)?)"


def extract_value(name, text):
    pattern = rf"(?<![A-Za-z0-9_]){re.escape(name)}(?![A-Za-z0-9_])\s*[:=]\s*{NUM_RE}"
    m = re.search(pattern, text)
    return m.group(1) if m else None


def is_watched(rel):
    p = rel.replace("\\", "/")
    return any(p == w or p.endswith("/" + w) for w in WATCHED_PATHS)


def changed_params(tool_name, tool_input):
    old_text = tool_input.get("old_string", "") or ""
    new_text = tool_input.get("new_string", "") or ""
    if tool_name == "Write":
        new_text = tool_input.get("content", "") or ""
        old_text = ""

    changed = []
    for name in PARAM_NAMES:
        new_val = extract_value(name, new_text)
        if new_val is None:
            continue
        old_val = extract_value(name, old_text)
        if old_val is None or old_val != new_val:
            changed.append((name, old_val, new_val))
    return changed


def session_has_fresh_read(session_id):
    if not STATE_FILE.exists():
        return False
    try:
        state = json.loads(STATE_FILE.read_text(encoding="utf-8"))
    except Exception:
        return False
    entry = state.get(session_id)
    if not entry or entry.get("consumed"):
        return False
    return True


def consume_read(session_id):
    try:
        state = json.loads(STATE_FILE.read_text(encoding="utf-8"))
    except Exception:
        return
    if session_id in state:
        state[session_id]["consumed"] = True
        STATE_FILE.write_text(json.dumps(state, indent=2), encoding="utf-8")


def main():
    try:
        data = json.load(sys.stdin)
    except Exception:
        return

    tool_name = data.get("tool_name", "")
    if tool_name not in ("Write", "Edit"):
        return

    tool_input = data.get("tool_input", {}) or {}
    fpath = tool_input.get("file_path") or tool_input.get("path")
    if not fpath:
        return
    rel = Path(fpath).as_posix()
    if not is_watched(rel):
        return

    full_new = tool_input.get("new_string") or tool_input.get("content") or ""
    if ESCAPE_MARKER in full_new:
        return

    changed = changed_params(tool_name, tool_input)
    if not changed:
        return

    session_id = data.get("session_id", "unknown")
    if session_has_fresh_read(session_id):
        consume_read(session_id)
        return

    lines = [f"  - {n}: {o!r} -> {v!r}" for n, o, v in changed]
    msg = (
        "BLOQUEADO: mudanca de parametro fisico sem leitura fresca de log.csv nesta sessao.\n\n"
        "Parametro(s) detectado(s):\n" + "\n".join(lines) + "\n\n"
        "CLAUDE.md Sec.10 (Proibicoes) exige ler log.csv (Sec.2.1: tendencias de mean_v, "
        "n_fast, a_marangoni, a_flag, mean_cs, max_cs, contrast_cs) e formular uma predicao "
        "quantitativa (Sec.2.3) ANTES de qualquer calibracao de parametro fisico.\n\n"
        "Leia log.csv primeiro (Read ou skill analisar-log-csv), ou se esta mudanca ja foi "
        "justificada por leitura/predicao discutida nesta conversa, adicione o comentario "
        f'"{ESCAPE_MARKER}" no trecho editado.'
    )
    print(msg, file=sys.stderr)
    sys.exit(2)


if __name__ == "__main__":
    main()
