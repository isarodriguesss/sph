#!/usr/bin/env python3
"""Hook PreToolUse (Edit|Write): bloqueia introducao de rugosidade de contorno
(Pass L) em arquivos .py enquanto o pre-requisito morfologico nao for
satisfeito.

Regra (CLAUDE.md Sec.1 Objetivo 1, Sec.2.2, Sec.10 Proibicoes): Pass L (paredes
com topografia irregular) so pode ser proposto apos a simulacao reproduzir a
morfologia dendritica de reference.jpg / reference_result.png painel (b) via
mecanismos hidrodinamicos puros. Rugosidade e refinamento fisico, nao muleta.

Dispara por palavra-chave (heuristica, nao semantica) em codigo novo tocando
geometria de contorno. Escotilha: comentario "# pass-l-aprovado: <motivo>"
libera (ex.: usuario decidiu explicitamente desbloquear apos validar frames
contra reference.jpg).

Fail-closed (exit 2). Sem dependencias externas.
"""

import json
import re
import sys
from pathlib import Path

ESCAPE_MARKER = "# pass-l-aprovado"

KEYWORD_RE = re.compile(
    r"rugos|roughness|rough_wall|wavy_wall|topograf|"
    r"wall_amplitude|boundary_amplitude",
    re.IGNORECASE,
)


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
    if not fpath or not fpath.endswith(".py"):
        return

    text = tool_input.get("new_string") or tool_input.get("content") or ""
    if not text or ESCAPE_MARKER in text:
        return

    if not KEYWORD_RE.search(text):
        return

    rel = Path(fpath).as_posix()
    msg = (
        f"BLOQUEADO: {rel} parece introduzir geometria de contorno rugosa (Pass L).\n\n"
        "CLAUDE.md Sec.1/Sec.2.2/Sec.10: Pass L (superficies rugosas) e "
        "PRE-REQUISITO BLOQUEADO ate a simulacao reproduzir a morfologia dendritica de "
        "reference.jpg / reference_result.png painel (b) via mecanismos hidrodinamicos "
        "puros (Marangoni + Flagelar + EOS). Rugosidade e refinamento fisico, nao muleta "
        "para motor insuficiente ou selecao competitiva ausente.\n\n"
        "Se os frames ja validam a morfologia (Protocolo Sec.11 cumprido) e o usuario "
        "decidiu explicitamente desbloquear, adicione o comentario "
        f'"{ESCAPE_MARKER}: <motivo>" no trecho editado.'
    )
    print(msg, file=sys.stderr)
    sys.exit(2)


if __name__ == "__main__":
    main()
