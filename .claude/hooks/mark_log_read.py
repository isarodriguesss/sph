#!/usr/bin/env python3
"""Hook PostToolUse (Read): marca no estado da sessao que log.csv foi lido.

Consumido por guard_param_change.py, que exige uma leitura FRESCA de log.csv
antes de cada mudanca de parametro fisico (CLAUDE.md Sec.10, proibicao
"Nao alterar parametros fisicos sem antes ler log.csv"; protocolo Sec.2.1).

Sem dependencias externas.
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

STATE_DIR = Path(__file__).resolve().parent.parent / ".state"
STATE_FILE = STATE_DIR / "log_read.json"


def main():
    try:
        data = json.load(sys.stdin)
    except Exception:
        return

    tool_input = data.get("tool_input", {}) or {}
    fpath = tool_input.get("file_path") or tool_input.get("path") or ""
    if not fpath.replace("\\", "/").endswith("log.csv"):
        return

    session_id = data.get("session_id", "unknown")

    STATE_DIR.mkdir(parents=True, exist_ok=True)
    state = {}
    if STATE_FILE.exists():
        try:
            state = json.loads(STATE_FILE.read_text(encoding="utf-8"))
        except Exception:
            state = {}

    state[session_id] = {
        "read_at": datetime.now(timezone.utc).isoformat(),
        "consumed": False,
    }
    STATE_FILE.write_text(json.dumps(state, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
