#!/usr/bin/env python3
"""Hook PreToolUse (Bash): lembra de pedir confirmacao antes de rodar
simulacoes SPH longas.

Nao bloqueia (a duracao real depende de total_sim_time/resolucao, que o hook
nao pode calcular com seguranca) -- so injeta um lembrete de contexto quando o
comando parece disparar uma simulacao completa (make run / python main.py).

Reforca estruturalmente o que ja e licao registrada em memoria (never
auto-background long SPH runs; ask first) para o caso de a memoria ter saido
do contexto por compactacao.

Sem dependencias externas.
"""

import json
import re
import sys

RUN_RE = re.compile(r"\bmake\s+run(_view)?\b|\bpython[3]?\s+.*main\.py\b")


def main():
    try:
        data = json.load(sys.stdin)
    except Exception:
        return

    if data.get("tool_name") != "Bash":
        return

    command = (data.get("tool_input", {}) or {}).get("command", "") or ""
    if not RUN_RE.search(command):
        return

    ctx = (
        "Lembrete: este comando parece disparar uma simulacao SPH completa "
        "(main.py). Simulacoes podem levar minutos a horas dependendo de "
        "total_sim_time/resolucao/Pass N. Confirme com a usuaria antes de "
        "rodar, e nunca inicie em background silenciosamente -- se for "
        "genuinamente longa, use run_in_background so apos combinar isso."
    )
    print(
        json.dumps(
            {
                "hookSpecificOutput": {
                    "hookEventName": "PreToolUse",
                    "additionalContext": ctx,
                }
            }
        )
    )


if __name__ == "__main__":
    main()
