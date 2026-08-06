#!/usr/bin/env bash
# Arquiva a rodada corrente (main_output + log.csv) em runs/<nome>/
set -euo pipefail
NAME="${1:?uso: tools/archive_run.sh <nome-da-rodada>}"
DEST="runs/${NAME}"
rm -rf "${DEST}"
mkdir -p "${DEST}"
[ -d main_output ] && cp -R main_output "${DEST}/"
[ -f log.csv ] && cp log.csv "${DEST}/"
echo "arquivado em ${DEST}"
