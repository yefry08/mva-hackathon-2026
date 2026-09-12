#!/usr/bin/env bash
# Driver de la corrida de 24h. Reanudable, con checkpoints y commits periodicos.
#   tmux new -s mva
#   ./orchestrator/run.sh
set -uo pipefail
cd "$(dirname "$0")/.."

: "${MAX_CONCURRENT:=15}"
: "${CONSOLIDATE_EVERY:=7200}"   # 2h
: "${COMMIT_EVERY:=1800}"        # 30min
: "${RUN_HOURS:=24}"

export MAX_CONCURRENT
START=$(date +%s)
DEADLINE=$((START + RUN_HOURS * 3600))
LAST_COMMIT=$START
LAST_CONSOLIDATE=$START

mkdir -p reports/workers

guard() {
  # Aborta si hay datos de paciente staged. No negociable.
  if git status --porcelain 2>/dev/null | grep -Eq '(^A|^M).*(data/|ref/|\.fastq|\.bam|\.vcf|\.cram)'; then
    echo "ABORT: datos de paciente en staging. Revisa .gitignore." >&2
    git reset >/dev/null 2>&1
    return 1
  fi
  return 0
}

echo "Inicio: $(date -u +%FT%TZ) | deadline: $(date -u -d "@$DEADLINE" +%FT%TZ 2>/dev/null || echo "+${RUN_HOURS}h")"

while [ "$(date +%s)" -lt "$DEADLINE" ]; do
  python3 orchestrator/dispatch.py --max-concurrent "$MAX_CONCURRENT"
  RC=$?
  NOW=$(date +%s)

  if [ $((NOW - LAST_CONSOLIDATE)) -ge "$CONSOLIDATE_EVERY" ]; then
    claude -p "Agente de convergencia. Lee reports/workers/*.json y orchestrator/state.json.
Consolida hallazgos, elimina duplicados, detecta contradicciones entre workers, y reescribe
reports/STATUS.md con: mejores candidatos y su confianza, evidencia ortogonal que los apoya,
bloqueos activos, y presupuesto de submissions restante. Si algun hallazgo cambia las
prioridades, agrega tareas nuevas al final de orchestrator/queue.jsonl con IDs no usados.
No envies nada ni hagas push." --permission-mode acceptEdits
    LAST_CONSOLIDATE=$NOW
  fi

  if [ $((NOW - LAST_COMMIT)) -ge "$COMMIT_EVERY" ]; then
    if guard; then
      git add -A ':!data' ':!ref' 2>/dev/null
      git commit -q -m "checkpoint $(date -u +%FT%TZ)" 2>/dev/null || true
    fi
    LAST_COMMIT=$NOW
  fi

  if [ $RC -eq 0 ]; then
    echo "Cola vacia. Esperando tareas nuevas del agente de convergencia (5 min)..."
    sleep 300
  else
    sleep 60
  fi
done

echo "Fin de ventana. Estado final en reports/STATUS.md"
