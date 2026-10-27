#!/bin/sh
# ==============================================================================
# Runs the whole evaluation study unattended (from the package root):
#   1. waits until --start (HH:MM, local time; default: start now)
#   2. the accuracy tasks of sim_accuracy.R (beta and sfh) on $WORKERS processes
#   3. collects them into inst/extdata/eval_{sfh,beta}_accuracy.rds
#   4. bench_resources.R, one cell at a time on the then idle machine
# caffeinate keeps the Mac awake. Everything resumes: rerunning the script skips
# finished replicates and cells. Logs go to benchmarks/output/eval/logs/.
#
# Usage: nohup sh benchmarks/eval/run_overnight.sh [HH:MM] > /dev/null 2>&1 &
# Run R CMD INSTALL . --preclean first: the scripts use the installed fastsae.
# ==============================================================================

START=${1:-}
WORKERS=${WORKERS:-6}
LOG=benchmarks/output/eval/logs
mkdir -p "$LOG"

if [ -n "$START" ]; then
  target=$(date -j -f "%Y-%m-%d %H:%M" "$(date +%Y-%m-%d) $START" +%s)
  [ "$target" -le "$(date +%s)" ] && target=$((target + 86400))
  echo "$(date '+%F %T') waiting until $(date -r "$target" '+%F %T')" >> "$LOG/progress.log"
  while [ "$(date +%s)" -lt "$target" ]; do sleep 30; done
fi

exec caffeinate -i -s sh -c '
  LOG='"$LOG"'
  echo "$(date "+%F %T") accuracy: start ('"$WORKERS"' workers)" >> $LOG/progress.log
  { Rscript benchmarks/eval/sim_accuracy.R --list=beta; Rscript benchmarks/eval/sim_accuracy.R --list=sfh; } > $LOG/tasks.txt
  xargs -P '"$WORKERS"' -I{} sh -c "Rscript benchmarks/eval/sim_accuracy.R --task={} >> $LOG/accuracy_tasks.log 2>&1 || echo \"FAILED {}\" >> $LOG/accuracy_tasks.log" < $LOG/tasks.txt
  echo "$(date "+%F %T") accuracy: done" >> $LOG/progress.log
  Rscript benchmarks/eval/sim_accuracy.R --collect >> $LOG/progress.log 2>&1
  echo "$(date "+%F %T") resources: start" >> $LOG/progress.log
  Rscript benchmarks/eval/bench_resources.R > $LOG/resources.log 2>&1
  echo "$(date "+%F %T") resources: done (exit $?)" >> $LOG/progress.log
'
