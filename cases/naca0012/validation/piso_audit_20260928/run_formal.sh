#!/usr/bin/env bash
cd /home/midway/BabelSim || exit 1
run_dir=cases/naca0012/validation/piso_audit_20260928/formal
mkdir -p "$run_dir"
start_epoch=$(date +%s)
OMP_NUM_THREADS=1 mpirun -np 4 build-petsc/babelsim-solve -case cases/naca0012 > "$run_dir/solver.log" 2>&1
run_exit=$?
end_epoch=$(date +%s)
printf '{"exit_code":%d,"wall_seconds":%d,"deltaT":0.01,"endTime":30,"ranks":4}\n' "$run_exit" "$((end_epoch-start_epoch))" > "$run_dir/status.json"
exit "$run_exit"
