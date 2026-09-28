#!/usr/bin/env bash
set +e
cd /home/midway/BabelSim
start_epoch=$(date +%s)
OMP_NUM_THREADS=1 mpirun -np 4 build-petsc/babelsim-solve -case cases/naca0012 > cases/naca0012/validation/formal_dt0p01_T30_mpi4_retry1/solver.log 2>&1
run_exit=$?
end_epoch=$(date +%s)
printf '{"exit_code":%d,"wall_seconds":%d,"requested_steps":3000,"deltaT":0.01,"endTime":30.0,"ranks":4,"mesh_sha256":"ef58e557a8fee0423e193bc5bcd715aec466c562b435c1cc03f07f267f54fe25"}\n' "$run_exit" "$((end_epoch - start_epoch))" > cases/naca0012/validation/formal_dt0p01_T30_mpi4_retry1/run_status.json.tmp
mv cases/naca0012/validation/formal_dt0p01_T30_mpi4_retry1/run_status.json.tmp cases/naca0012/validation/formal_dt0p01_T30_mpi4_retry1/run_status.json
exit "$run_exit"
