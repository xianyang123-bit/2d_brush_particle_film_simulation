#!/usr/bin/env bash
# Uses the HOOMD environment already installed by install_and_run_film.sh.
set -euo pipefail
task_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python_bin="${XDG_DATA_HOME:-$HOME/.local/share}/gnp-hoomd/env/bin/python"
if [[ ! -x "$python_bin" ]]; then
    echo "Installed HOOMD interpreter not found: $python_bin"
    exit 1
fi
if [[ "${1:-}" != --worker ]]; then
    "$python_bin" -c 'import hoomd, gsd.hoomd, matplotlib; print("HOOMD", hoomd.version.version)'
    run_dir="$(mktemp -d "$task_dir/film400_long_XXXXXXXX")"
    cp -- "$task_dir/gnp_film.py" "$task_dir/analyze_film.py" "$task_dir/run_film400_long.sh" "$run_dir/"
    printf '%s\n' "$run_dir" > "$task_dir/latest_film400_run.txt"
    printf 'STARTING\n' > "$run_dir/status.txt"
    nohup bash "$run_dir/run_film400_long.sh" --worker "$run_dir" > "$run_dir/run.log" 2>&1 < /dev/null &
    printf '%s\n' "$!" > "$run_dir/pid.txt"
    printf 'Started background run, PID %s\nResults: %s\nLog: %s/run.log\n' "$!" "$run_dir" "$run_dir"
    exit 0
fi
run_dir="$2"
trap 'printf "FAILED at line %s\n" "$LINENO" > "$run_dir/status.txt"' ERR
printf 'RUNNING\n' > "$run_dir/status.txt"
"$python_bin" -u "$run_dir/gnp_film.py" \
    --nside 20 --gap 4 --kT 1.0 --dt 0.0005 \
    --tau-thermostat 0.2 --tau-barostat 50 \
    --relax-steps 100000 \
    --prep-pressure 0.1 --prep-steps 500000 \
    --pressure 0 --npt-steps 1000000 \
    --equil-steps 1000000 --production-steps 1000000 \
    --period 10000 --checkpoint-period 100000 \
    --out "$run_dir/film400"
printf 'ANALYZING\n' > "$run_dir/status.txt"
"$python_bin" "$run_dir/analyze_film.py" "$run_dir/film400"
printf 'COMPLETE\n' > "$run_dir/status.txt"
