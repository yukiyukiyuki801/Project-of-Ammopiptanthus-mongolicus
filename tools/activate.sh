# Source from the repository root after installing the documented dependencies.
export AMMO_REVISION_RUN="$PWD/runs/revision_run_01"
export AMMO_DIR="$PWD/ammo_analysis_package_20261014/ammo"
export AMMO_ANALYSIS_OUT="$AMMO_REVISION_RUN/outputs"
export AMMO_PYTHON="${AMMO_PYTHON:-python3}"
export PYTHONPATH="$PWD/public_dependencies/software/pyRona${PYTHONPATH:+:$PYTHONPATH}"
