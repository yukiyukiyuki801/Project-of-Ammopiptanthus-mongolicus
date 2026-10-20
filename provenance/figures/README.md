# Latest saved-data figure sources

All selected drawing sources preserve their current data, geometry and aesthetics. Runtime wrapper `tools/render_saved_figure.R` evaluates the same plotting code with optional QA calls guarded. Parent must place inventory JSON at `provenance/figure_inventory.json` and restore data paths.

Optional QA: set `AMMO_FIGURE_QA_DIR` to an installed nature-figure skill directory and `AMMO_PYTHON` to a Python interpreter. No external QA source is bundled. With helpers absent, only optional QA expressions and assertions are skipped; a warning and JSON status **NOT_RUN** are emitted. No dummy exit status or synthetic PASS is supplied. Saved-result validity checks remain active. Individual original renderers can still be used with their original CLI and an installed helper directory.

Wrapper usage: `Rscript tools/render_saved_figure.R S6 runs/revision_run_01/outputs/oct16_revision/all_driver_models_postprocess revision_oct16_night/figures/output_S6`. These are drawing commands, not new analytical runs. They have not been executed during repository preparation.

The staged S14/S15 renderer has its obsolete S13 build/export branch removed. S14/S15 drawing expressions are unchanged. Latest S13 uses the separate all-driver selected-model validation renderer in inventory.
