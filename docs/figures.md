# Figure and table reconstruction

The checked-in `figures/Main` and `figures/Supplement` contain the exact latest six main and twenty supplementary figures, each as PDF and PNG. Their SHA256 values and original render paths are in `provenance/figure_inventory.json`. `tables/Supplementary_tables.xlsx` contains the current28 individually numbered tables. The inventory preserves their numbering/source associations.

Restore data first with `python3 tools/prepare_data.py`, then source `tools/activate.sh`. Use the wrapper below so missing optional local figure-QA tools are marked NOT_RUN rather than required for drawing or falsely marked PASS. Numerical-result validation guards remain active. The wrapper does not change the figure's data, geometry or palette.

```bash
Rscript tools/render_saved_figure.R 6
Rscript tools/render_saved_figure.R S6 runs/revision_run_01/outputs/oct16_revision/all_driver_models_postprocess revision_oct16_night/figures/output_S6
Rscript tools/render_saved_figure.R S13 runs/revision_run_01/outputs/oct16_revision/all_driver_models_postprocess runs/revision_run_01/outputs/oct16_revision/all_driver_models_postprocess/figure_S13
Rscript tools/render_saved_figure.R S20 data/mapping_plot_config.json
```

Other current renderers take no arguments unless listed below. A shared renderer may draw more than one paired figure. Output files retain the executed source names; use the source-stem column to identify the requested figure. Published filenames under `figures/` are supplied separately, so running a shared renderer does not silently replace the released final figures with an older intermediate panel.

| Released figure | Drawing source | Argument notes | Executed output stem |
|---|---|---|---|
| 1 | `scripts/revision/remake_main_1_4_original.R` | Default/no arguments | `runs/revision_run_01/outputs/figure_revision_original_design/main_1_4/Figure_1` |
| 2 | `scripts/revision/remake_main_2_6_original.R` | Default/no arguments | `runs/revision_run_01/outputs/figure_revision_original_design/main_2_6/Figure2_original_design` |
| 3 | `scripts/revision/remake_main_3_5_original.R` | Default/no arguments | `runs/revision_run_01/outputs/figure_revision_original_design/main_3_5/Figure_3_original_design` |
| 4 | `scripts/revision/remake_main_1_4_original.R` | Default/no arguments | `runs/revision_run_01/outputs/figure_revision_original_design/main_1_4/Figure_4` |
| 5 | `revision_oct18/figure_updates/render_main5.R` | Default/no arguments | `revision_oct18/figure_updates/main5/Figure_5_original_design` |
| 6 | `revision_oct18_new_structure/figure6_opacity/render_Figure6.R` | Default/no arguments | `revision_oct18_new_structure/figure6_opacity/Figure_6` |
| S1 | `revision_oct18_discussion/figure_niche/render_S1_asia.R` | Default/no arguments | `revision_oct18_discussion/figure_niche/S1/S1_original_design` |
| S2 | `revision_oct18/figure_updates/render_S2.R` | Default/no arguments | `revision_oct18/figure_updates/S2/S2_original_design` |
| S3 | `scripts/revision/remake_S1_S4_original.R` | Default/no arguments | `runs/revision_run_01/outputs/figure_revision_original_design/supplement_S1_S4/S3_original_design` |
| S4 | `scripts/revision/remake_S1_S4_original.R` | Default/no arguments | `runs/revision_run_01/outputs/figure_revision_original_design/supplement_S1_S4/S4_original_design` |
| S5 | `scripts/revision/remake_S5_original.R` | Default/no arguments | `runs/revision_run_01/outputs/figure_revision_original_design/supplement_S5/Figure_S5` |
| S6 | `revision_oct16_night/figures/render_S6_combined.R` | Two arguments: postprocess directory and output directory | `revision_oct16_night/figures/output_S6/Figure_S6` |
| S7 | `scripts/revision/remake_S7_S8_original.R` | Default/no arguments | `runs/revision_run_01/outputs/figure_revision_original_design/supplement_S7_S8/Figure_S7` |
| S8 | `scripts/revision/remake_S7_S8_original.R` | Default/no arguments | `runs/revision_run_01/outputs/figure_revision_original_design/supplement_S7_S8/Figure_S8` |
| S9 | `scripts/revision/remake_S9_S10_original.R` | Default/no arguments | `runs/revision_run_01/outputs/figure_revision_original_design/supplement_S9_S10/Figure_S9` |
| S10 | `scripts/revision/remake_S9_S10_original.R` | Default/no arguments | `runs/revision_run_01/outputs/figure_revision_original_design/supplement_S9_S10/Figure_S10` |
| S11 | `scripts/revision/remake_S11_S12_original.R` | Default/no arguments | `runs/revision_run_01/outputs/figure_revision_original_design/supplement_S11_S12/Figure_S11` |
| S12 | `scripts/revision/remake_S11_S12_original.R` | Default/no arguments | `runs/revision_run_01/outputs/figure_revision_original_design/supplement_S11_S12/Figure_S12` |
| S13 | `revision_oct16/driver_combinations/figures/render_driver_validation.R` | Two arguments: postprocess directory and output directory | `runs/revision_run_01/outputs/oct16_revision/all_driver_models_postprocess/figure_S13/Figure_S13` |
| S14 | `scripts/revision/render_requested_S13_S15.R` | Default/no arguments | `runs/revision_run_01/outputs/figure_revision_original_design/supplement_S13_S15/Figure_S14` |
| S15 | `scripts/revision/render_requested_S13_S15.R` | Default/no arguments | `runs/revision_run_01/outputs/figure_revision_original_design/supplement_S13_S15/Figure_S15` |
| S16 | `revision_oct18/figure_updates/render_supp_pies.R` | Default/no arguments | `revision_oct18/figure_updates/distance_pies/Figure_S16` |
| S17 | `revision_oct18/figure_updates/render_supp_pies.R` | Default/no arguments | `revision_oct18/figure_updates/distance_pies/Figure_S17` |
| S18 | `revision_oct18/figure_updates/render_supp_pies.R` | Default/no arguments | `revision_oct18/figure_updates/distance_pies/Figure_S18` |
| S19 | `revision_oct18/figure_updates/render_supp_pies.R` | Default/no arguments | `revision_oct18/figure_updates/distance_pies/Figure_S19` |
| S20 | `revision_oct18/figure_updates/render_S20.R` | data/mapping_plot_config.json | `revision_oct18/figure_updates/S20/Figure_S20_reference_PCA` |

Figure6 uses the final approved green opacity0.50. The S14/S15 renderer excludes its obsolete S13 branch; the final S13 comes from the all-driver nested-validation source. Optional plot-object files that are not needed by the drawing code are distinguished from required inputs in the inventory.

Optional QA is enabled by setting `AMMO_FIGURE_QA_DIR` to a separately installed nature-figure helper directory. `AMMO_PYTHON` selects the interpreter. Without those helpers, the wrapper warns and records `NOT_RUN`; it skips only helper-dependent checks and their metadata, never substitutes a successful return status. Release preparation checked code parsing, source hashes and the optional-helper control path but did not rerender every figure in a fresh environment.
