# Complete data and source archive for the seven-species GNINA project

**Archive delivery date:** 30 October 2026.

The recovered source materials were consolidated without modifying coordinates, searches, scores, or conclusions from the deposited control analyses. The primary library contains 40 named compounds, 51 declared stereochemical structures, and 59 primary SDF inputs. Direct counts identified 24,780 primary candidate poses and 1,200 native-control poses. Diagnostic, historical, and quarantined files remain explicitly identified and must not be added to the primary analysis.

## Requested files and locations

| Content | Path relative to this directory |
| --- | --- |
| Structure manifest | `dados/Restinga_GNINA_Alvos_V4/structures.csv` |
| Primary docking-input manifest | `dados/Restinga_GNINA_Alvos_V4/primary_docking_inputs.csv` |
| 59 primary SDF inputs | Exact paths are given in the manifest `input_path` column; files are under `dados/Restinga_GNINA_Alvos_V4/prepared/ligands`. |
| Configurations and protocols | `dados/Restinga_GNINA_Alvos_V4/protocol.json`, `protocol_round_b.json`, `dados/Restinga_GNINA_Multialvo_V3/protocol.json`, and `targets/*/config.json`. |
| Primary V3 outputs | `dados/Restinga_GNINA_Multialvo_V3/results/COX2_5IKR`, `ACHE_4EY7`, and `CB2_5ZTY/docking`; restrict to primary entries using the manifest. |
| Primary V4 outputs | `dados/Restinga_GNINA_Alvos_V4/results/SEH_3ANS`, `KEAP1_4XMB`, `PARP1_7KK4`, and `TYRS_1JIL/docking`. |
| Native-ligand controls | `results/*/redocking` within each of the two data directories. |
| Editable PyMOL sessions | `fontes_figuras/Figures_7_Plants_GNINA_V3/sessions`; 12 native controls and 7 candidate complexes. |
| Sources for Figures 3, 4, S1, and S2 | `render_pymol.py` and `compose_pymol_figures.py` in `Figures_7_Plants_GNINA_V3`. |
| Sources for Figures 5 and S3–S7 | `build_molecule_figures.py`, `figure_manifest.json`, and `molecules/*.svg` in the same directory; 51 SVG and 51 PNG files. |
| Figure 2 | `fontes_figuras/Figura_2_reconstituida`; see the provenance statement below. |
| Supplementary Data S8 and S9 | `informacoes_suplementares`. |
| Earlier verification records | `auditoria_verificada` and the original figure `verification_2026-10-27` files. |
| Methodological responses | `documentacao/Respostas_metodologicas.md`. |
| PDB-level rationale | `documentacao/Selecao_documentada_12_alvos.csv` and the original selection records. |

A consolidated SDF batch, additional historical outputs, and diagnostic copies may be present in addition to the 59 primary inputs. The inclusion unit is each row of `primary_docking_inputs.csv`, not the unfiltered count of all `*.sdf` files. The compressed V3 delivery is retained under `previous_delivery` as a provenance record and was separately expanded in this archive for direct access.

## Figure 2 provenance

The historical Figure 2 plotting script was not located in the recovered archived packages. `rebuild_figure2.py` is a later reconstruction dated 30 October 2026 and is explicitly identified as such in the code and metadata. It reads the two-decimal medians reported in Table S1, calculates protein-wise ranks with average ties, and compares all 224 values with the previously verified matrix. The 224 values matched. This verification does not assert pixel-identical reproduction of the former figure's typography, dimensions, or rendering.

The manuscript PNG is retained as `Figure_2_Original_From_Manuscript.png`, together with source-file identification and its hash. New editable SVG/PDF outputs and the reconstructed PNG are stored in `reconstructed`. All remaining figure sources and the 19 PyMOL sessions are unaltered copies of previously delivered files.

## Use of the deposited sources

A local extraction may be stored at `/home/juan/Supplementary_Data_7_Plants_GNINA_Completo`. Figure scripts use paths relative to their own directory. Run `rebuild_figure2.py` with the dependencies in `requirements-figura2.txt`. Historical PyMOL and RDKit versions are recorded in `requirements-visualization.txt` and `visualization_environment.json` in the figure directory. The reconstruction environment does not replace the historical GNINA execution environment.

For rank verification, the reconstruction script accepts `--input`, `--expected`, and `--output-dir`. To inspect a session, open the corresponding `.pse` file in PyMOL; it preserves both objects and camera settings. To reproduce molecular figures, execute `build_molecule_figures.py` in the recorded environment. `render_pymol.py` regenerates views from coordinates, and `compose_pymol_figures.py` assembles the panels. These scripts may overwrite generated outputs and manifests; they should therefore be run in a working copy while this deposited archive remains unmodified.

Historical docking commands retain their original absolute paths as execution records. They are not portable commands ready to execute from the current directory. The GNINA executable is not redistributed; its documented identity and historical dependencies are included in the archive. No new docking searches, molecular-dynamics simulations, or biological assays were performed while preparing this delivery.

## Delivery checks

`VERIFICACAO_COMPLETUDE.json` records the preserved original copies, verification of the 59 input hashes, direct pose and session counts, and rank reconstruction. `COX2_verificacao_regra_da_caixa.json` records that the rule `max(22 Å, native extent + 12 Å)` reproduces the inherited 22 × 22 × 22 Å box, including its centre. `verificar_caixa_cox2.py` reproduces only this geometric calculation using the Python standard library.

`SHA256SUMS.txt` covers every file in this package except the checksum file itself. Historical checksums and execution records are retained without modification in their original folders. The IC documents under `documentacao/comparacao_IC_Docking_OE` are examples used for comparing score ordering and aggregation; they are not part of the seven-plant screening analysis and do not constitute a complete audit of that separate study.
