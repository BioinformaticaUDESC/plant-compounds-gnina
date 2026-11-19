# Complementary ADMET analysis: local execution

This Ubuntu/Linux package uses Conda. The analysis is exploratory and complementary to the manuscript's docking study. It does not establish that any substance is safe, effective, or suitable for human use.

## Running the workflow

Extract the ZIP archive and, in a terminal, run:

```bash
cd ~/Downloads/ADMET_Plantas_Local_v1
bash install_cpu.sh
conda activate admet_plantas_v1
bash run_local.sh
bash pack_results.sh
```

Adjust only the `cd` command if the archive was extracted elsewhere. The installer requires internet access and creates a separate `admet_plantas_v1` environment. If that name already exists, it stops rather than overwriting the environment. Inference uses CPU, four threads, and zero workers; it does not depend on V100 GPUs and does not modify GPU drivers or CUDA. The model package follows the official ADMET-AI distribution. No SMILES strings are sent to an external prediction server.

Each execution creates a new directory under `results/` and a log under `logs/`; prior results are not overwritten. If an error occurs, the process stops and `protocol.json` records `failed`. Partial files do not constitute a completed execution.

`pack_results.sh` packages the most recent complete inference as `ADMET_Resultados_Local_*.zip`. This is the file to provide when incorporating the analysis, supplementary tables, and manuscript comments. This module does not alter the manuscript.

To validate identities and calculate descriptors without ADMET inference:

```bash
python run_admet.py --mode descriptors --out results/descriptors_only
```

This mode does **not** generate absorption, distribution, metabolism, excretion, or toxicity predictions.

## Preserved library and input identity

The data derive from `Supplementary_Data_7_Plants_GNINA_Final.zip`, original directory `Restinga_GNINA_Alvos_V4`. The package includes the original CSV files and the 59 SDF files with their original hashes. Historical species-association files are retained for provenance and were not re-reviewed in this module; any association used in narrative text should follow the current manuscript.

The library contains 40 compounds, 51 declared structures, and 59 conformational docking inputs. The primary analytical unit is a chemical structure, not a conformer or random seed. The 59 inputs are not 59 independent molecules. The 51 structures retain their identifiers, stereochemistry, and microstates. Before analysis, the program verifies the SMILES, molecular formula, charge, InChIKey, and SHA-256 of every SDF file. Identity discrepancies stop the workflow.

The 32 volatile compounds and eight non-volatile extract comparators retain a separate `panel` field in all tables. The analysis does not demonstrate sample occurrence, properties of an essential oil, or mixture effects. No compounds were added by taxonomic inference.

Rotundic acid is represented with charge −1 in the docking library, and the primary result retains that representation. An additional sensitivity analysis applies RDKit Uncharger only where neutralisation changes the SMILES. This representation neither replaces the primary input nor models pKa, physiological pH, or microstate distributions. One additional structure is expected in the sensitivity analysis.

## Implemented methods and interpretation

**RDKit 2025.9.5** calculates average molecular weight, cLogP, TPSA, hydrogen-bond donors and acceptors, rotatable bonds, heavy atoms, Csp3 fraction, QED, Lipinski violations (MW >500, cLogP >5, HBA >10, HBD >5), and PAINS/BRENK/NIH alerts. These are structural descriptors and filters. Alerts do not demonstrate toxicity; Lipinski and QED values were not used to exclude natural products or revise docking selection.

**ADMET-AI 2.0.1** with **Chemprop 2.2.2** on CPU obtains 41 trained endpoints using `ADMETModel(include_physchem=False, drugbank_path=None, num_workers=0)`. The workflow preserves continuous outputs, verifies each SMILES alignment, and archives model-weight hashes and software versions. Core dependencies are pinned and the complete resolved environment is recorded in `environment_freeze.txt`. The ADMET environment is independent of the original docking environment.

Version 2 was retrained and does not numerically reproduce the version-1 model described in the 2024 paper or the version-1-based server. Cite the methodological publication and explicitly identify version 2.0.1.

The panel covers intestinal absorption, bioavailability, solubility, permeability, P-gp, blood–brain-barrier penetration, protein binding, volume of distribution, CYPs, half-life, clearance, and toxicity-related endpoints including hERG, AMES, and DILI. The official endpoint file retains the vendor definitions, species, units, dataset sizes, and reference metrics. Metadata and units remain on the vendor scale; the workflow does not convert outputs to doses or concentrations. Vendor metrics are not an independent performance evaluation for this panel, and individual endpoints require endpoint-specific interpretation.

Classification scores indicate the positive class used for each training set; they are not calibrated clinical probabilities. High values have different meanings across endpoints. There is no universal 0.5 threshold and no global “safe/toxic” label. Regression predictions are not truncated to physical bounds, to avoid concealing problematic values. LD50 remains on the model's logarithmic scale and is not converted to a human dose. Absence of alerts does not demonstrate safety.

No validated applicability domain or individual uncertainty estimate was established for these structures. Minimum-to-maximum intervals across declared structures describe representational spread, not confidence intervals. The workflow does not calculate a global ADMET score, combine endpoints with docking energies, or promote/exclude compounds using a threshold selected after observing results. The analysis should therefore be described as complementary and post hoc, while preserving the historical docking-selection record.

## Output files

| File | Content |
| --- | --- |
| `protocol.json` | Run status, versions, options, counts, and limitations |
| `validated_structures.csv` | 51 structures with the canonical SMILES used by the model |
| `input_identity_audit.csv` | Audit of all 59 SDF inputs |
| `docking_input_mapping.csv` | Mapping to docking identifiers |
| `physchem_by_structure.csv` | Structural descriptors and alerts |
| `admet_by_structure.csv` | 51 structures × 41 predicted endpoints |
| `admet_long_with_units.csv` | Predictions with endpoint definitions, units, and species |
| `admet_by_compound_ranges.csv` | Per-compound and per-panel ranges, without averaging stereoisomers |
| `charge_sensitivity_inputs.csv` | Neutralised representations differing from primary structures |
| `charge_sensitivity_predictions.csv` | Comparison of original and neutralised forms |
| `raw_model_predictions.csv` | Raw output for unique SMILES, without DrugBank percentiles |
| `vendor_endpoint_metadata*.csv` | Original vendor dictionary |
| `model_weights_sha256.json` | Identity of model weights actually used |
| `environment_freeze.txt` | Fully resolved execution environment |
| `input_checksums.json`, `output_sha256.json` | Input/output traceability |

Only a `protocol.json` file marked `completed_full_inference` supports the statement that ADMET inference was completed. Descriptor-only mode records a distinct status. Validation outputs distributed with the package are identified in `VALIDACAO.md` and do not substitute for a local execution.

## References

- Swanson K. et al. *ADMET-AI: a machine learning ADMET platform for evaluation of large-scale chemical libraries.* Bioinformatics. 2024;40:btae416. https://doi.org/10.1093/bioinformatics/btae416
- Implementation used: https://github.com/swansonk14/admet_ai/tree/v_2.0.1
- Version-2 changes: https://github.com/swansonk14/admet_ai/releases/tag/v_2.0.0
- RDKit: https://www.rdkit.org/docs/

ADMET-AI metadata are reproduced under their MIT licence in `reference/ADMET_AI_LICENSE.txt`. The scripts in this package are new workflow scripts, not historical docking scripts.
