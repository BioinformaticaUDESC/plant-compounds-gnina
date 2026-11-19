# GNINA multi-target V3: reproducibility record

This CPU-based exploratory study was run on 25 October 2026. The accompanying PDF and `reports/RESULTADOS_V3.md` provide the technical summary. Project-specific GC–MS data were not available; relationships among compounds, species, and matrices are literature-derived.

## Scope and native controls

The library contains 40 literature-derived compounds (32 volatile constituents and 8 extract comparators), represented by 51 declared stereochemical structures. Four systems were assessed by native-ligand redocking with five seeds. The prespecified criterion required a symmetry-aware heavy-atom RMSD ≤2 Å for the first Energy-ranked pose, in the fixed protein reference frame, for at least 4/5 seeds. The starting reference ligand was prepared independently of the crystallographic conformation.

- Human COX-2, 5IKR, mefenamate/ID8: **pass**, 5/5; median RMSD 0.481 Å.
- Human acetylcholinesterase, 4EY7 chain B, protonated donepezil/E20: **pass**, 5/5; median RMSD 1.070 Å.
- Human CB2/T4-lysozyme crystallographic construct, 5ZTY, AM10257/9JU: **pass**, 5/5; median RMSD 0.537 Å. The antagonist-bound construct does not establish agonism for any candidate.
- Human 5-LOX stable construct, 6N2W chain B, NDGA/30Z and Fe²⁺: **fail**, 0/5; median RMSD 3.549 Å. Crystal-like poses were sampled but were not first by Energy; no candidate was screened against this target under the protocol.

For accepted targets, the library was screened with seeds 42, 142, and 242, with up to 20 poses per input. These repetitions vary the search, not the rigid protein structure. Scores from different proteins are not a common affinity scale and do not establish efficacy, mixture activity, selectivity, toxicity, or safety.

## Curated primary library

The initial 60-input set contained `caryophyllene_oxide_conf03`, whose three-dimensional stereochemical identity differed from the declared structure. Its files and outputs are retained in quarantine. A corrected, identity-verified geometry was tested in the accepted systems but remained a separate diagnostic because it was not conformationally independent (0.5 Å threshold) from a retained correct conformer.

The final primary library comprises **59 verified conformational inputs**, retaining all 40 compounds and 51 declared structures. Primary synthesis and contact tables use this library only.

- `docking_inputs.csv` and `prepared/all_inputs.sdf`: historical 60-input batch, including the quarantined input.
- `primary_docking_inputs.csv` and `prepared/primary_inputs.sdf`: curated 59-input primary library.
- `all_attempted_docking_inputs.csv`: all 61 attempted inputs, including the historical discrepancy and the correction diagnostic.
- `correction_docking_inputs.csv`: corrected diagnostic outside the primary set.
- `reports/input_identity_qc.csv` and `reports/input_curation.json`: expected/observed identities, exclusion rationale, and correction selection.

## Protocol and quality checks

GNINA 1.3.3 was fingerprinted by SHA-256 and executed on CPU (four CPUs per process; at most two simultaneous GNINA processes) using `--cnn_scoring rescore --cnn fast --pose_sort_order Energy --exhaustiveness 16 --num_modes 20 --no_gpu`. Energy ranking was specified before this analysis. CNN-based reordering within exported poses is a secondary diagnostic, not an independent search.

Deposited receptor heavy-atom coordinates were retained; hydrogens were modelled at pH 7.4. Loops were not rebuilt and no heavy-atom relaxation was performed. Waters, glycans, and additives were excluded. The structural limitations of the COX-2, 5-LOX, and CB2 constructs are retained in `target_rationale.json`.

The audit checks chemical identity, finite coordinates, non-metal contacts <2 Å, distance/radius ratios <0.70, positive energy, box escape >0.1 Å, mapped bond-length changes >0.15 Å, and severe intramolecular overlaps. Metal contacts are reported separately. An absence of flags does not establish pose correctness, and residue proximity ≤4 Å is not automatically assigned as a hydrogen bond or mechanism.

## Contents and reproduction

- `reports/`: technical PDF, Markdown synthesis, compound/structure summaries, ranking comparisons, recurrent contacts, and PNG/SVG figures.
- `results/<target>/redocking/`: all native-control runs, including failures.
- `results/<target>/docking/`: historical and correction-diagnostic runs, including commands, logs, and completion hashes.
- `results/<target>/primary_docking_poses.csv`: poses from the curated primary library only.
- `targets/<target>/view_redocking.pml`: PyMOL inspection in the unaligned reference frame.
- `protocol.json`, `correction_protocol.json`, `method_checks.json`, and `native_identity_checks.json`: parameter and validation records.
- `SHA256SUMS.txt`: delivery inventory.

To prepare the CPU environment:

```bash
cd "$HOME/Restinga_GNINA_Multialvo_V3"
INSTALL_REPORT=1 bash 00_preparar_ambiente.sh
```

To regenerate delivered tables and figures only:

```bash
bash run.sh audit --reuse-audited
bash run.sh report
```

For a clean rerun, generate a new directory containing only the 59 primary inputs:

```bash
.venv/bin/python prepare_clean_reproduction.py \
  --destination "$HOME/Restinga_GNINA_V3_Reproducao_Limpa"
```

The tool refuses an existing destination and does not automatically start simulations. The clean reproduction excludes quarantined and diagnostic inputs. Historical-batch reproduction requires its explicit option because the current preflight rejects the quarantined input.

Before prioritising compounds, integrate GC–MS data, retention indices, and matrix/species-specific abundances; confirm vouchers; resolve missing composition data; and test hypotheses in target-appropriate assays. Docking alone cannot identify the most potent or safest plant.
