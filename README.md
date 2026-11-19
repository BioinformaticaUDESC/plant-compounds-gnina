# Plant compounds: GNINA docking, redocking, and ADMET

## Reproducibility archive for peer review

This repository provides the computational records underlying the revised manuscript, *Compound–Target Hypotheses for the Reported Bioactivities of Six Plant Species Occurring in Brazil: A Redocking-Gated, Size-Corrected In Silico Screen of 40 Constituents*.

It is designed to support transparent editorial and peer-review assessment. The archive preserves the computational inputs, programmatic workflows, execution records, output structures, audit files, and figure sources used in the study. It should be read as a reproducibility record—not as evidence of experimentally confirmed molecular binding, biological activity, selectivity, pharmacokinetics, or safety.

| Study component | Archived record |
| --- | ---: |
| Plant species addressed in the manuscript | 6 |
| Literature-derived compounds | 40 (32 volatile constituents and 8 non-volatile extract comparators) |
| Declared chemical structures / primary SDF inputs | 51 / 59 |
| Systems tested by native-ligand redocking / eligible for screening | 12 / 7 |
| Primary candidate poses / native-control poses | 24,780 / 1,200 |
| ADMET predictions | 2,091 (41 endpoints × 51 structures) |

## Computational workflow and deposited records

| Stage | Primary records |
| --- | --- |
| **1. Literature curation and chemical library definition** | [Species–source mapping](Supplementary_Data_7_Plants_GNINA_Completo/dados/Restinga_GNINA_Alvos_V4/literature_species_links.csv), [structure manifest](Supplementary_Data_7_Plants_GNINA_Completo/dados/Restinga_GNINA_Alvos_V4/structures.csv), [primary docking-input manifest](Supplementary_Data_7_Plants_GNINA_Completo/dados/Restinga_GNINA_Alvos_V4/primary_docking_inputs.csv), and [V3 workflow record](Supplementary_Data_7_Plants_GNINA_Completo/dados/Restinga_GNINA_Multialvo_V3/README.md). |
| **2. Native-ligand redocking and candidate docking** | [V3 protocol](Supplementary_Data_7_Plants_GNINA_Completo/dados/Restinga_GNINA_Multialvo_V3/protocol.json), [V4 protocol](Supplementary_Data_7_Plants_GNINA_Completo/dados/Restinga_GNINA_Alvos_V4/protocol.json), [V4 round-B protocol](Supplementary_Data_7_Plants_GNINA_Completo/dados/Restinga_GNINA_Alvos_V4/protocol_round_b.json), and [V4 workflow record](Supplementary_Data_7_Plants_GNINA_Completo/dados/Restinga_GNINA_Alvos_V4/README.md). Pose recovery was assessed using the highest-ranked Energy pose and fixed-frame heavy-atom RMSD. Five systems did not satisfy the prespecified criterion and were not advanced to candidate screening. |
| **3. Audit, score adjustment, and sensitivity analysis** | [Verified audit records](Supplementary_Data_7_Plants_GNINA_Completo/auditoria_verificada/), [archived commands](Supplementary_Data_7_Plants_GNINA_Completo/informacoes_suplementares/Supplementary_Data_S9_Archived_Commands.csv), and [revision sensitivity analyses](sensitivity-r2/). |
| **4. ADMET analysis** | [Local execution guide](admet/ADMET_Plantas_Local_v1/README.md), [CPU reference run](admet/ADMET_Plantas_Local_v1/validation_reference/), [reference protocol](admet/ADMET_Plantas_Local_v1/validation_reference/protocol.json), and [frozen environment](admet/ADMET_Plantas_Local_v1/validation_reference/environment_freeze.txt). |
| **5. Figures and coordinate-level evidence** | [Historical figure scripts and source files](Supplementary_Data_7_Plants_GNINA_Completo/fontes_figuras/), [revision annotation scripts](figures-r2/), and [minimum residue-distance table](figures-r2/residue_proximities.csv). The archive includes rendered PNG/SVG outputs and 19 original PyMOL sessions. |

## Accessing and verifying the complete archive

For a complete local copy, clone the repository or use **Code → Download ZIP**. To restore and verify the deposited archive from the repository root:

```sh
python3 archive/restore_archive.py --output Plant_Compounds_GNINA_GitHub_Deposit_R2.zip
unzip Plant_Compounds_GNINA_GitHub_Deposit_R2.zip
(cd plant-compounds-gnina && sha256sum -c --quiet SHA256SUMS)
```

The restoration script validates the archive parts and the reconstructed ZIP against the [SHA-256 manifest](archive/manifest.json). The deposited ZIP contains 1,403 files, including primary SDF inputs, prepared receptors, docking poses, execution logs, frozen environments, analysis scripts, and the PyMOL sessions. The text files mirrored in the repository facilitate browsing; commands intended for reproduction should be executed from the restored archive, which contains the referenced molecular inputs. The GNINA binary and model weights are identified by recorded fingerprints and are not redistributed.

## Scope and interpretation

Native-ligand redocking evaluates pose recovery for each deposited cognate ligand under the recorded protocol. The adjusted rankings, the eleven reported compound–target pairs, and the ADMET outputs are exploratory prioritization results. They do not constitute independent validation of prospective hit selection, molecular binding, functional activity, selectivity, exposure, or safety.

The manuscript evaluates six species with attributable constituent data. The historical archive retains an earlier seven-species project label in some paths and filenames; this nomenclature is retained solely to preserve provenance. The Figure 2 plotting script was not available in the recovered historical materials; its later reconstruction is clearly labelled and documented in the source archive.

## References and archive guide

[Complete bibliography](REFERENCES.md) · [Earlier reference audit](Supplementary_Data_7_Plants_GNINA_Completo/auditoria_verificada/reference_audit.csv) · [Detailed source-archive guide](Supplementary_Data_7_Plants_GNINA_Completo/README.md)
