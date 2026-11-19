# Plant compounds · GNINA docking and redocking

**Reproducibility repository** for *Compound-Target Hypotheses for the Reported Bioactivities of Six Plant Species Occurring in Brazil: A Redocking-Gated, Size-Corrected In Silico Screen of 40 Constituents*

This project used GNINA to generate hypotheses about molecular targets of literature-reported plant constituents. It preserves the computational inputs, outputs, analysis code, audit records and figure sources used in the revised study.

| Study scope | Archived record |
| --- | ---: |
| Plant species in the manuscript | 6 |
| Compounds | 40 (32 volatiles and 8 extract comparators) |
| Declared chemical structures / SDF inputs | 51 / 59 |
| Systems tested by native redocking / eligible for screening | 12 / 7 |
| Primary candidate poses / native control poses | 24,780 / 1,200 |
| ADMET predictions | 2,091 (41 endpoints × 51 structures) |

## What was run

| Step | Code and execution records |
| --- | --- |
| **1 · Literature and library curation** | [Species–source mapping](Supplementary_Data_7_Plants_GNINA_Completo/dados/Restinga_GNINA_Alvos_V4/literature_species_links.csv), [structure manifest](Supplementary_Data_7_Plants_GNINA_Completo/dados/Restinga_GNINA_Alvos_V4/structures.csv), [59 primary docking inputs](Supplementary_Data_7_Plants_GNINA_Completo/dados/Restinga_GNINA_Alvos_V4/primary_docking_inputs.csv), [V3 workflow](Supplementary_Data_7_Plants_GNINA_Completo/dados/Restinga_GNINA_Multialvo_V3/README.md) |
| **2 · Native redocking and candidate screening** | [V3 protocol](Supplementary_Data_7_Plants_GNINA_Completo/dados/Restinga_GNINA_Multialvo_V3/protocol.json), [V4 protocol](protocol.json), [V4 round B](Supplementary_Data_7_Plants_GNINA_Completo/dados/Restinga_GNINA_Alvos_V4/protocol_round_b.json), [V4 workflow](Supplementary_Data_7_Plants_GNINA_Completo/dados/Restinga_GNINA_Alvos_V4/README.md). The criterion used the highest-ranked Energy pose and fixed-frame heavy-atom RMSD. Five systems did not pass. |
| **3 · Score adjustment and sensitivity** | [Data audit](Supplementary_Data_7_Plants_GNINA_Completo/auditoria_verificada/), [archived commands](Supplementary_Data_7_Plants_GNINA_Completo/informacoes_suplementares/Supplementary_Data_S9_Archived_Commands.csv), [R2 sensitivity code and results](sensitivity-r2/) |
| **4 · ADMET** | [Execution guide](admet/ADMET_Plantas_Local_v1/README_PT.md), [CPU reference run](admet/ADMET_Plantas_Local_v1/validation_reference/), [reference protocol](admet/ADMET_Plantas_Local_v1/validation_reference/protocol.json) and [environment](admet/ADMET_Plantas_Local_v1/validation_reference/environment_freeze.txt) |
| **5 · Figures and coordinates** | [Historical figure scripts](Supplementary_Data_7_Plants_GNINA_Completo/fontes_figuras/), [R2 annotation scripts](figures-r2/), [minimum residue distances](figures-r2/residue_proximities.csv). The archive contains the PNG/SVG outputs and original PyMOL sessions. |

## Download and verify the complete archive

Use **Code → Download ZIP**, then extract the repository ZIP, or clone the repository. From the repository root:

```sh
python3 archive/restore_archive.py --output Plant_Compounds_GNINA_GitHub_Deposit_R2.zip
unzip Plant_Compounds_GNINA_GitHub_Deposit_R2.zip
(cd plant-compounds-gnina && sha256sum -c --quiet SHA256SUMS)
```

The restore script checks parts and the reconstructed ZIP against the [SHA-256 manifest](archive/manifest.json). The ZIP contains all 1,403 deposited files, including the primary SDFs, prepared receptors, poses, logs, environments, scripts and 19 PyMOL sessions. Text files mirrored above are for browsing; run reproduction commands from the extracted archive, which includes referenced molecular inputs. The GNINA binary and model weights are referenced by fingerprint, not bundled.

## Interpretation and references

The adjusted rankings, 11 reported pairs and ADMET outputs are exploratory. Native ligand redocking checks pose recovery for the deposited ligand; this study did not independently validate prospective hit selection, binding, function, selectivity or safety. The manuscript covers six species with attributable constituent data; the historical archive retains the wider seven-species project name. Figure 2 has an explicitly documented script reconstruction.

[Full bibliography](REFERENCES.md) · [Earlier reference audit](Supplementary_Data_7_Plants_GNINA_Completo/auditoria_verificada/reference_audit.csv) · [Source archive guide](Supplementary_Data_7_Plants_GNINA_Completo/README.md)
