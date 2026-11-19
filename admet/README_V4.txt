# Supplementary ADMET data, version 4

The reference inference was executed on CPU in the recorded preparation environment. No run from Juan's local workstation was received. The results provided in this V4 delivery correspond to the execution recorded in `ADMET_Plantas_Local_v1/validation_reference/protocol.json`.

The package includes the original data, 59 SDF inputs, 51 chemical identities, scripts, requirements, 2,091 primary predictions, 41 charge-sensitivity outputs, molecular descriptors, and execution records. `V4_Derived_Tables` contains manuscript and supplementary tables derived from the underlying CSV files at presentation precision; the original CSV files retain full numerical precision.

V4 uses ADMET predictions as a post hoc, complementary analysis and does not alter the eleven compound–target pairs selected by docking. Volatile and non-volatile panels remain separate. The eleven regression results outside vendor-declared limits are explicitly reported in Table S13. These scores do not provide experimental validation of exposure, efficacy, or safety.

For local reproduction, see [the ADMET execution guide](ADMET_Plantas_Local_v1/README.md). GPU-driver installation or modification is not required.
