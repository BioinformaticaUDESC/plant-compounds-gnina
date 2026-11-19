# GNINA multi-target V3: technical results

The library comprises 40 literature-derived compounds, represented by 51 declared structures and 59 verified conformational inputs. Four targets were assessed by native-ligand redocking; three satisfied the prespecified criterion for candidate screening. Results are exploratory for isolated molecules. No project-specific chemical analysis confirmed the occurrence of any compound in study samples.

## Native-ligand redocking

| Target | Top-ranked Energy pose with RMSD ≤2 Å | Median RMSD (Å) | Status |
| --- | ---: | ---: | --- |
| ACHE_4EY7 | 5/5 | 1.070 | Pass |
| ALOX5_6N2W | 0/5 | 3.549 | Fail |
| CB2_5ZTY | 5/5 | 0.537 | Pass |
| COX2_5IKR | 5/5 | 0.481 | Pass |

The primary analysis contains 10,620 candidate-docking poses and 400 native-control poses. A total of 10,980 candidate-docking poses were exported: 180 were quarantined because of a stereochemical identity discrepancy and 180 belonged to a separately retained methodological correction. Candidate screening for accepted targets used seeds 42, 142, and 242, with up to 20 poses per input and seed.

One inherited caryophyllene-oxide input had a three-dimensional stereochemical identity that differed from its declared structure. It was excluded from the primary analysis. The verified replacement was retained as a separate diagnostic because it was not sufficiently conformationally independent from two already retained correct conformers. The primary library therefore contains 59 verified inputs.

For 5-LOX, poses within 2 Å of the crystal structure were sampled, but the first Energy-ranked pose failed in each of the five seeds. The target was excluded from candidate screening under the recorded protocol, without retrospectively selecting a pose to pass the control. This outcome does not demonstrate absence of biological activity.

## Botanical coverage

| Species | Selected compounds | Archive note |
| --- | ---: | --- |
| *Varronia curassavica* | 2 | Literature selection not exhaustive |
| *Campomanesia littoralis* | 0 | No sufficiently supported nominal chemical assignment; no generic *C. xanthocarpa* substitution |
| *Myrcia palustris* | 3 | Literature selection not exhaustive |
| *Eugenia catharinae* | 8 | Extract comparators only; volatile composition pending |
| *Schinus terebinthifolia* | 8 | Literature selection not exhaustive |
| *Psidium guajava* | 11 | Literature selection not exhaustive |
| *Lippia alba* | 19 | Literature selection not exhaustive |

No project-specific GC–MS data were available. *E. catharinae* is represented by extract comparators; indexed fragments from its essential-oil table were retained as literature leads rather than as a confirmed complete composition. *C. littoralis* remains without a chemical assignment.

## Methodological record

GNINA 1.3.3 was executed on CPU with CNN fast rescoring, Energy as the prespecified primary ranking, exhaustiveness 16, and up to 20 poses. The same target-specific scoring protocol was used for redocking and candidate screening. CNN fast was selected to make the expanded screen computationally tractable; its values should not be compared numerically with standard-ensemble values from earlier deliveries.

Native-control inputs were prepared independently of crystal-ligand coordinates. The prespecified criterion was a symmetry-aware heavy-atom RMSD ≤2 Å for the first Energy-ranked pose in at least four of five seeds, evaluated in the fixed protein reference frame without ligand alignment. Chemical identities of controls were checked against RCSB CCD records, including the documented protonation state of donepezil.

Receptors were rigid; deposited heavy-atom coordinates were retained and hydrogens were modelled at pH 7.4. Waters and glycans were excluded. The CB2 construct retains the mutant T4-lysozyme fusion without membrane or lipids and was deposited in an antagonist-bound state. The 5-LOX model retains iron under an empirical treatment and the deposited disorder of the site-covering region.

The audit evaluates chemical identity, finite coordinates, non-metal contacts below 2 Å, distances below 70% of the summed van der Waals radii, positive empirical energy, box escape >0.1 Å, bond-length change >0.15 Å, and serious intramolecular overlap between atoms separated by more than three bonds. Metal contacts are recorded separately. Absence of an audit flag is a limited quality check and does not establish binding. Residue proximities ≤4 Å are geometric observations and were not automatically labelled as hydrogen bonds or mechanisms.

## Interpretation

The expansion provides traceable structural hypotheses and enables assessment of sensitivity to random seed, conformer, and pose ranking. A failing native control was retained and prevented screening for that target. The data do not provide a potency ranking across plants or proteins and do not estimate mixture activity.

Docking scores do not establish anti-inflammatory, antimicrobial, anticancer, antioxidant, or antidepressant activity. Predicted CB2 binding does not distinguish agonism from antagonism. Calculated physicochemical descriptors are not experimental assays or validated predictions of toxicity, bioavailability, or safety.

## Archive contents

The PDF contains the visual and methodological summary. Complete tables are under `reports/` and `results/<target>/`. Chemical identities and literature sources are recorded in `structures.csv` and `literature_species_links.csv`. `target_rationale.json` documents target-selection reasoning and structural caveats. Target-specific PyMOL scripts use the same unaligned reference frame.
