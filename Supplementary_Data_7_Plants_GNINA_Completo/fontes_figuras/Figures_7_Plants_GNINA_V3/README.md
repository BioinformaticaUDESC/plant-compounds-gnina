# PyMOL and chemical structure figures for the seven-plant GNINA study

This package visualizes completed V3/V4 docking results. No docking search, score calculation, molecular dynamics or biological experiment was performed for this figure update.

## Contents

- figures: twelve publication plates for manuscript Figures 3–5 and supplementary Figures S1–S7. Figures 3 and 4 each have two parts.
- panels: twelve native-control views and seven candidate views, ray traced in PyMOL at 1200 × 900 pixels, with 300 dpi metadata.
- sessions: nineteen editable PyMOL PSE sessions.
- structures: the twelve prepared receptors, native-reference ligands, seed-101 top1 redocking poses, and seven selected seed-42 candidate poses.
- molecules: individual PNG and vector SVG depictions for all 51 declared stereochemical structures, representing 40 named compounds.
- figure_manifest.json: selection rules, exact input identifiers, individual-pose scores, native RMSDs, source fingerprints and software versions.
- molecule_depictions.json: complete correspondence between chemical figure IDs C01–C51, declared structures, InChIKeys, and original isomeric SMILES.
- session_coordinate_audit.json: all nineteen sessions checked against exported original heavy-atom coordinates.

## Presentation rules

Native overlays consistently show Energy top1 at seed 101, including the failed TyrS seed and all five system-level failures. The complete acceptance rule remains RMSD <= 2 Å in at least four of five seeds. A single displayed pose cannot change a system-level outcome.

Candidate views show the compound with the lowest unrounded median Energy in each volatile panel, using its lexicographically first primary input and seed 42. No pose was chosen by visual appearance or lowest RMSD. Candidate and native ligands have different chemistry; the native serves only as a site reference. Displayed individual-pose energies can differ from the pooled medians.

Teal denotes experimental-reference ligand carbons; magenta denotes redocked carbons; amber denotes candidate carbons. Oxygen and nitrogen use red and blue element colors. Protein context and surface are gray. No electrostatic-potential calculation or automatic hydrogen-bond interpretation is implied.

Molecular drawings preserve the declared chemical graph and stereochemistry. RDKit coordinates were generated in 2D for display, with CoordGen used for the eucalyptol cage to avoid overlapping atom labels. Every depiction passed its declared InChIKey check. Drawings do not establish chemical occurrence or stereoisomer composition in plant samples.

## Reopening and reproducing

Extract in a folder under HOME, for example ~/Figures_7_Plants_GNINA_V3. Open any sessions/*.pse file directly in PyMOL. The PSE contains the displayed objects and camera view.

To reproduce images with Python, use the versions in requirements-visualization.txt and run render_pymol.py, build_molecule_figures.py and compose_pymol_figures.py from this folder. The rendering scripts use local relative paths. Source input preparation and the full GNINA results remain in the separately supplied V4 study archive.

The rendered camera is a display transformation. Original molecular coordinates are retained; there is no independent ligand fitting. The session audit maximum coordinate difference is below 0.00002 Å and reflects single-precision storage.

## Independent authenticity verification on 27 October 2026

The verification_2026-10-27 folder records the live PubChem comparison for all 40 parent CIDs, the 51 declared chemical structures, all 59 primary input SDF files, original coordinate checks for all 19 PyMOL scenes, and pixel-identical re-renderings from all 19 saved PSE sessions. All twelve native reference ligands also match deposited PDB heavy-atom positions exactly. Original panel images and sessions are preserved. The accompanying manuscript revision changes only Table 2 to two horizontal rules without shading.
