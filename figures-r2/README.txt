IJMS 4725052 — Figures 3 and 7 (R2)

figures/ contains the four submission plates (2472 x 2333 pixels, 300 dpi).
File suffix A contains panels A–D; suffix B contains panels E–G.
panels/ contains the 14 individual annotated panels (1200 x 900 pixels).
residue_proximities.csv contains the 187 panel/residue observations in
Table S17, with distances in Å. "labelled" = 1 for the three nearest
residues shown in each panel. The distances describe geometry only.
figure_label_audit.json records source hashes, pixel anchors and the
coordinate comparison for each panel.

To reproduce the annotations, set SOURCE from this directory:
  SOURCE=../Supplementary_Data_7_Plants_GNINA_Completo/fontes_figuras/Figures_7_Plants_GNINA_V3
Run from the directory containing this readme, with Python and NumPy,
SciPy, Pillow, Matplotlib, DejaVu Sans and the PyMOL Python module
(the archived sessions were checked using PyMOL 3.2.0a0.post5859):

  python annotate_panels.py "$SOURCE" annotated
  cp compose_plates.py annotated/
  python annotated/compose_plates.py

The original PyMOL panel pixels are preserved except for the overlaid
labels and leader lines. The original PyMOL session camera is used to
project the label anchors. The manuscript's Figure 7 corresponds to
the candidate panels archived as Figure 4. Neither the structures nor
the displayed poses or seeds were changed. No new docking runs were
performed for this revision.
