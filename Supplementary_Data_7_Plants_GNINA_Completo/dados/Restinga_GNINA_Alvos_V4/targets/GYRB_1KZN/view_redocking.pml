# Execute from project root: @targets/GYRB_1KZN/view_redocking.pml
# Actual coordinates; no align, fit, or super command.
reinitialize
load targets/GYRB_1KZN/receptor_H.pdb, receptor
load targets/GYRB_1KZN/native_reference.sdf, crystal
load targets/GYRB_1KZN/native_top1_seed101.sdf, redocking
hide everything
remove hydro
show cartoon, receptor
show sticks, crystal or redocking
show sticks, receptor within 4 of crystal
color gray80, receptor
color marine, crystal
color tv_orange, redocking
set stick_radius, 0.16
set cartoon_transparency, 0.55
bg_color white
orient crystal
zoom crystal, 8
set ray_opaque_background, on
