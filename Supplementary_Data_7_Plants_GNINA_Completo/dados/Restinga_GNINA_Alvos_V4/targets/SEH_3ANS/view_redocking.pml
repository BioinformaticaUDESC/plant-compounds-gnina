# Execute from project root: @targets/SEH_3ANS/view_redocking.pml
# Actual coordinates; no align, fit, or super command.
reinitialize
load targets/SEH_3ANS/receptor_H.pdb, receptor
load targets/SEH_3ANS/native_reference.sdf, crystal
load targets/SEH_3ANS/native_top1_seed101.sdf, redocking
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
