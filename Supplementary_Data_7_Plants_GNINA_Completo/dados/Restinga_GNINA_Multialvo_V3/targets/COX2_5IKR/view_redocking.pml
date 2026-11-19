reinitialize
load targets/COX2_5IKR/receptor_5IKR_AB_COH_H.pdb, receptor
load targets/COX2_5IKR/native_reference.sdf, reference
load results/COX2_5IKR/redocking/COX2_5IKR_native_seed101.sdf, redocking
hide everything
show cartoon, receptor
color gray80, receptor
show sticks, reference or redocking
color cyan, reference
color magenta, redocking
frame 1
select site, byres (receptor within 5 of reference)
show sticks, site
zoom reference, 8
bg_color white
# No alignment or fitting: common deposited receptor frame.
