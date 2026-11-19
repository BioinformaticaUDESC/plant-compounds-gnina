reinitialize
load targets/ACHE_4EY7/receptor_H.pdb, receptor
load targets/ACHE_4EY7/native_reference.sdf, reference
load results/ACHE_4EY7/redocking/ACHE_4EY7_native_seed101.sdf, redocking
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
