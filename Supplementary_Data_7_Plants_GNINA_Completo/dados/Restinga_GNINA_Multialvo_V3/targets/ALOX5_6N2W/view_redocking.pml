reinitialize
load targets/ALOX5_6N2W/receptor_H.pdb, receptor
load targets/ALOX5_6N2W/native_reference.sdf, reference
load results/ALOX5_6N2W/redocking/ALOX5_6N2W_native_seed101.sdf, redocking
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
