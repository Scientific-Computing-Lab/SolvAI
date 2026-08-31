reinitialize
load nma_water_shell.pdb, system
select solute, resn UNL
select shellwater, byres ((not solute) within 4.5 of solute)
hide everything, all
show sticks, solute
show surface, solute
show sticks, shellwater
show spheres, shellwater and elem O
set stick_radius, 0.10, solute
set stick_radius, 0.045, shellwater
set sphere_scale, 0.14, shellwater and elem O
set transparency, 0.62, solute
set surface_color, gray80, solute
color gray35, solute and elem C
color marine, solute and elem N
color firebrick, elem O
color gray85, elem H
bg_color white
set ray_opaque_background, on
set antialias, 2
set orthoscopic, 1
set depth_cue, 0
orient solute
zoom shellwater, 1.0
ray 1400, 1050
png nma_water_shell_pymol.png, dpi=300
quit
