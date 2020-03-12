# Forces-on-carangiform-swimmers
Calculate surface pressures and swimming forces on the bodies of carangiform swimmers (bluegill and trout). Used alongside John Dabiri's queen2 (http://dabirilab.com/software/), these prepare masks for the queen2 and process queen2's pressure field output to calculate forces. Associated with Lucas et al. (Submitted) & related data. For details on methods and validation, see Lucas et al. (2017) PLoS ONE (https://doi.org/10.1371/journal.pone.0189225).

Required Matlab scripts:

https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections

https://www.mathworks.com/matlabcentral/fileexchange/34874-interparc

https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab

Required Python script (my translation of the Matlab script; see my InterX repository)

https://github.com/kelseynlucas/InterX


Descriptions of scripts in this repository:

Matlab scripts

	fish_masks m files - used to make masks of the fish body for queen2 from manually-digitized outlines of fishes

	fish_midline_fctn.m - called by other scripts; automatically finds the midline of the fish from the fish's body outline

	height_dict_body_only.m - called by other scripts; measures body depth over body length from a lateral-view outline of the fish

	pressureintegrator_for_fish_force_calc m files - calculates the forces acting on the fish's body from queen2 pressure fields

Python scripts

	get_fish_kinematics_v5.py - reads midlines extracted by fish_midline_fctn, smooths, and calculates kinematics parameters

	assemble_fish_data_for_stats_v3.py - reads in data produced by other codes and compiles it for later statistical analyses

	assemble_fish_data_for_stats_fctns_v2.py - custom data aggregation functions called by the above


