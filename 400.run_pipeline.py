########################
# PAPER PROJECT SCRIPT #
########################

# meta script to run all scripts of this project (project code 400)

execfile('401.prepare_images.py')
execfile('402.mask_images.py')
execfile('403.moment_maps.py')
execfile('404.regrid_masks.py')
execfile('405.separate_disk_other.py')
execfile('406.mask_moments.py')
execfile('407.mask_mass.py')
execfile('408.plot_disk_nondisk.py')
execfile('409.mass_outflow_rate.py')
execfile('410.get_statistics.py')
