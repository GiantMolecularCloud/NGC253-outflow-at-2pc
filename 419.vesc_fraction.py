########################
# CASA ANALYSIS SCRIPT #
########################

# Get fraction of gas that is above/below escape velocity.

###################################################################################################

# import required modules
execfile('NGC253/project_info.py')


###################################################################################################

# load sub project info
execfile(os.path.join(projectdir, 'info.py'))
ratedir = ratedir+'.improved3'

vesc = 500*u.km/u.s             # according to Walter +17


###################################################################################################
# get escape fraction
###################################################################################################

for dataset in datasets:
    mass = fits.getdata(os.path.join(ratedir, dataset['cube']+'.ppV_mask_major.non-disk.5.0s.mass_cube.fits'))
    velo = fits.getdata(os.path.join(ratedir, dataset['cube']+'.velo.fits'))

    m_esc = np.nansum(mass[velo>vesc.value])
    m_non = np.nansum(mass[velo<vesc.value])
    m_tot = np.nansum(mass)

    f_esc = m_esc/m_tot
    f_non = m_non/m_tot

    print(dataset['line']+"   escape fraction by mass: "+str(f_esc))
    print(dataset['line']+"   not escaping fraction:   "+str(f_non))


###################################################################################################
#
###################################################################################################
