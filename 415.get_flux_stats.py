########################
# CASA ANALYSIS SCRIPT #
########################

# Get estimates for various quantities for disk (inside disk mask) and the rest (mainly outflows).

###################################################################################################

# import required modules
execfile('scripts/casa_imports.py')
execfile('NGC253/project_info.py')


###################################################################################################

# load sub project info
execfile(os.path.join(projectdir, 'info.py'))


###################################################################################################

# get flux/luminosity statistics
################################

fluxes       = {}
luminosities = {}
for dataset in datasets:
    fluxes[dataset['line']]       = {}
    luminosities[dataset['line']] = {}

    head = fits.getheader(os.path.join(datdir, dataset['cube']+'.fits'))
    l1pc = (distance*np.sin((u.Quantity(str(abs(head['cdelt1']))+head['cunit1']).to(u.radian)).value)).to(u.parsec)
    l2pc = (distance*np.sin((u.Quantity(str(abs(head['cdelt2']))+head['cunit2']).to(u.radian)).value)).to(u.parsec)
    Apix = l1pc*l2pc
    b1pc = (distance*np.sin(((abs(head['bmaj'])*u.degree).to(u.radian)).value)).to(u.parsec)
    b2pc = (distance*np.sin(((abs(head['bmin'])*u.degree).to(u.radian)).value)).to(u.parsec)
    Abeam = 1.13309*b1pc*b2pc
    beam_per_pix = Apix/Abeam

    for ax_type in ['major','minor']:
        fluxes[dataset['line']][ax_type]       = {}
        luminosities[dataset['line']][ax_type] = {}
        for SNR in SNRs:

            mom0 = os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.'+kin_type+'.'+str(SNR)+'s.mom0')
            fluxstat = imstat(imagename=mom0, axes=-1, algorithm='classic', clmethod='auto')
            fluxes[dataset['line']][ax_type][SNR]       = fluxstat['sum'][0]*u.K*u.km/u.s
            luminosities[dataset['line']][ax_type][SNR] = fluxstat['sum'][0]*u.K*u.km/u.s*Apix

fnpickle(fluxes, 'fluxes.pickle')
fnpickle(luminosities, 'luminosities.pickle')


###################################################################################################

# get mass statistics
#####################

masses = {}
for dataset in datasets:
    masses[dataset['line']] = {}

    for ax_type in ['major','minor']:
        masses[dataset['line']][ax_type] = {}
        for SNR in SNRs:

            mfile = os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.'+kin_type+'.'+str(SNR)+'s.mass')
            mass  = imstat(imagename=mfile, mask='', axes=-1, algorithm='classic', clmethod='auto')['sum'][0]*u.Msun

            if ( kin_type == 'non-disk'):
                try:
                    M_real    = imstat(imagename=mfile, mask='"'+os.path.join(ratedir, 'nondisk.'+dataset['line']+'.real_outflow_new.mask.fits"'), stretch=True, axes=-1, algorithm='classic', clmethod='auto')['sum'][0]*u.Msun
                    M_bubble  = imstat(imagename=mfile, mask='"'+os.path.join(ratedir, 'nondisk.'+dataset['line']+'.superbubble_new.mask.fits"'), stretch=True, axes=-1, algorithm='classic', clmethod='auto')['sum'][0]*u.Msun
                    M_codisk  = imstat(imagename=mfile, mask='"'+os.path.join(ratedir, 'nondisk.'+dataset['line']+'.cospatial_disk_new.mask.fits"'), stretch=True, axes=-1, algorithm='classic', clmethod='auto')['sum'][0]*u.Msun
                except:
                    M_real   = np.nan*u.Msun
                    M_bubble = np.nan*u.Msun
                    M_codisk = np.nan*u.Msun

            masses[dataset['line']][ax_type][SNR] = [mass, M_real, M_bubble, M_codisk]

fnpickle(masses, 'masses.pickle')


###################################################################################################
