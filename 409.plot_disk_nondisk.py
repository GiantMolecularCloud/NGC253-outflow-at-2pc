##########################
# PYTHON ANALYSIS SCRIPT #
##########################

# Plot the separated disk and other components.

###################################################################################################

# import required modules
execfile('scripts/plotting_imports.py')
execfile('NGC253/project_info.py')


###################################################################################################

# load sub project info
execfile(os.path.join(projectdir, 'info.py'))

major_slices = fnunpickle('slices_major.pickle')
minor_slices = fnunpickle('slices_minor.pickle')


###################################################################################################

# plot the separated moment maps
################################

os.system('mkdir -p '+os.path.join(plotdir, 'separated_moments'))

ap._colorbar_fontsize = 12.
ap._velo_fontsize     = 12.
props   = {'boxstyle': "round", 'facecolor': "w", 'edgecolor': "black", 'linewidth': 0.5, 'alpha': 0.8}

for dataset in datasets:
    for SNR in [3.0, 5.0, 10.0, 50.0, 75.0, 100.0]:
        for ax_type in ['major','minor']:
            infiles = [os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.'+kin_type+'.'+str(SNR)+'s.mom'+str(mom)+'.fits') for kin_type in ['disk','non-disk'] for mom in [0,1,2]]
            contours = [[os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.disk.'+str(SNR)+'s.mom0.fits'), [100,250,500,1000,2000,4000], 'black'],
                        [os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.disk.'+str(SNR)+'s.mom1.fits'), [150,200,250,300,350], 'black'],
                        [os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.disk.'+str(SNR)+'s.mom2.fits'), [20,40,60,80,100], 'black'],
                        [os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.mom0.fits'), [100,250,500,1000,2000,4000], 'black'],
                        [os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.mom1.fits'), [150,200,250,300,350], 'black'],
                        [os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.mom2.fits'), [20,40,60,80,100], 'black']]

            ap.aplpy_moments_grid(infiles,
                scaling  = [[10, 8000, r'flux [K\,km\,s$^{-1}$]',  'linear', ap.viridis_cropped],
                           [100, 400,  r'v$_{rad}$ [km\,s$^{-1}$]', 'linear', 'jet'],
                           [0,   80,   r'$\sigma$ [km\,s$^{-1}$]', 'linear', 'jet']
                          ],
                labels   = [None, 'disk', None, None, 'non-disk', None],
                label_kwargs = {'bbox': props},
                figsize  = (10.0,7.5),
                recenter = [SkyCoord('00h47m33.10s -25d17m19.68s'), 2.0*u.arcmin, 2.0*u.arcmin],
                contours = contours,
                scalebar = [29.45*u.arcsec, r'500\,pc', 'bottom'],
                beam     = 'bottom left',
                out      = plotdir+'separated_moments/'+dataset['cube']+'.ppV_mask_'+ax_type+'.'+str(SNR)+'s.moments.pdf'
                )


###################################################################################################
