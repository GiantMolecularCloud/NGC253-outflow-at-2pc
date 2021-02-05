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

# get moments
#############

def masked_moments(inps):

    dataset  = inps[0]
    SNR      = inps[1]      # 3.0 / 5.0 / 10.0 / 25.0 / 50.0 / 75.0 / 100.0
    kin_type = inps[2]      # disk / non-disk
    ax_type  = inps[3]      # major / minor

    masked_cube = os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.'+kin_type+'.'+str(SNR)+'s.image')

    for mom in [0,1,2,4]:
        mom_map = os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.'+kin_type+'.'+str(SNR)+'s.mom'+str(mom))
        os.system('rm -rf '+mom_map)
        immoments(imagename = masked_cube,
            outfile    = mom_map,
            moments    = [mom],
            axis       = 'spectral',
            region     = '',
            chans      = '',
            includepix = [0,1000]       # is needed for moment 4 despite the map being positive already
            )
        exportfits(imagename = mom_map,
            fitsimage    = mom_map+'.fits',
            dropstokes   = True,
            dropdeg      = True,
            overwrite    = True
            )


###################################################################################################

# execute in parallel
#####################

inps = []
for dataset in datasets:
    for SNR in SNRs:
        for kin_type in ['disk','non-disk']:
            for ax_type in ['major','minor']:
                inps.append([dataset, SNR, kin_type, ax_type])

pool = Pool(30)
pool.map(masked_moments, inps)


###################################################################################################

# mask 0.0 values
#################

mom_files = [os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.'+kin_type+'.'+str(SNR)+'s.mom0.fits') for ax_type in ['major','minor'] for kin_type in ['disk','non-disk'] for dataset in datasets for SNR in SNRs]

for mom_file in tqdm(mom_files):
    mom = fits.open(mom_file)[0]
    mom.data[mom.data==0.0] = np.nan
    fits.writeto(mom_file, data=mom.data, header=mom.header, overwrite=True)


####################################################################################################

# get log10
###########

mom_files = [os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.'+kin_type+'.'+str(SNR)+'s.mom0.fits') for ax_type in ['major','minor'] for kin_type in ['disk','non-disk'] for dataset in datasets for SNR in SNRs]
for mom_file in tqdm(mom_files):
    mom = fits.open(mom_file)[0]
    mom.data = np.log10(mom.data)
    fits.writeto(mom_file.replace('.fits','.log.fits'), data=mom.data, header=mom.header, overwrite=True)

# also get log maps of non-separated moments
mom_files = [os.path.join(momdir, dataset['cube']+'.mask_'+str(SNR)+'s.mom0.fits') for dataset in datasets for SNR in SNRs]
for mom_file in tqdm(mom_files):
    mom = fits.open(mom_file)[0]
    mom.data = np.log10(mom.data)
    fits.writeto(mom_file.replace('.fits','.log.fits'), data=mom.data, header=mom.header, overwrite=True)


###################################################################################################
