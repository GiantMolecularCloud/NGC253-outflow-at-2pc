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

# get column density
####################

def get_column_density(inps):

    dataset  = inps[0]
    SNR      = inps[1]
    kin_type = inps[2]
    ax_type  = inps[3]
    masked_mom0 = os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.'+kin_type+'.'+str(SNR)+'s.mom0')
    coldens_map = os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.'+kin_type+'.'+str(SNR)+'s.coldens')

    # bmin = imhead(imagename=masked_mom0, mode='get', hdkey='bmin')['value']
    # bmaj = imhead(imagename=masked_mom0, mode='get', hdkey='bmaj')['value']
    # pix_scale = imhead(imagename=masked_mom0, mode='get', hdkey='cdelt1')['value']/(2*np.pi)*360*60*60
    # pix_per_beam = 1.13309*bmin*bmaj/pix_scale**2

    os.system('rm -rf '+coldens_map)
    immath(imagename = masked_mom0,
        outfile  = coldens_map,
        mode     = 'evalexpr',
        expr     = '0.5e20*IM0/'+str(dataset['Xco_corr'])
#        expr     = '0.5e20*IM0/'+str(dataset['Xco_corr'])+'/'+str(pix_per_beam)
        )
    imhead(imagename = coldens_map,
        mode     = 'put',
        hdkey    = 'bunit',
        hdvalue  = 'cm^-2'
        )
    exportfits(imagename = coldens_map,
        fitsimage    = coldens_map+'.fits',
        dropstokes   = True,
        dropdeg      = True,
        overwrite    = True
        )


###################################################################################################

# get masses
############

def get_masses(inps):

    dataset  = inps[0]
    SNR      = inps[1]
    kin_type = inps[2]
    ax_type  = inps[3]
    coldens_map = os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.'+kin_type+'.'+str(SNR)+'s.coldens')
    mass_map    = os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.'+kin_type+'.'+str(SNR)+'s.mass')

    os.system('rm -rf '+mass_map)
    cdelt = imhead(imagename = coldens_map,
                   mode      = 'get',
                   hdkey     = 'cdelt1'
                  )
    pix_scale = (np.abs(cdelt['value'])*u.rad).to(u.arcsec)
    area_cm2 = ( distance.to(u.cm)*np.tan(pix_scale) )**2
    factor  = ((area_cm2 * (1.36*2.0*u.u).to(u.kg) )/(1*u.Msun)).to(1.0*u.cm**2)       # atomic weight of H2: 2    # 1.36 to account for helium

    immath(imagename = coldens_map,
           mode      = 'evalexpr',
           expr      = 'IM0*'+str(factor.value),
           outfile   = mass_map
          )
    imhead(imagename = mass_map,
           mode      ='put',
           hdkey     = 'bunit',
           hdvalue   = 'M0'
          )
    exportfits(imagename = mass_map,
        fitsimage    = mass_map+'.fits',
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

# column density
pool = Pool(30)
pool.map(get_column_density, inps)

# mass
pool = Pool(30)
pool.map(get_masses, inps)


###################################################################################################
