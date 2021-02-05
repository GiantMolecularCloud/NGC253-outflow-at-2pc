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
# define overlap area
###################################################################################################

# regrid images to common pixel scale
imregrid(imagename = os.path.join(datdir, 'NGC253.CO_1-0'),
         template  = os.path.join(datdir, 'NGC253.CO_3-2'),
         output    = os.path.join(datdir, 'NGC253.CO_1-0.regrid_CO_3-2'),
         asvelocity = True
        )
imregrid(imagename = os.path.join(datdir, 'NGC253.CO_2-1'),
         template  = os.path.join(datdir, 'NGC253.CO_3-2'),
         output    = os.path.join(datdir, 'NGC253.CO_2-1.regrid_CO_3-2'),
         asvelocity = True
        )
exportfits(imagename = os.path.join(datdir, 'NGC253.CO_1-0.regrid_CO_3-2'),
           fitsimage = os.path.join(datdir, 'NGC253.CO_1-0.regrid_CO_3-2.fits'),
           velocity  = True,
           optical   = True,
           dropdeg   = True,
           history   = False
          )
exportfits(imagename = os.path.join(datdir, 'NGC253.CO_2-1.regrid_CO_3-2'),
           fitsimage = os.path.join(datdir, 'NGC253.CO_2-1.regrid_CO_3-2.fits'),
           velocity  = True,
           optical   = True,
           dropdeg   = True,
           history   = False
          )

# load images
co10 = fits.open(os.path.join(datdir, 'NGC253.CO_1-0.regrid_CO_3-2.fits'))[0]
co21 = fits.open(os.path.join(datdir, 'NGC253.CO_2-1.regrid_CO_3-2.fits'))[0]
co32 = fits.open(os.path.join(datdir, 'NGC253.CO_3-2.fits'))[0]

# get where image is defined
co10_def = ~np.isnan(co10.data[0])
co21_def = ~np.isnan(co21.data[0])
co32_def = ~np.isnan(co32.data[0])
fits.writeto(datdir+'/test4.fits', data=co32_def.astype(float), overwrite=True)

# overlap region
overlap = np.logical_and(co10_def, np.logical_and(co21_def, co32_def)).astype(float)
fits.writeto(os.path.join(datdir, 'NGC253.CO_overlap.fits'), data=overlap, header=co32.header, overwrite=True)

# regrid to original datasets
importfits(fitsimage = os.path.join(datdir, 'NGC253.CO_overlap.fits'),
           imagename = os.path.join(datdir, 'NGC253.CO_overlap.image')
          )
imregrid(imagename = os.path.join(datdir, 'NGC253.CO_overlap.image'),
         template  = os.path.join(momdir, 'NGC253.CO_1-0.mask_5.0s.mom0'),
         output    = os.path.join(datdir, 'NGC253.CO_overlap.CO_1-0.image'),
         asvelocity = True
        )
imregrid(imagename = os.path.join(datdir, 'NGC253.CO_overlap.image'),
         template  = os.path.join(momdir, 'NGC253.CO_2-1.mask_5.0s.mom0'),
         output    = os.path.join(datdir, 'NGC253.CO_overlap.CO_2-1.image'),
         asvelocity = True
        )
exportfits(imagename = os.path.join(datdir, 'NGC253.CO_overlap.CO_1-0.image'),
           fitsimage = os.path.join(datdir, 'NGC253.CO_overlap.CO_1-0.fits'),
           velocity  = True,
           optical   = True,
           dropdeg   = True,
           history   = False
          )
exportfits(imagename = os.path.join(datdir, 'NGC253.CO_overlap.CO_2-1.image'),
           fitsimage = os.path.join(datdir, 'NGC253.CO_overlap.CO_2-1.fits'),
           velocity  = True,
           optical   = True,
           dropdeg   = True,
           history   = False
          )
os.system('cp -r '+os.path.join(datdir, 'NGC253.CO_overlap.fits')+' '+os.path.join(datdir, 'NGC253.CO_overlap.CO_3-2.fits'))
importfits(fitsimage = os.path.join(datdir, 'NGC253.CO_overlap.CO_3-2.fits'),
           imagename = os.path.join(datdir, 'NGC253.CO_overlap.CO_3-2.image')
          )


###################################################################################################
# get flux/luminosity statistics
###################################################################################################

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

    for kin_type in ['disk','non-disk']:
        fluxes[dataset['line']][kin_type]       = {}
        luminosities[dataset['line']][kin_type] = {}
        for ax_type in ['major','minor']:
            fluxes[dataset['line']][kin_type][ax_type]       = {}
            luminosities[dataset['line']][kin_type][ax_type] = {}
            for SNR in SNRs:

                mom0 = os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.'+kin_type+'.'+str(SNR)+'s.mom0')
                mask = os.path.join(datdir, 'NGC253.CO_overlap.'+dataset['line']+'.image')
                fluxstat = imstat(imagename=mom0, mask='"'+mask+'"', stretch=True, axes=-1, algorithm='classic', clmethod='auto')
                fluxes[dataset['line']][kin_type][ax_type][SNR]       = fluxstat['sum'][0]*u.K*u.km/u.s
                luminosities[dataset['line']][kin_type][ax_type][SNR] = fluxstat['sum'][0]*u.K*u.km/u.s*Apix

fnpickle(fluxes, 'fluxes_equal_area.pickle')
fnpickle(luminosities, 'luminosities_equal_area.pickle')


###################################################################################################
# print stats
###################################################################################################

for dataset in datasets:
    print(dataset['line'])
    print('{:.2e}'.format(luminosities[dataset['line']]['disk']['major'][5.0]))
    print('{:.2e}'.format(luminosities[dataset['line']]['non-disk']['major'][5.0]))


###################################################################################################
# get area
###################################################################################################

execfile('scripts/various_helpers.py')

olim     = os.path.join(datdir, 'NGC253.CO_overlap.CO_3-2.fits')
overlap  = fits.open(olim)[0]
pix_area = (np.abs(get_axis_info(olim,1)[2]) * np.abs(get_axis_info(olim,2)[2])).to(u.arcsec**2)
npix     = np.nansum(overlap.data)
overlap_area = pix_area * npix
print(overlap_area)
print(angle_to_parsec(np.sqrt(overlap_area), 'NGC253')**2)


###################################################################################################
#
###################################################################################################
