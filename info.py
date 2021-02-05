######################
# OUTFLOW COMPARISON #
######################

# This script defines some details, e.g. which datasets are available and projection info of NGC253.


####################################################################################################

# directories
#############

subprojectdir = 'XXXXX'
rawdir  = os.path.join(subprojectdir, '01.raw_data')
datdir  = os.path.join(subprojectdir, '02.data')
momdir  = os.path.join(subprojectdir, '03.moments')
pvdir   = os.path.join(subprojectdir, '04.pVs')
sepdir  = os.path.join(subprojectdir, '05.disk_non-disk_mask')
dnddir  = os.path.join(subprojectdir, '06.disk_non-disk')
ratedir = os.path.join(subprojectdir, '07.mass_outflow_rate')
resultdir = os.path.join(subprojectdir, '08.results')
varmdir = os.path.join(subprojectdir, '09.mask_variation')
plotdir = 'XXXXX'



# available tracers
###################

datasets = [{'line': 'CO_1-0', 'cube': 'NGC253.CO_1-0', 'rms': 0.075*u.K, 'Xco_corr': 1.00, 'restfreq': 115.27120180*u.GHz, 'plot_color': 'aqua'},
            {'line': 'CO_2-1', 'cube': 'NGC253.CO_2-1', 'rms': 0.029*u.K, 'Xco_corr': 0.80, 'restfreq': 230.53800000*u.GHz, 'plot_color': 'dodgerblue'},
            {'line': 'CO_3-2', 'cube': 'NGC253.CO_3-2', 'rms': 0.36*u.K, 'Xco_corr': 0.67, 'restfreq': 345.79598990*u.GHz, 'plot_color': 'midnightblue'}
           ]
#datasets = [{'line': 'CO_3-2', 'cube': 'NGC253.CO_3-2', 'rms': 0.36*u.K, 'Xco_corr': 0.67, 'restfreq': 345.79598990*u.GHz, 'plot_color': 'midnightblue'}]


# NGC253 details
################

center       = kin_center
disk_PA      = 55.*u.degree


# slicing details
#################

slicewidth   = 5.0*u.arcsec
major_length = 50.*u.arcsec
minor_length = 50.*u.arcsec
nslices_major = 2*int(minor_length.value/slicewidth.value)-1
nslices_minor = 2*int(major_length.value/slicewidth.value)-1
SNRs = [3.0, 5.0, 10.0, 25.0, 50.0, 75.0, 100.0]


# keep debug info?
##################

debug = True


# plot setup
############

if 'aplpy_plotting' in sys.modules:
    ap.tick_label_xformat = 'hh:mm:ss.s'
    ap.tick_label_yformat = 'dd:mm:ss.s'
    ap.ticks_xspacing = Angle('0 0 1.0', unit='hourangle')
    ap.ticks_yspacing = 10.0*u.arcsec
    ap.ticks_minor_frequency = 5

if 'easy_aplpy' in sys.modules:
    easy_aplpy.settings.tick_label_xformat = 'hh:mm:ss.s'
    easy_aplpy.settings.tick_label_yformat = 'dd:mm:ss.s'
    easy_aplpy.settings.ticks_xspacing = Angle('0 0 1.0', unit='hourangle')
    easy_aplpy.settings.ticks_yspacing = 10.0*u.arcsec
    easy_aplpy.settings.ticks_minor_frequency = 5
