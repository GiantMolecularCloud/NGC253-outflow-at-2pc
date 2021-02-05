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

major_slices  = fnunpickle('pickle_jar/slices_major.pickle')
minor_slices  = fnunpickle('pickle_jar/slices_minor.pickle')

fluxes        = fnunpickle('pickle_jar/fluxes.pickle')
luminosities  = fnunpickle('pickle_jar/luminosities.pickle')
masses        = fnunpickle('pickle_jar/masses.pickle')
outflow_rates = fnunpickle('pickle_jar/outflow_rates.5.0.pickle')
energies      = fnunpickle('pickle_jar/energies.pickle')
momenta       = fnunpickle('pickle_jar/momenta.pickle')


###################################################################################################

stats = [{'line': dataset['line'], 'SNR': SNR, 'kin_type': kin_type, 'ax_type': ax_type} for dataset in datasets for SNR in SNRs for kin_type in ['disk','non-disk'] for ax_type in ['major','minor']]

def stat_idx(stats, dataset, SNR, kin_type, ax_type):
    #return stats.index({'line': dataset['line'], 'SNR': SNR, 'kin_type': kin_type, 'ax_type': ax_type})

    for idx,stat in enumerate(stats):
        if ( stat['line'] == dataset['line'] ):
            if ( stat['SNR'] == SNR ):
                if ( stat['kin_type'] == kin_type ):
                    if ( stat['ax_type'] == ax_type ):
                        return idx


###################################################################################################

# get fluxes, masses, ... in disk and other
###########################################

for dataset in datasets:
    for SNR in SNRs:
        for kin_type in ['disk','non-disk']:
            for ax_type in ['major','minor']:

                print("Getting statistics: "+dataset['line']+" "+str(SNR)+"sigma "+kin_type+" "+ax_type, end='\r')

                # some quantities are only available for non-disk. They get masked in this huge if/else block:
                if ( kin_type == 'non-disk'):
                    OR        = outflow_rates[dataset['line']][ax_type][SNR][0]*u.Msun/u.yr
                    OR_real   = outflow_rates[dataset['line']][ax_type][SNR][3]*u.Msun/u.yr
                    OR_bubble = outflow_rates[dataset['line']][ax_type][SNR][4]*u.Msun/u.yr
                    E         = energies[dataset['line']][ax_type][SNR][0]*1e51*u.erg
                    E_real    = energies[dataset['line']][ax_type][SNR][3]*1e51*u.erg
                    E_bubble  = energies[dataset['line']][ax_type][SNR][4]*1e51*u.erg
                    P         = momenta[dataset['line']][ax_type][SNR][0]*u.Msun*u.km/u.s
                    P_real    = momenta[dataset['line']][ax_type][SNR][3]*u.Msun*u.km/u.s
                    P_bubble  = momenta[dataset['line']][ax_type][SNR][4]*u.Msun*u.km/u.s
                else:
                    M_real    = np.nan*u.Msun
                    M_bubble  = np.nan*u.Msun
                    M_codisk  = np.nan*u.Msun
                    OR        = np.nan*u.Msun/u.yr
                    OR_real   = np.nan*u.Msun/u.yr
                    OR_bubble = np.nan*u.Msun/u.yr
                    E         = np.nan*u.erg
                    E_real    = np.nan*u.erg
                    E_bubble  = np.nan*u.erg
                    P         = np.nan*u.Msun*u.km/u.s
                    P_real    = np.nan*u.Msun*u.km/u.s
                    P_bubble  = np.nan*u.Msun*u.km/u.s
                idx = stat_idx(stats, dataset, SNR, kin_type, ax_type)

                # store values
                stats[idx]['flux']      = fluxes[dataset['line']][ax_type][SNR]
                stats[idx]['lum']       = luminosities[dataset['line']][ax_type][SNR]
                stats[idx]['M']         = masses[dataset['line']][ax_type][SNR][0]
                stats[idx]['M_real']    = masses[dataset['line']][ax_type][SNR][1]
                stats[idx]['M_bubble']  = masses[dataset['line']][ax_type][SNR][2]
                stats[idx]['M_codisk']  = masses[dataset['line']][ax_type][SNR][3]
                stats[idx]['OR']        = OR
                stats[idx]['OR_real']   = OR_real
                stats[idx]['OR_bubble'] = OR_bubble
                stats[idx]['E']         = E
                stats[idx]['E_real']    = E_real
                stats[idx]['E_bubble']  = E_bubble
                stats[idx]['P']         = P
                stats[idx]['P_real']    = P_real
                stats[idx]['P_bubble']  = P_bubble


###################################################################################################

# write to disk
###############

# machine readable format
fnpickle(stats, 'statistics.pickle')

# human readable format
f = open("statistics.txt","w")
f.write('FLUXES (K km/s)\n~~~~~~~~~~~~~\n')
f.write("{:<6}".format('mask')+"{:>4}".format('SNR')+"{:^26}".format('CO (1-0)')+"  "+"{:^26}".format('CO (2-1)')+"  "+"{:^26}".format('CO (3-2)')+"\n")
f.write("{:<6}".format('')+"{:<4}".format('')+("{:>10}".format('disk')+"{:>10}".format('non-disk')+"{:>6}".format('ratio')+"  ")*3+"\n")
for SNR in SNRs:
    f.write("{:<6}".format('major')+"{:4.0f}".format(SNR))
    for dataset in datasets:
        idx_disk = stat_idx(stats, dataset, SNR, 'disk', 'major')
        idx_nond = stat_idx(stats, dataset, SNR, 'non-disk', 'major')
        disk = stats[idx_disk]['flux'].value
        nond = stats[idx_nond]['flux'].value
        f.write("{:10.2e}".format(disk))
        f.write("{:10.2e}".format(nond))
        f.write("{:6.3f}".format(nond/disk))
        f.write("  ")
    f.write("\n")
for SNR in SNRs:
    f.write("{:<6}".format('minor')+"{:4.0f}".format(SNR))
    for dataset in datasets:
        idx_disk = stat_idx(stats, dataset, SNR, 'disk', 'minor')
        idx_nond = stat_idx(stats, dataset, SNR, 'non-disk', 'minor')
        disk = stats[idx_disk]['flux'].value
        nond = stats[idx_nond]['flux'].value
        f.write("{:10.2e}".format(disk))
        f.write("{:10.2e}".format(nond))
        f.write("{:6.3f}".format(nond/disk))
        f.write("  ")
    f.write("\n")
f.write("\n")

f.write('LUMINOSITIES (K km/s pc^2)\n~~~~~~~~~~~~~\n')
f.write("{:<6}".format('mask')+"{:>4}".format('SNR')+"{:^26}".format('CO (1-0)')+"  "+"{:^26}".format('CO (2-1)')+"  "+"{:^26}".format('CO (3-2)')+"\n")
f.write("{:<6}".format('')+"{:<4}".format('')+("{:>10}".format('disk')+"{:>10}".format('non-disk')+"{:>6}".format('ratio')+"  ")*3+"\n")
for SNR in SNRs:
    f.write("{:<6}".format('major')+"{:4.0f}".format(SNR))
    for dataset in datasets:
        idx_disk = stat_idx(stats, dataset, SNR, 'disk', 'major')
        idx_nond = stat_idx(stats, dataset, SNR, 'non-disk', 'major')
        disk = stats[idx_disk]['lum'].value
        nond = stats[idx_nond]['lum'].value
        f.write("{:10.2e}".format(disk))
        f.write("{:10.2e}".format(nond))
        f.write("{:6.3f}".format(nond/disk))
        f.write("  ")
    f.write("\n")
for SNR in SNRs:
    f.write("{:<6}".format('minor')+"{:4.0f}".format(SNR))
    for dataset in datasets:
        idx_disk = stat_idx(stats, dataset, SNR, 'disk', 'minor')
        idx_nond = stat_idx(stats, dataset, SNR, 'non-disk', 'minor')
        disk = stats[idx_disk]['lum'].value
        nond = stats[idx_nond]['lum'].value
        f.write("{:10.2e}".format(disk))
        f.write("{:10.2e}".format(nond))
        f.write("{:6.3f}".format(nond/disk))
        f.write("  ")
    f.write("\n")
f.write("\n")

f.write('MASSES (Msun)\n~~~~~~~~~~~~~\n')
f.write("{:<6}".format('mask')+"{:>4}".format('SNR')+"{:^50}".format('CO (1-0)')+"{:^50}".format('CO (2-1)')+"{:^50}".format('CO (3-2)')+"\n")
f.write("{:<6}".format('')+"{:<4}".format('')+"{:>10}".format('disk')+"{:>10}".format('non-disk')+"{:>10}".format('outflow')+"{:>10}".format('bubble')+"{:>10}".format('codisk')+"{:>10}".format('disk')+"{:>10}".format('non-disk')+"{:>10}".format('outflow')+"{:>10}".format('bubble')+"{:>10}".format('codisk')+"{:>10}".format('disk')+"{:>10}".format('non-disk')+"{:>10}".format('outflow')+"{:>10}".format('bubble')+"{:>10}".format('codisk')+"\n")
for SNR in SNRs:
    f.write("{:<6}".format('major')+"{:4.0f}".format(SNR))
    for dataset in datasets:
        for kin_type in ['disk','non-disk']:
            idx = stat_idx(stats, dataset, SNR, kin_type, 'major')
            f.write("{:10.2e}".format(stats[idx]['M'].value))
            if (kin_type=='non-disk'):
                f.write("{:10.2e}".format(stats[idx]['M_real'].value))
                f.write("{:10.2e}".format(stats[idx]['M_bubble'].value))
                f.write("{:10.2e}".format(stats[idx]['M_codisk'].value))
    f.write("\n")
for SNR in SNRs:
    f.write("{:<6}".format('minor')+"{:4.0f}".format(SNR))
    for dataset in datasets:
        for kin_type in ['disk','non-disk']:
            idx = stat_idx(stats, dataset, SNR, kin_type, 'minor')
            f.write("{:10.2e}".format(stats[idx]['M'].value))
            if (kin_type=='non-disk'):
                f.write("{:10.2e}".format(stats[idx]['M_real'].value))
                f.write("{:10.2e}".format(stats[idx]['M_bubble'].value))
                f.write("{:10.2e}".format(stats[idx]['M_codisk'].value))
    f.write("\n")
f.write("\n")

f.write('OUTFLOW RATE (Msun/yr)\n~~~~~~~~~~~~~\n')
f.write("{:<6}".format('mask')+"{:>4}".format('SNR')+"{:^15}".format('CO (1-0)')+"{:^15}".format('CO (2-1)')+"{:^15}".format('CO (3-2)')+"\n")
f.write("{:<6}".format('')+"{:>4}".format('')+"{:^5}".format('total')+"{:^5}".format('oflow')+"{:^5}".format('bubl')+"{:^5}".format('total')+"{:^5}".format('oflow')+"{:^5}".format('bubl')+"{:^5}".format('total')+"{:^5}".format('oflow')+"{:^5}".format('bubl')+"\n")
for SNR in SNRs:
    f.write("{:<6}".format('major')+"{:4.0f}".format(SNR))
    for dataset in datasets:
        idx = stat_idx(stats, dataset, SNR, 'non-disk', 'major')
        f.write("{:5.1f}".format(stats[idx]['OR'].value))
        f.write("{:5.1f}".format(stats[idx]['OR_real'].value))
        f.write("{:5.1f}".format(stats[idx]['OR_bubble'].value))
    f.write("\n")
for SNR in SNRs:
    f.write("{:<6}".format('minor')+"{:4.0f}".format(SNR))
    for dataset in datasets:
        idx = stat_idx(stats, dataset, SNR, 'non-disk', 'minor')
        f.write("{:5.1f}".format(stats[idx]['OR'].value))
        f.write("{:5.1f}".format(stats[idx]['OR_real'].value))
        f.write("{:5.1f}".format(stats[idx]['OR_bubble'].value))
    f.write("\n")

f.write("\n")
f.write('ENERGY (erg)\n~~~~~~~~~~~~~\n')
f.write("{:<6}".format('mask')+"{:>4}".format('SNR')+"{:^30}".format('CO (1-0)')+"{:^30}".format('CO (2-1)')+"{:^30}".format('CO (3-2)')+"\n")
f.write("{:<6}".format('')+"{:>4}".format('')+"{:^10}".format('total')+"{:^10}".format('outflow')+"{:^10}".format('sbubble')+"{:^10}".format('total')+"{:^10}".format('outflow')+"{:^10}".format('sbubble')+"{:^10}".format('total')+"{:^10}".format('outflow')+"{:^10}".format('sbubble')+"\n")
for SNR in SNRs:
    f.write("{:<6}".format('major')+"{:4.0f}".format(SNR))
    for dataset in datasets:
        idx = stat_idx(stats, dataset, SNR, 'non-disk', 'major')
        f.write("{:10.2e}".format(stats[idx]['E'].value))
        f.write("{:10.2e}".format(stats[idx]['E_real'].value))
        f.write("{:10.2e}".format(stats[idx]['E_bubble'].value))
    f.write("\n")
for SNR in SNRs:
    f.write("{:<6}".format('minor')+"{:4.0f}".format(SNR))
    for dataset in datasets:
        idx = stat_idx(stats, dataset, SNR, 'non-disk', 'minor')
        f.write("{:10.2e}".format(stats[idx]['E'].value))
        f.write("{:10.2e}".format(stats[idx]['E_real'].value))
        f.write("{:10.2e}".format(stats[idx]['E_bubble'].value))
    f.write("\n")

f.write("\n")
f.write('MOMENTUM (Msun km/s)\n~~~~~~~~~~~~~\n')
f.write("{:<6}".format('mask')+"{:>4}".format('SNR')+"{:^30}".format('CO (1-0)')+"{:^30}".format('CO (2-1)')+"{:^30}".format('CO (3-2)')+"\n")
f.write("{:<6}".format('')+"{:>4}".format('')+"{:^10}".format('total')+"{:^10}".format('outflow')+"{:^10}".format('sbubble')+"{:^10}".format('total')+"{:^10}".format('outflow')+"{:^10}".format('sbubble')+"{:^10}".format('total')+"{:^10}".format('outflow')+"{:^10}".format('sbubble')+"\n")
for SNR in SNRs:
    f.write("{:<6}".format('major')+"{:4.0f}".format(SNR))
    for dataset in datasets:
        idx = stat_idx(stats, dataset, SNR, 'non-disk', 'major')
        f.write("{:10.2e}".format(stats[idx]['P'].value))
        f.write("{:10.2e}".format(stats[idx]['P_real'].value))
        f.write("{:10.2e}".format(stats[idx]['P_bubble'].value))
    f.write("\n")
for SNR in SNRs:
    f.write("{:<6}".format('minor')+"{:4.0f}".format(SNR))
    for dataset in datasets:
        idx = stat_idx(stats, dataset, SNR, 'non-disk', 'minor')
        f.write("{:10.2e}".format(stats[idx]['P'].value))
        f.write("{:10.2e}".format(stats[idx]['P_real'].value))
        f.write("{:10.2e}".format(stats[idx]['P_bubble'].value))
    f.write("\n")

f.close()
os.system('cp statistics.txt '+os.path.join(plotdir, 'statistics.txt'))


####################################################################################################

# write table for paper
#######################

# human readable format
f = open("statistics_paper.txt","w")

f.write('LUMINOSITIES (K km/s)\n~~~~~~~~~~~~~\n')
f.write("{:<6}".format('mask')+"{:>4}".format('SNR')+"{:^26}".format('CO (1-0)')+"  "+"{:^26}".format('CO (2-1)')+"  "+"{:^26}".format('CO (3-2)')+"\n")
f.write("{:<6}".format('')+"{:<4}".format('')+("{:>10}".format('disk')+"{:>10}".format('non-disk')+"{:>6}".format('f_non')+"  ")*3+"\n")
for SNR in SNRs:
    f.write("{:<6}".format('major')+"{:4.0f}".format(SNR) +' & ')
    for dataset in datasets:
        idx_disk = stat_idx(stats, dataset, SNR, 'disk', 'major')
        idx_nond = stat_idx(stats, dataset, SNR, 'non-disk', 'major')
        disk = stats[idx_disk]['lum'].value
        nond = stats[idx_nond]['lum'].value
        f.write('$'+r' \times 10^'.join(("{:7.1e}".format(disk)).split('e+0'))+'$' +' & ')
        f.write('$'+r' \times 10^'.join(("{:7.1e}".format(nond)).split('e+0'))+'$' +' & ')
        f.write("{:4.1f}".format(nond/(nond+disk)*100) +' & ')
    f.write("\\\\ \n")
f.write("\n")

f.write('MASSES (Msun)\n~~~~~~~~~~~~~\n')
f.write("{:<6}".format('mask')+"{:>4}".format('SNR')+"{:>10}".format('disk')+"{:>10}".format('non-disk')+"{:>10}".format('outflow')+"{:>10}".format('bubble')+"{:>10}".format('codisk')+"\n")
for dataset in datasets:
    f.write(r'\sidehead{\co} '+dataset['line'] +'\n')
    for SNR in SNRs:
        f.write("{:<6}".format('major')+"{:4.0f}".format(SNR) +' & ')
        for kin_type in ['disk','non-disk']:
            idx = stat_idx(stats, dataset, SNR, kin_type, 'major')
            f.write('$'+r' \times 10^'.join(("{:7.1e}".format(stats[idx]['M'].value)).split('e+0'))+'$' +' & ')
            if (kin_type=='non-disk'):
                f.write('$'+r' \times 10^'.join(("{:7.1e}".format(stats[idx]['M_real'].value)).split('e+0'))+'$' +' & ')
                f.write('$'+r' \times 10^'.join(("{:7.1e}".format(stats[idx]['M_bubble'].value)).split('e+0'))+'$' +' & ')
                f.write('$'+r' \times 10^'.join(("{:7.1e}".format(stats[idx]['M_codisk'].value)).split('e+0'))+'$' +' & ')
        f.write("\n")
    f.write("\n")
f.write("\n")

f.write('OUTFLOW RATE (Msun/yr)\n~~~~~~~~~~~~~\n')
f.write("{:<6}".format('mask')+"{:>4}".format('SNR')+"{:^21}".format('CO (1-0)')+"{:^21}".format('CO (2-1)')+"{:^21}".format('CO (3-2)')+"\n")
f.write("{:<6}".format('')+"{:<4}".format('')+"{:>10}".format('total')+"{:>10}".format('outflow')+"{:>10}".format('bubble')+"{:>10}".format('total')+"{:>10}".format('outflow')+"{:>10}".format('bubble')+"{:>10}".format('total')+"{:>10}".format('outflow')+"{:>10}".format('bubble')+"\n")
for SNR in SNRs:
    f.write("{:<6}".format('major')+"{:4.0f}".format(SNR) +' & ')
    for dataset in datasets:
        idx_nond  = stat_idx(stats, dataset, SNR, 'non-disk', 'major')
        OR        = stats[idx_nond]['OR'].value
        OR_real   = stats[idx_nond]['OR_real'].value
        OR_bubble = stats[idx_nond]['OR_bubble'].value
        f.write("{:5.1f}".format(OR) +' & ')
        f.write("{:5.1f}".format(OR_real) +' & ')
        f.write("{:5.1f}".format(OR_bubble) +' & ')
    f.write("\\\\ \n")
f.write("\n")


f.write('ENERGY (erg)\n~~~~~~~~~~~~~\n')
f.write("{:<6}".format('mask')+"{:>4}".format('SNR')+"{:^30}".format('CO (1-0)')+"{:^30}".format('CO (2-1)')+"{:^30}".format('CO (3-2)')+"\n")
f.write("{:<6}".format('')+"{:<4}".format('')+"{:>10}".format('total')+"{:>10}".format('outflow')+"{:>10}".format('bubble')+"{:>10}".format('total')+"{:>10}".format('outflow')+"{:>10}".format('bubble')+"{:>10}".format('total')+"{:>10}".format('outflow')+"{:>10}".format('bubble')+"\n")
for SNR in SNRs:
    f.write("{:<6}".format('major')+"{:4.0f}".format(SNR) +' & ')
    for dataset in datasets:
        idx = stat_idx(stats, dataset, SNR, 'non-disk', 'major')
        f.write('$'+r' \times 10^5'.join(("{:7.1e}".format(stats[idx]['E'].value)).split('e+5'))+'$' +' & ')
        f.write('$'+r' \times 10^5'.join(("{:7.1e}".format(stats[idx]['E_real'].value)).split('e+5'))+'$' +' & ')
        f.write('$'+r' \times 10^5'.join(("{:7.1e}".format(stats[idx]['E_bubble'].value)).split('e+5'))+'$' +' & ')
    f.write("\n")
f.write("\n")

f.write('MOMENTUM (Msun km/s)\n~~~~~~~~~~~~~\n')
f.write("{:<6}".format('mask')+"{:>4}".format('SNR')+"{:^30}".format('CO (1-0)')+"{:^30}".format('CO (2-1)')+"{:^30}".format('CO (3-2)')+"\n")
f.write("{:<6}".format('')+"{:<4}".format('')+"{:>10}".format('total')+"{:>10}".format('outflow')+"{:>10}".format('bubble')+"{:>10}".format('total')+"{:>10}".format('outflow')+"{:>10}".format('bubble')+"{:>10}".format('total')+"{:>10}".format('outflow')+"{:>10}".format('bubble')+"\n")
for SNR in SNRs:
    f.write("{:<6}".format('major')+"{:4.0f}".format(SNR) +' & ')
    for dataset in datasets:
        idx = stat_idx(stats, dataset, SNR, 'non-disk', 'major')
        f.write('$'+r' \times 10^'.join(("{:7.1e}".format(stats[idx]['P'].value)).split('e+0'))+'$' +' & ')
        f.write('$'+r' \times 10^'.join(("{:7.1e}".format(stats[idx]['P_real'].value)).split('e+0'))+'$' +' & ')
        f.write('$'+r' \times 10^'.join(("{:7.1e}".format(stats[idx]['P_bubble'].value)).split('e+0'))+'$' +' & ')
    f.write("\n")
f.write("\n")

f.close()
os.system('cp statistics_paper.txt '+os.path.join(plotdir, 'statistics_paper.txt'))


####################################################################################################
