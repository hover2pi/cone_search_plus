import numpy as np
from astropy.table import Table, vstack
import glob
import matplotlib
from matplotlib import pyplot as plt
import sys
import ephem
import subprocess
import requests

def equatorial_deg_to_ecliptic_deg(ra, dec):
    """Convert RA and declination as decimal degrees to
    ecliptic longitude and latitude in decimal degrees"""
    eq = ephem.Equatorial(ra * ephem.degree, dec * ephem.degree)
    ec = ephem.Ecliptic(eq)
    return ec.lon / ephem.degree, ec.lat / ephem.degree

def make_reduced_table(in_table, daily_min, base_simbad_url, match_radius):
    reduc_status = False
    twomass_IDs = []
    N_cand = len(in_table)
    rr = 0
    while rr < N_cand:
        simbadurl = basesimbadurl + "%f"%in_table['RA'][rr]  + "%20" + "%f"%in_table['Dec'][rr] + "%20radius=" + "%ds"%matchradius
        f = requests.get(simbadurl)
        queryout = f.text
        splitlines = queryout.split('\n')
        targetline = False
        for ll, line in enumerate(splitlines):
            if '2MASS J' in line and ( np.abs( float(line.split('|')[2]) - in_table['K'][rr] ) <= 0.01 or \
                                       in_table['qual'][rr][2] >= 'C' ):
                targetline = line
                break
            if ll == len(splitlines)-1:
                print("WARNING: no SIMBAD match for RA Dec = %f %f and Kmag = %.2f in %d arcsec cone"%(in_table['RA'][rr],in_table['Dec'][rr],in_table['K'][rr],matchradius))
        if targetline:
            otype = targetline.split('|')[1]
            names = targetline.split('|')[0]
            simbadKmag = float(targetline.split('|')[2])
            splitnames = names.split(',')
            N_names = len(splitnames)
            for name in splitnames:
                if '2MASS' in name: twomass_name = name
                elif 'HD' in name: hd_name = name
                elif 'HIP' in name: hip_name = name
                elif 'GSC' in name: gsc_name = name
            if '**' not in otype:
                reduc_table = Table(in_table[rr])
                twomass_IDs.append(twomass_name)
                rr = rr + 1
                break
            else:
                print 'Rejecting %s %s from all-latitude reduced list due to double star classification in SIMBAD'%(names,otype)
        rr = rr + 1
    while np.sum(reduc_table['avail'],axis=0).min() < daily_min and rr < N_cand:
        simbadurl = basesimbadurl + "%f"%in_table['RA'][rr]  + "%20" + "%f"%in_table['Dec'][rr] + "%20radius=" + "%ds"%matchradius
        f = requests.get(simbadurl)
        queryout = f.text
        splitlines = queryout.split('\n')
        targetline = False
        for ll, line in enumerate(splitlines):
            if '2MASS J' in line and ( np.abs( float(line.split('|')[2]) - in_table['K'][rr] ) <= 0.01 or \
                                       in_table['qual'][rr][2] >= 'C' ):
                targetline = line
                break
            if ll == len(splitlines)-1:
                print("WARNING: no SIMBAD match for RA Dec = %f %f and Kmag = %.2f in %d arcsec cone"%(in_table['RA'][rr],in_table['Dec'][rr],in_table['K'][rr],matchradius))
        if targetline:
            otype = targetline.split('|')[1]
            names = targetline.split('|')[0]
            simbadKmag = float(targetline.split('|')[2])
            splitnames = names.split(',')
            N_names = len(splitnames)
            for name in splitnames:
                if '2MASS' in name: twomass_name = name
                elif 'HD' in name: hd_name = name
                elif 'HIP' in name: hip_name = name
                elif 'GSC' in name: gsc_name = name
            if '**' not in otype:
                reduc_table.add_row(in_table[rr])
                twomass_IDs.append(twomass_name)
            else:
                print 'Rejecting %s %s from reduced list due to double star classification in SIMBAD'%(names,otype)
        rr = rr + 1
    if rr == N_cand:
        print('WARNING: Could not reduce the target list to satisfy %d min targets per day. Will not store the attempt.' % daily_min) 
    else:
        twomass_col = Table.Column(twomass_IDs, '2MASS')
        reduc_table.add_column(twomass_col)
        reduc_table.sort(keys='eb')
        reduc_table.reverse()
        reduc_status = True
    return reduc_table, reduc_status

targets_fname = sys.argv[-2]
avail_fname = sys.argv[-1]

min_targets_per_day = 2
min_targets_per_day_hemi = 1

avail = np.load(avail_fname)
avail_1yr = avail[:,:365]

targets = Table.read(targets_fname, format='ascii', names=['RA', 'Dec', 'J', 'H', 'K', 'qual', 'idx'])
N_targets_full = len(targets)

assert N_targets_full == avail.shape[0]

el = []
eb = []
for row in targets:
    l, b = equatorial_deg_to_ecliptic_deg(row['RA'], row['Dec'])
    el.append(l)
    eb.append(b)
el_col = Table.Column(el, 'el')
eb_col = Table.Column(eb, 'eb')
avail_col = Table.Column(avail_1yr, 'avail')
N_days_per_target = np.sum(avail_1yr, axis=1)
N_days_col = Table.Column(N_days_per_target, 'N_days')
targets.add_column(el_col, index=2)
targets.add_column(eb_col, index=3)
targets.add_column(avail_col) 
targets.add_column(N_days_col)

targets_eN = targets[targets['eb'] >= 0]
targets_eS = targets[targets['eb'] < 0]
N_targets_eN = len(targets_eN)
N_targets_eS = len(targets_eS)

targets.sort(keys='N_days')
targets.reverse()
targets_eN.sort(keys='N_days')
targets_eN.reverse()
targets_eS.sort(keys='N_days')
targets_eS.reverse()

total_avail_all = np.sum(avail_1yr, axis=0)
min_avail_all = total_avail_all.min()
min_avail_eS = np.sum(targets_eS['avail'],axis=0).min()
min_avail_eN = np.sum(targets_eN['avail'],axis=0).min()

print("Starting from a list of %d target stars." % N_targets_full)
print("On any day of the year, at least %d targets are available at all ecliptic latitudes." % min_avail_all)
print("On any day of the year, at least %d targets are available in the northern ecliptic hemisphere." % min_avail_eN)
print("On any day of the year, at least %d targets are available in the southern ecliptic hemisphere." % min_avail_eS)

#basesimbadurl = "http://simbad.u-strasbg.fr/simbad/sim-script?script="
basesimbadurl = "http://simbad.cfa.harvard.edu/simbad/sim-script?script="
basesimbadurl = basesimbadurl + "format%20object%20%22%25IDLIST%28A;1,NAME,2MASS,HD,HIP,GSC%29%20|%25OTYPELIST%20|%20%25FLUXLIST%28K;F%29%22"
basesimbadurl = basesimbadurl + "%0Aoutput%20console=off%20script=off%0Aquery%20coo%20"
matchradius = 2

print('Reducing all-latitude list...')
reduc_targets, reduc_status = make_reduced_table(targets, min_targets_per_day, basesimbadurl, matchradius)
print('Reducing northern latitude list...')
reduc_targets_eN, reduc_status_eN = make_reduced_table(targets_eN, min_targets_per_day_hemi, basesimbadurl, matchradius)
print('Reducing southern latitude list...')
reduc_targets_eS, reduc_status_eS = make_reduced_table(targets_eS, min_targets_per_day_hemi, basesimbadurl, matchradius)

if reduc_status:
    print('\nReduced list, all latitudes (%d stars for min daily avail. %d stars)'%(len(reduc_targets),min_targets_per_day))
    print reduc_targets
    reduc_targets_full = reduc_targets['RA', 'Dec', 'J', 'H', 'K', 'qual', 'idx']
    reduc_targets_RADec = reduc_targets['RA', 'Dec']
    reduc_targets_twomass = reduc_targets['2MASS']
    reduc_targets_RADec.meta = {} 
    reduc_targets_twomass.meta = {} 
    reduc_targets_fname = targets_fname + '_reduc'
    reduc_targets_RADec_fname = reduc_targets_fname + '_RADec'
    reduc_targets_twomass_fname = reduc_targets_fname + '_2MASS'
    reduc_targets_full.write(reduc_targets_fname, format='ascii.no_header')
    reduc_targets_RADec.write(reduc_targets_RADec_fname, format='ascii.no_header')
    reduc_targets_twomass.tofile(reduc_targets_twomass_fname,sep="\n")
    print('Wrote all-latitude reduced target lists to\n%s,\n%s,\n%s'%(reduc_targets_fname, reduc_targets_RADec_fname, reduc_targets_twomass_fname))
if reduc_status_eN:
    print('\nReduced list, northern latitudes (%d stars for min daily avail. %d stars per hem.)'%(len(reduc_targets_eN),min_targets_per_day_hemi))
    print reduc_targets_eN
    reduc_targets_eN_full = reduc_targets_eN['RA', 'Dec', 'J', 'H', 'K', 'qual', 'idx']
    reduc_targets_eN_RADec = reduc_targets_eN['RA', 'Dec']
    reduc_targets_eN_twomass = reduc_targets_eN['2MASS']
    reduc_targets_eN_RADec.meta = {}
    reduc_targets_eN_twomass.meta = {}
    reduc_targets_eN_fname = targets_fname + '_north_reduc'
    reduc_targets_eN_RADec_fname = reduc_targets_eN_fname + '_RADec'
    reduc_targets_eN_twomass_fname = reduc_targets_eN_fname + '_2MASS'
    reduc_targets_eN_full.write(reduc_targets_eN_fname, format='ascii.no_header')
    reduc_targets_eN_RADec.write(reduc_targets_eN_RADec_fname, format='ascii.no_header')
    reduc_targets_eN_twomass.tofile(reduc_targets_eN_twomass_fname,sep="\n")
    print('Wrote northern latitude reduced target lists to\n%s,\n%s,\n%s'%(reduc_targets_eN_fname, reduc_targets_eN_RADec_fname, reduc_targets_eN_twomass_fname))
if reduc_status_eS:
    print ('\nReduced list, southern latitudes (%d stars for min daily avail. %d stars per hem.)'%(len(reduc_targets_eS),min_targets_per_day_hemi))
    print reduc_targets_eS
    reduc_targets_eS_full = reduc_targets_eS['RA', 'Dec', 'J', 'H', 'K', 'qual', 'idx']
    reduc_targets_eS_RADec = reduc_targets_eS['RA', 'Dec']
    reduc_targets_eS_twomass = reduc_targets_eS['2MASS']
    reduc_targets_eS_RADec.meta = {} 
    reduc_targets_eS_twomass.meta = {} 
    reduc_targets_eS_fname = targets_fname + '_south_reduc'
    reduc_targets_eS_RADec_fname = reduc_targets_eS_fname + '_RADec'
    reduc_targets_eS_twomass_fname = reduc_targets_eS_fname + '_2MASS'
    reduc_targets_eS_full.write(reduc_targets_eS_fname, format='ascii.no_header')
    reduc_targets_eS_RADec.write(reduc_targets_eS_RADec_fname, format='ascii.no_header')
    reduc_targets_eS_twomass.tofile(reduc_targets_eS_twomass_fname,sep="\n")
    print('Wrote southern latitude reduced target lists to\n%s,\n%s,\n%s'%(reduc_targets_eS_fname, reduc_targets_eS_RADec_fname, reduc_targets_eS_twomass_fname))
if reduc_status_eN and reduc_status_eS:
    reduc_targets_eNS = vstack([reduc_targets_eN, reduc_targets_eS])
    print ('\nReduced list, combined northern+southern (%d stars for min daily avail. %d stars per hem.)'%(len(reduc_targets_eNS),min_targets_per_day_hemi))
    reduc_targets_eNS_full = reduc_targets_eNS['RA', 'Dec', 'J', 'H', 'K', 'qual', 'idx']
    reduc_targets_eNS_RADec = reduc_targets_eNS['RA', 'Dec']
    reduc_targets_eNS_twomass = reduc_targets_eNS['2MASS']
    reduc_targets_eNS_RADec.meta = {}
    reduc_targets_eNS_twomass.meta = {} 
    reduc_targets_eNS_fname = targets_fname + '_NS_reduc'
    reduc_targets_eNS_RADec_fname = reduc_targets_eNS_fname + '_RADec'
    reduc_targets_eNS_twomass_fname = reduc_targets_eNS_fname + '_2MASS'
    reduc_targets_eNS_full.write(reduc_targets_eNS_fname, format='ascii.no_header')
    reduc_targets_eNS_RADec.write(reduc_targets_eNS_RADec_fname, format='ascii.no_header')
    reduc_targets_eNS_twomass.tofile(reduc_targets_eNS_twomass_fname,sep="\n")
    print('Wrote combined northern+southern reduced target lists to\n%s,\n%s,\n%s'%(reduc_targets_eNS_fname, reduc_targets_eNS_RADec_fname, reduc_targets_eNS_twomass_fname))

#plt.figure(figsize=(8,8))
#plt.imshow(reduc_targets['avail'], aspect='auto', interpolation='nearest')
##plt.imshow(targets_eS['avail'], aspect='auto', interpolation='nearest')
#plt.show()
