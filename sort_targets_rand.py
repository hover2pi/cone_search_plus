'''
Neil Zimmerman
February 2016

SUMMARY 

Taking a list of targets generated by list_builder.py as input,
along with separately computed availability table (produced with availability_checker.py),
form a reduced subset of favored commissioning target candidates.

DESCRIPTION

The algorithm to form this reduced subset goes as follows:

First, the stars are ranked by duration of observing window,
over all ecliptic latitudes, and then separately for northern ecliptic latitudes 
and souther ecliptic latitudes, resulting in three tables.
Beginning with the star having the longest observing window, copy stars into the
reduced subset until the list contains sufficient year-round coverage that some 
minimum number of targets are available on any given calendar date.
This daily minimum number is hard-coded below in the main body, specified separately for
all ecliptic latitudes, and for the ecliptic hemispheres:

min_targets_per_day
and
min_targets_per_day_hemi

The subset meeting the all-latitude daily minimum availability is stored in several
table versions with different organizations and data fields:
{targets_fname}_reduc  --  in the original query_2mass output format: RA, Dec, J, H, K, qual, idx
{targets_fname}_reduc_RADec -- RA and Dec only
{targets_fname}_reduc_2MASS -- 2MASS IDs only, based on SIMBAD coordinate query matching

The table file name is augmented to reflect the specific additional criteria based on
ecliptic latitude considerations:
{tarets_fname}_north_reduc -- Northern ecliptic latitude subset meeting min_targets_per_day_hemi
{tarets_fname}_south_reduc -- Southern ecliptic latitude subset meeting min_targets_per_day_hemi
{tarets_fname}_NS_reduc -- Combined Southern and Northern subsets that simultaneously
                           satisfy min_targets_per_day_hemi in each ecliptic hemisphere.

The last variant of these tables is the most valuable one from the JWST commissioning perspective.
Note that there is no guarantee of success, since there may not be a wide enough input target
distribution to satisfy these criteria. If it fails to produce a given type of reduced taget
list, the script will give a warning and will not write the aborted subset.

USAGE

Two command-line arguments are required.
1. The target list produced by list_builder.py (by Joseph Long and Neil Zimmerman)
2. The availability table produced by availability_checker.py (by Joseph Long and Neil Zimmerman)

EXAMPLE

python2 sort_targets.py ../target_lists/initial_image_mosaic_R45 ../initial_image_mosaic_R45_avail.npy

'''

import numpy as np
from astropy.table import Table, vstack
import glob
import matplotlib
from matplotlib import pyplot as plt
import sys
import os
import multiprocessing
import ephem
import subprocess
import requests
from urllib import urlencode
import astropy.io.votable
import astroquery.simbad
import StringIO
import warnings

def equatorial_deg_to_ecliptic_deg(ra, dec):
    """Convert RA and declination as decimal degrees to
    ecliptic longitude and latitude in decimal degrees"""
    eq = ephem.Equatorial(ra * ephem.degree, dec * ephem.degree)
    ec = ephem.Ecliptic(eq)
    return ec.lon / ephem.degree, ec.lat / ephem.degree

def best_rand_subset(targets, daily_min, N_subset, N_tries, elat_weight_pow):
    reduc_status = False
    N_cand = len(targets)
    weights = np.power(np.abs(targets['eb'].data), elat_weight_pow)
    weights = weights / np.sum(weights)
    best_N_gaps = 365
    best_N_under = 365
    best_sum_avail = 0
    best_subset_ind = None
    avail_array = targets['avail'].data

    for t in range(N_tries):
        subset_ind = np.random.choice(N_cand, N_subset, replace=False, p=weights)
        avail_totals = np.sum(avail_array[subset_ind,:], axis=0)
        N_gaps = np.where(avail_totals == 0)[0].shape[0]
        N_under = np.where(avail_totals < daily_min)[0].shape[0]
        sum_avail = np.sum(avail_totals)
        if N_gaps < best_N_gaps: # first priority is to reduce the gaps
            best_subset_ind = subset_ind.copy()
            best_N_gaps = N_gaps
            best_N_under = N_under
            best_sum_avail = sum_avail # note, this specific metric may actually get worse
            print('Sample #%9d: gap day total %3d, underfill day total %3d, sum avail %5d' % (t+1, best_N_gaps, best_N_under, best_sum_avail))
        elif (N_under < best_N_under and N_gaps <= best_N_gaps): # second priority
            best_subset_ind = subset_ind.copy()
            best_N_gaps = N_gaps
            best_N_under = N_under
            best_sum_avail = sum_avail # note, this specific metric may actually get worse
            print('Sample #%9d: gap day total %3d, underfill day total %3d, sum avail %5d' % (t+1, best_N_gaps, best_N_under, best_sum_avail))
        elif (sum_avail > best_sum_avail and N_under <= best_N_under and N_gaps <= best_N_gaps): # third priority
            best_subset_ind = subset_ind.copy()
            best_N_gaps = N_gaps
            best_N_under = N_under
            best_sum_avail = sum_avail
            print('Sample #%9d: gap day total %3d, underfill day total %3d, sum avail %5d' % (t+1, best_N_gaps, best_N_under, best_sum_avail))

    best_subset = targets[best_subset_ind]
    best_subset.sort(keys='eb')
    best_subset.reverse()
    return best_subset, best_N_gaps, best_N_under, best_sum_avail

def best_rand_subset_mp(targets, daily_min, N_subset, N_tries, elat_weight_pow, best_N_gaps, best_N_under, best_sum_avail, best_subset_ind):
    proc_name = multiprocessing.current_process().name
    print proc_name, 'launching'
    N_cand = len(targets)
    weights = np.power(np.abs(targets['eb'].data), elat_weight_pow)
    weights = weights / np.sum(weights)
    avail_array = targets['avail'].data
    np.random.seed()

    for t in range(N_tries):
        subset_ind = np.random.choice(N_cand, N_subset, replace=False, p=weights)
        avail_totals = np.sum(avail_array[subset_ind,:], axis=0)
        N_gaps = np.where(avail_totals == 0)[0].shape[0]
        N_under = np.where(avail_totals < daily_min)[0].shape[0]
        sum_avail = np.sum(avail_totals)
        if N_gaps < best_N_gaps.value: # first priority is to reduce the gaps
            best_subset_ind[:] = subset_ind.copy()
            best_N_gaps.value = N_gaps
            best_N_under.value = N_under
            best_sum_avail.value = sum_avail # note, this specific metric may actually get worse
            print('%s sample #%9d: gap day total %3d, underfill day total %3d, sum avail %5d' %\
                  (proc_name, t+1, best_N_gaps.value, best_N_under.value, best_sum_avail.value))
        elif (N_under < best_N_under.value and N_gaps <= best_N_gaps.value): # second priority
            best_subset_ind[:] = subset_ind.copy()
            best_N_gaps.value = N_gaps
            best_N_under.value = N_under
            best_sum_avail.value = sum_avail # note, this specific metric may actually get worse
            print('%s sample #%9d: gap day total %3d, underfill day total %3d, sum avail %5d' %\
                  (proc_name, t+1, best_N_gaps.value, best_N_under.value, best_sum_avail.value))
        elif (sum_avail > best_sum_avail and N_under <= best_N_under and N_gaps <= best_N_gaps): # third priority
            best_subset_ind[:] = subset_ind.copy()
            best_N_gaps.value = N_gaps
            best_N_under.value = N_under
            best_sum_avail.value = sum_avail # note, this specific metric may actually get worse
            print('%s sample #%9d: gap day total %3d, underfill day total %3d, sum avail %5d' %\
                  (proc_name, t+1, best_N_gaps.value, best_N_under.value, best_sum_avail.value))
    print proc_name, 'finished'

def make_reduced_table(in_table, daily_min, max_length):
    reduc_status = False
    N_cand = len(in_table)
    max_length = np.min([N_cand, max_length])
    reduc_table = Table(in_table[0])
    rr = 1
    while np.sum(reduc_table['avail'],axis=0).min() < daily_min and rr < N_cand and len(reduc_table) < max_length:
        reduc_table.add_row(in_table[rr])
        rr = rr + 1
    if np.sum(reduc_table['avail'],axis=0).min() < daily_min:
        print('    WARNING: Could not reduce this list while satisfying >= %d targets per day.' % (daily_min)) 
        if len(reduc_table) < max_length:
            print('    Reached end of %d-star candidate list before filling up to max allowed list length of %d stars; incomplete result has %d stars' %\
                  (N_cand, max_length, len(reduc_table)))
        elif np.sum(reduc_table['avail'],axis=0).min() < daily_min and len(reduc_table) == max_length:
            print('    Incomplete reduced list filled up to max allowed length of %d stars.' % (max_length))
            reduc_status = True
    else:
        print('    Reduction was successful.')
        reduc_status = True
    reduc_table.sort(keys='eb')
    reduc_table.reverse()
    return reduc_table, reduc_status

targets_fname = sys.argv[-2]
avail_fname = sys.argv[-1]

min_targets_per_day = 2
min_targets_per_day_hemi = 2
max_reduc_length = 5
max_reduc_length_hemi = 5
N_rand_samp = int(1e7)
elat_weight_pow = 1.5
N_proc = 30

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

twomass_IDs = []
print("Retrieving 2MASS Point Source Catalog IDs from IPAC IRSA...")
irsa_query_vot = astropy.io.votable.tree.VOTableFile()
resource = astropy.io.votable.tree.Resource()
irsa_query_vot.resources.append(resource)
table = astropy.io.votable.tree.Table(irsa_query_vot)
resource.tables.append(table)
table.fields.extend([
        astropy.io.votable.tree.Field(irsa_query_vot, name="ra", datatype="double", arraysize="*"),
        astropy.io.votable.tree.Field(irsa_query_vot, name="dec", datatype="double", arraysize="*")])
table.create_arrays(N_targets_full)
for rr, row in enumerate(targets):
    table.array[rr] = (row['RA'], row['Dec']) 

#irsa_query_vot_fname = "irsa_query.xml"
irsa_query_vot_fname = os.path.join(os.getcwd(), "irsa_query.xml")
irsa_query_vot.to_xml(irsa_query_vot_fname)
irsa_response_vot_fname = os.path.join(os.getcwd(), "irsa_response.xml")

irsa_search_rad = 1./3600 # in degree units
curl_cmd = "curl -o \"%s\" -F \"UPLOAD=my_table,param:table\" -F \"table=@%s\" -F \"QUERY=SELECT fp_psc.designation FROM fp_psc WHERE CONTAINS(POINT(\'J2000\',ra,dec), CIRCLE(\'J2000\',TAP_UPLOAD.my_table.ra, TAP_UPLOAD.my_table.dec, %.5f))=1\" http://irsa.ipac.caltech.edu/TAP/sync" % (irsa_response_vot_fname, irsa_query_vot_fname, irsa_search_rad)
print curl_cmd
subprocess.call(curl_cmd, shell=True)

irsa_response_vot = astropy.io.votable.parse_single_table( irsa_response_vot_fname, columns=['designation'] )
assert( len(irsa_response_vot.array.data) == N_targets_full )
for rr, row in enumerate(targets):
    irsa_ID = "2MASS J%s"%(irsa_response_vot.array.data['designation'][rr])
    twomass_IDs.append(irsa_ID)

el_col = Table.Column(el, 'el')
eb_col = Table.Column(eb, 'eb')
twomass_col = Table.Column(twomass_IDs, '2MASS')
avail_col = Table.Column(avail_1yr, 'avail')
N_days_per_target = np.sum(avail_1yr, axis=1)
N_days_col = Table.Column(N_days_per_target, 'N_days')
targets.add_column(el_col, index=2)
targets.add_column(eb_col, index=3)
targets.add_column(avail_col) 
targets.add_column(N_days_col)
targets.add_column(twomass_col, index=0)

# Remove stars with multiplicity flags in SIMBAD object types.
warnings.simplefilter("ignore", UserWarning) # Don't complain about 2MASS stars that lack SIMBAD entries

customSimbad = astroquery.simbad.Simbad()
customSimbad.add_votable_fields('otype(N)')
customSimbad.add_votable_fields('id(2MASS)')
simbad_result = customSimbad.query_objects(targets['2MASS'].data)
for ss, simbad_row in enumerate(simbad_result):
    otype_part2 = int(simbad_row['OTYPE_N'].split('.')[1])
    if otype_part2 != 13:
        continue
    else:
        rr = np.where(targets['2MASS'] == simbad_row['ID_2MASS'])[0][0]
        print 'Rejecting %s, object type code %s from reduced list due to multiplicity classification in SIMBAD'%\
              (targets['2MASS'][rr], simbad_row['OTYPE_N'])
        targets.remove_row(rr)

print('%d candidates remain after removing stars with multiplicity flags.'%(len(targets)))

targets_eN = targets[targets['eb'] >= 0]
targets_eS = targets[targets['eb'] < 0]
N_targets_eN = len(targets_eN)
N_targets_eS = len(targets_eS)

dens_el_eS = [] # number of stars within a range of ecliptic longitudes
sort_metric_eS = []
for row in targets_eS:
    hw = 15.
    el_upper = (row['el'] + hw) % 360
    el_lower = (row['el'] - hw) % 360
    if el_upper > el_lower:
        N_el_neighbs = len(targets_eS[ np.all([[el_lower < targets_eS['el']],[targets_eS['el'] < el_upper]],axis=0).ravel() ])
    else:
        N_el_neighbs = len(targets_eS[ np.any([[el_lower < targets_eS['el']],[targets_eS['el'] < el_upper]],axis=0).ravel() ])
    dens_el_eS.append(N_el_neighbs)
    sort_metric_eS.append(row['N_days']/N_el_neighbs**(0.05))
el_dens_eS_col = Table.Column(dens_el_eS, 'el_dens')
sort_metric_eS_col = Table.Column(sort_metric_eS, 'sort_metric')
targets_eS.add_column(el_dens_eS_col)
targets_eS.add_column(sort_metric_eS_col)

dens_el_eN = [] # number of stars within a range of ecliptic longitudes
sort_metric_eN = []
for row in targets_eN:
    hw = 15.
    el_upper = (row['el'] + hw) % 360
    el_lower = (row['el'] - hw) % 360
    if el_upper > el_lower:
        N_el_neighbs = len(targets_eN[ np.all([[el_lower < targets_eN['el']],[targets_eN['el'] < el_upper]],axis=0).ravel() ])
    else:
        N_el_neighbs = len(targets_eN[ np.any([[el_lower < targets_eN['el']],[targets_eN['el'] < el_upper]],axis=0).ravel() ])
    dens_el_eN.append(N_el_neighbs)
    sort_metric_eN.append(row['N_days']/N_el_neighbs**(0.05))
el_dens_eN_col = Table.Column(dens_el_eN, 'el_dens')
sort_metric_eN_col = Table.Column(sort_metric_eN, 'sort_metric')
targets_eN.add_column(el_dens_eN_col)
targets_eN.add_column(sort_metric_eN_col)

targets.sort(keys='N_days')
targets.reverse()
#targets_eN.sort(keys='N_days')
targets_eN.sort(keys='sort_metric')
targets_eN.reverse()
#targets_eS.sort(keys='N_days')
targets_eS.sort(keys='sort_metric')
targets_eS.reverse()

total_avail_all = np.sum(avail_1yr, axis=0)
min_avail_all = total_avail_all.min()
min_avail_eS = np.sum(targets_eS['avail'],axis=0).min()
min_avail_eN = np.sum(targets_eN['avail'],axis=0).min()

print("On any day of the year, at least %d targets are available at all ecliptic latitudes." % min_avail_all)
print("On any day of the year, at least %d targets are available in the northern ecliptic hemisphere." % min_avail_eN)
print("On any day of the year, at least %d targets are available in the southern ecliptic hemisphere." % min_avail_eS)

print('Reducing all-latitude list...')
reduc_targets, reduc_status = make_reduced_table(targets, min_targets_per_day, max_reduc_length)
print('Reducing northern latitude list...')
reduc_targets_eN, reduc_status_eN = make_reduced_table(targets_eN, min_targets_per_day_hemi, max_reduc_length_hemi)
print('Reducing southern latitude list...')
reduc_targets_eS, reduc_status_eS = make_reduced_table(targets_eS, min_targets_per_day_hemi, max_reduc_length_hemi)
print('Number of gap days in Jan reduced southern list = %d'%(np.where( np.sum(reduc_targets_eS['avail'], axis=0) == 0 )[0].shape[0]))
print('Number of underfill days in Jan reduced southern list = %d'%(np.where( np.sum(reduc_targets_eS['avail'], axis=0) < min_targets_per_day_hemi )[0].shape[0]))

print('Searching for best subset of %d stars for southern latitude list...'%max_reduc_length_hemi)

if N_proc > 1:
    global_N_gaps = multiprocessing.Value('I', 365)
    global_N_under = multiprocessing.Value('I', 365)
    global_sum_avail = multiprocessing.Value('L', 0)
    global_subset_ind = multiprocessing.Array('i', max_reduc_length_hemi)
    # max_N_proc = multiprocessing.cpu_count()
    job_list = [ multiprocessing.Process( target=best_rand_subset_mp,\
                                          args=(targets_eS, min_targets_per_day_hemi, max_reduc_length_hemi,\
                                                N_rand_samp, elat_weight_pow, global_N_gaps, global_N_under,\
                                                global_sum_avail, global_subset_ind) ) for j in range(N_proc) ]
    for j in job_list:
        j.start()
    for j in job_list:
        j.join()
    
    best_reduc_targets_eS_mp = targets_eS[global_subset_ind[:]]
    best_reduc_targets_eS_mp.sort(keys='eb')
    best_reduc_targets_eS_mp.reverse()
    reduc_targets_eS = best_reduc_targets_eS_mp
else:
    best_reduc_targets_eS, best_N_gaps, best_N_under, best_sum_avail = best_rand_subset(targets_eS, min_targets_per_day_hemi,\
                                                                                        max_reduc_length_hemi, N_rand_samp, elat_weight_pow)
    reduc_targets_eS = best_reduc_targets_eS

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
