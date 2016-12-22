#!/usr/bin/env python2.7

'''
Neil Zimmerman
January 2016

SUMMARY 

Taking a list of targets generated by list_builder.py as input,
along with separately computed availability table (produced with availability_checker.py),
form a reduced subset of favored commissioning target candidates.

DESCRIPTION

The algorithm to form this reduced subset goes as follows:

First, the stars are ranked by duration of observing window, over all ecliptic
latitudes, and then separately for northern ecliptic latitudes and souther
ecliptic latitudes, resulting in three tables.  Beginning with the star having
the longest observing window, copy stars into the reduced subset until the list
contains sufficient year-round coverage that some minimum number of targets are
available on any given calendar date.  This daily minimum number can be
specified as a command line argument, and specified separately for all ecliptic
latitudes, and for the ecliptic hemispheres:

--min_per_day=
and
--min_per_day_hemi=

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

Optional arguments:
--min_per_day=        Integer, the goal minimum number of available targets on any given calendar date
--min_per_day_hemi=   Integer, the goal minimum number of available targets per ecliptic hemisphere on any given calendar date 
--max_length=         Integer, the maximum length of the final reduced lists
--nowrite             Flag, when specified no output files are written.

EXAMPLE

python2 sort_targets.py ../target_lists/fine_phasing_cvz ../fine_phasing_cvz_avail.npy --min_per_day_hemi=4 --max_length=8 --nowrite

//////////////////////////////////////////////////////////////////////////////////////////////////////

08-23-2016 -- unified version for all target lists, parses input arguments to
determine how the list is reduced.

'''

import numpy as np
from astropy.table import Table, vstack
import glob
import matplotlib
#from matplotlib import pyplot as plt
import sys
import os
import ephem
import subprocess
import requests
import urllib
import astropy.io.votable
import datetime
#import astroquery.simbad
#import StringIO
import pdb
import argparse
import logging

def equatorial_deg_to_ecliptic_deg(ra, dec):
    """Convert RA and declination as decimal degrees to
    ecliptic longitude and latitude in decimal degrees"""
    eq = ephem.Equatorial(ra * ephem.degree, dec * ephem.degree)
    ec = ephem.Ecliptic(eq)
    return ec.lon / ephem.degree, ec.lat / ephem.degree

def make_reduced_table(in_table, daily_min, max_length):
    basesimbadurl = "http://simbad.u-strasbg.fr/simbad/sim-script?script="
    #basesimbadurl = "http://simbad.cfa.harvard.edu/simbad/sim-script?script="
    basesimbadurl = basesimbadurl + "format%20object%20%22%25IDLIST%28A;1,NAME,2MASS,HD,HIP,GSC%29%20|%25OTYPELIST%20|%20%25FLUXLIST%28K;F%29%22"
    basesimbadurl = basesimbadurl + "%0Aoutput%20console=off%20script=off%0Aquery%20id%20"
   
    requests_log = logging.getLogger("requests.packages.urllib3")
    requests_log.setLevel(logging.WARNING)

    reduc_status = False
    N_cand = len(in_table)
    rr = 0
    while rr < N_cand:
        #simbadurl = basesimbadurl + "%f"%in_table['RA'][rr]  + "%20" + "%f"%in_table['Dec'][rr] + "%20radius=" + "%ds"%matchradius
        #simbadurl = basesimbadurl + "%s"%urllib.parse.quote_plus(in_table['2MASS'][rr])
        simbadurl = basesimbadurl + "%s"%urllib.quote_plus(in_table['2MASS'][rr])
        f = requests.get(simbadurl)
        queryout = f.text
        splitlines = queryout.split('\n')
        targetline = False
        for ll, line in enumerate(splitlines):
            if '*' in line:
                targetline = line
                break
            elif 'not found' in line:
                logging.info("No SIMBAD match for %s"%(in_table['2MASS'][rr]))
                break
        if targetline:
            otype = targetline.split('|')[1]
            if '**' not in otype:
                reduc_table = Table(in_table[rr])
                rr = rr + 1
                break
            else:
                logging.info('Rejecting %s %s from reduced list due to double star classification in SIMBAD'%(in_table['2MASS'][rr],otype))
        else: # No SIMBAD entry, but assume ok
            reduc_table = Table(in_table[rr])
            rr = rr + 1
            break
        rr = rr + 1
    try:
        reduc_table
    except NameError:
        logging.warning('    No candidates in the %d-star input list qualify for the reduced list.'%(N_cand))
    else:
        while np.sum(reduc_table['avail'],axis=0).min() < daily_min and rr < N_cand and len(reduc_table) < max_length:
            simbadurl = basesimbadurl + "%s"%urllib.parse.quote_plus(in_table['2MASS'][rr])
            f = requests.get(simbadurl)
            queryout = f.text
            splitlines = queryout.split('\n')
            targetline = False
            for ll, line in enumerate(splitlines):
                if '*' in line:
                    targetline = line
                    break
                elif 'not found' in line:
                    logging.info("No SIMBAD match for %s"%(in_table['2MASS'][rr]))
                    break
            if targetline:
                otype = targetline.split('|')[1]
                if '**' not in otype:
                    reduc_table.add_row(in_table[rr])
                else:
                    logging.info('Rejecting %s %s from reduced list due to double star classification in SIMBAD'%(in_table['2MASS'][rr],otype))
            else: # No SIMBAD entry, but assume ok
                reduc_table.add_row(in_table[rr])
            rr = rr + 1
        if np.sum(reduc_table['avail'],axis=0).min() < daily_min:
            logging.info('    Could not reduce this list while satisfying >= %d targets per day.' % (daily_min))
            if len(reduc_table) < max_length:
                logging.info('    Reached end of %d-star candidate list before filling up to max allowed list length of %d stars; incomplete result has %d stars' %\
                      (N_cand, max_length, len(reduc_table)))
            elif np.sum(reduc_table['avail'],axis=0).min() < daily_min and len(reduc_table) == max_length:
                logging.info('    Incomplete reduced list filled up to max allowed length of %d stars.' % (max_length))
        else:
            logging.info('    Reduction was successful.')
        reduc_table.sort(keys='eb')
        reduc_table.reverse()
        reduc_status = True
    return reduc_table, reduc_status

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Taking a list of targets generated by list_builder.py as input, along with separately computed availability table (produced with availability_checker.py) form a reduced subset of favored commissioning target candidates.")
    parser.add_argument("--nowrite", help="Don't write any files, only display the tables.", action="store_true")
    parser.add_argument("--writeall", help="Write all variants of reduced lists, per hemisphere and full sphere.", action="store_true")
    parser.add_argument("targets", help="Target list file, produced by list_builder script.")
    parser.add_argument("availability", help="Availability table, produced by availability_checker script.")
    parser.add_argument("--logfilepath", type=str, default=None, help="Log file name")
    parser.add_argument("--max_length", type=int, help="Maximum length of reduced lists", default=10)
    parser.add_argument("--min_per_day", type=int, help="Goal for minimum number of observable stars on any given calendar day", default=3)
    parser.add_argument("--min_per_day_hemi", type=int, help="Goal for minimum number of observable stars per ecliptic hemisphere on any given calendar day", default=3)

    args = parser.parse_args()

    if args.logfilepath is None:
        log_fname = os.path.join( os.path.abspath(os.path.join(os.path.dirname(args.targets), '..')),
                                  "ote_targets_{:s}.log".format(datetime.datetime.now().strftime("%Y-%m-%d")) )
    else:
        log_fname = args.logfilepath
    #logging.basicConfig(filename=log_fname, level=logging.DEBUG)
    logger = logging.getLogger()
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    logger.addHandler(ch)
    
    targets_fname = args.targets
    avail_fname = args.availability
    min_targets_per_day = args.min_per_day
    min_targets_per_day_hemi = args.min_per_day_hemi
    max_reduc_length = args.max_length
    
    #///// EARLY COMMISSIONING
    #min_targets_per_day = 2
    #min_targets_per_day_hemi = 2
    #max_reduc_length = 10
    
    #///// GLOBAL ALIGNMENT
    #min_targets_per_day = 21
    #min_targets_per_day_hemi = 3
    #max_reduc_length = 21
    
    #///// COARSE PHASING
    #min_targets_per_day = 2
    #min_targets_per_day_hemi = 2
    #max_reduc_length = 15
    
    #///// FINE PHASING
    #min_targets_per_day = 5
    #min_targets_per_day_hemi = 5
    #max_reduc_length = 10
    
    use_sort_metric = False
    
    targets = Table.read(targets_fname, format='ascii', names=['RA', 'Dec', 'J', 'H', 'K', 'qual', 'idx'])
    
    avail = np.load(avail_fname)
    avail_1yr = avail[:,:365]
    avail_col = Table.Column(avail_1yr, 'avail')
    N_days_per_target = np.sum(avail_1yr, axis=1)
    N_days_col = Table.Column(N_days_per_target, 'N_days')
    targets.add_column(avail_col)
    targets.add_column(N_days_col)
    
    targets.sort(keys='RA')
    N_targets_full = len(targets)
    
    assert N_targets_full == avail.shape[0]
    
    el = []
    eb = []
    for row in targets:
        l, b = equatorial_deg_to_ecliptic_deg(row['RA'], row['Dec'])
        el.append(l)
        eb.append(b)
    
    twomass_IDs = []
    logging.info("Retrieving 2MASS Point Source Catalog IDs from IPAC IRSA...")
    irsa_query_vot = astropy.io.votable.tree.VOTableFile()
    resource = astropy.io.votable.tree.Resource()
    irsa_query_vot.resources.append(resource)
    table = astropy.io.votable.tree.Table(irsa_query_vot)
    resource.tables.append(table)
    table.fields.extend([ astropy.io.votable.tree.Field(irsa_query_vot, name="ra", datatype="double"),
                          astropy.io.votable.tree.Field(irsa_query_vot, name="dec", datatype="double") ])
    table.create_arrays(N_targets_full)
    for rr, row in enumerate(targets):
        table.array[rr] = (row['RA'], row['Dec'])
    
    #irsa_query_vot_fname = "irsa_query.xml"
    irsa_query_vot_fname = os.path.join(os.getcwd(), "irsa_query.xml")
    irsa_query_vot.to_xml(irsa_query_vot_fname)
    irsa_response_vot_fname = os.path.join(os.getcwd(), "irsa_response.xml")
    
    irsa_search_rad = 1./3600 # in degree units
    #curl_cmd = "curl -o \"%s\" -F \"UPLOAD=my_table,param:table\" -F \"table=@%s\" -F \"QUERY=SELECT fp_psc.designation FROM fp_psc WHERE CONTAINS(POINT(\'J2000\',ra,dec), CIRCLE(\'J2000\',TAP_UPLOAD.my_table.ra, TAP_UPLOAD.my_table.dec, %.5f))=1\" http://irsa.ipac.caltech.edu/TAP/sync" % (irsa_response_vot_fname, irsa_query_vot_fname, irsa_search_rad)
    curl_cmd = "curl -o \"%s\" -F \"UPLOAD=my_table,param:table\" -F \"table=@%s\" -F \"QUERY=SELECT fp_psc.designation,fp_psc.ra,fp_psc.dec FROM fp_psc WHERE CONTAINS(POINT(\'J2000\',ra,dec), CIRCLE(\'J2000\',TAP_UPLOAD.my_table.ra, TAP_UPLOAD.my_table.dec, %.5f))=1\" http://irsa.ipac.caltech.edu/TAP/sync" % (irsa_response_vot_fname, irsa_query_vot_fname, irsa_search_rad)
    #curl_cmd = "curl -o \"%s\" -F \"UPLOAD=my_table,param:table\" -F \"table=@%s\" -F \"QUERY=SELECT fp_psc.designation,fp_psc.ra,fp_psc.dec FROM fp_psc WHERE CONTAINS(POINT(\'J2000\',ra,dec), CIRCLE(\'J2000\',TAP_UPLOAD.my_table.ra, TAP_UPLOAD.my_table.dec, %.5f))=1\" http://irsa.ipac.caltech.edu/TAP/sync >> %s" % (irsa_response_vot_fname, irsa_query_vot_fname, irsa_search_rad, log_fname)
    logging.info(curl_cmd)
    subprocess.call(curl_cmd, shell=True)
    
    #irsa_response_vot = astropy.io.votable.parse_single_table( irsa_response_vot_fname, columns=['designation'] )
    irsa_response_vot = astropy.io.votable.parse_single_table( irsa_response_vot_fname, columns=['designation','ra','dec'] )
    assert( len(irsa_response_vot.array.data) == N_targets_full )
    irsa_response_vot.array.data.sort(order='ra')
    for rr, row in enumerate(targets):
        irsa_ID = "2MASS J%s"%(irsa_response_vot.array.data['designation'][rr])
        twomass_IDs.append(irsa_ID)
    
    el_col = Table.Column(el, 'el')
    eb_col = Table.Column(eb, 'eb')
    twomass_col = Table.Column(twomass_IDs, '2MASS')
    targets.add_column(el_col, index=2)
    targets.add_column(eb_col, index=3)
    targets.add_column(twomass_col, index=0)
    
    #pdb.set_trace()
    
    #black_list = ['2MASS J21093180+6829268']
    black_list = []
    for rr, row in enumerate(targets):
        if row['2MASS'] in black_list:
            targets.remove_row(rr)
            logging.info("Removed externally flagged star %s from candidate list"%(row['2MASS']))
    
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
    if use_sort_metric:
        targets_eN.sort(keys='sort_metric')
        targets_eS.sort(keys='sort_metric')
    else:
        targets_eN.sort(keys='eb')
        targets_eS.sort(keys='eb')
    targets_eN.reverse()
    #targets_eS.reverse()
    
    total_avail_all = np.sum(avail_1yr, axis=0)
    min_avail_all = total_avail_all.min()
    min_avail_eS = np.sum(targets_eS['avail'],axis=0).min()
    min_avail_eN = np.sum(targets_eN['avail'],axis=0).min()
    
    logging.info("Starting from a list of %d candidate stars." % N_targets_full)
    logging.info("On any day of the year, at least %d candidates are available at all ecliptic latitudes." % min_avail_all)
    logging.info("On any day of the year, at least %d candidates are available in the northern ecliptic hemisphere." % min_avail_eN)
    logging.info("On any day of the year, at least %d candidates are available in the southern ecliptic hemisphere." % min_avail_eS)
    
    base_simbad_url = "http://simbad.u-strasbg.fr/simbad/sim-script?script="
    #base_simbad_url = "http://simbad.cfa.harvard.edu/simbad/sim-script?script="
    
    logging.info('Reducing all-latitude list...')
    reduc_targets, reduc_status = make_reduced_table(targets, min_targets_per_day, max_reduc_length)
    min_avail_reduc_all = np.sum(reduc_targets['avail'],axis=0).min()
    logging.info("In the reduced all-latitude list, on any day of the year, at least %d candidates are available." % min_avail_reduc_all) 
    
    if N_targets_eN > 0:
        logging.info('Reducing northern latitude list...')
        reduc_targets_eN, reduc_status_eN = make_reduced_table(targets_eN, min_targets_per_day_hemi, max_reduc_length/2)
        min_avail_reduc_eN = np.sum(reduc_targets_eN['avail'],axis=0).min()
        logging.info("In the reduced northern latitude list, on any day of the year, at least %d candidates are available." % min_avail_reduc_eN) 
    else:
        reduc_targets_eN = None
        reduc_status_eN = False
    
    if N_targets_eS > 0:
        logging.info('Reducing southern latitude list...')
        reduc_targets_eS, reduc_status_eS = make_reduced_table(targets_eS, min_targets_per_day_hemi, max_reduc_length/2)
        min_avail_reduc_eS = np.sum(reduc_targets_eS['avail'],axis=0).min()
        logging.info("In the reduced southern latitude list, on any day of the year, at least %d candidates are available." % min_avail_reduc_eS) 
    else:
        reduc_targets_eS = None
        reduc_status_eS = False
    
    if reduc_status:
        logging.info('\nReduced list, all latitudes (%d stars for min daily avail. %d stars)'%(len(reduc_targets),min_targets_per_day))
        logging.info(reduc_targets)
        reduc_targets_full = reduc_targets['RA', 'Dec', 'J', 'H', 'K', 'qual', 'idx']
        reduc_targets_RADec = reduc_targets['RA', 'Dec']
        reduc_targets_twomass = reduc_targets['2MASS']
        reduc_targets_apt = reduc_targets['2MASS', 'RA', 'Dec', 'J', 'H', 'K']
        reduc_targets_RADec.meta = {} 
        reduc_targets_twomass.meta = {} 
        reduc_targets_apt.meta = {} 
        reduc_targets_fname = targets_fname + '_reduc'
        reduc_targets_RADec_fname = reduc_targets_fname + '_RADec'
        reduc_targets_twomass_fname = reduc_targets_fname + '_2MASS'
        reduc_targets_apt_fname = reduc_targets_fname + '_apt.csv'
        if not args.nowrite and args.writeall:
            #reduc_targets_full.write(reduc_targets_fname, format='ascii.no_header', overwrite=True)
            #reduc_targets_RADec.write(reduc_targets_RADec_fname, format='ascii.no_header', overwrite=True)
            #reduc_targets_apt.write(reduc_targets_apt_fname, format='ascii', delimiter=',', overwrite=True)
            reduc_targets_full.write(reduc_targets_fname, format='ascii.no_header')
            reduc_targets_RADec.write(reduc_targets_RADec_fname, format='ascii.no_header')
            reduc_targets_apt.write(reduc_targets_apt_fname, format='ascii', delimiter=',')
            reduc_targets_twomass.tofile(reduc_targets_twomass_fname,sep="\n")
            logging.info('Wrote all-latitude reduced target lists to\n%s,\n%s,\n%s,\n%s'%(reduc_targets_fname, reduc_targets_RADec_fname, reduc_targets_twomass_fname, reduc_targets_apt_fname))
    if reduc_status_eN:
        logging.info('\nReduced list, northern latitudes (%d stars for min daily avail. %d stars per hem.)'%(len(reduc_targets_eN),min_targets_per_day_hemi))
        logging.info(reduc_targets_eN)
        reduc_targets_eN_full = reduc_targets_eN['RA', 'Dec', 'J', 'H', 'K', 'qual', 'idx']
        reduc_targets_eN_RADec = reduc_targets_eN['RA', 'Dec']
        reduc_targets_eN_twomass = reduc_targets_eN['2MASS']
        reduc_targets_eN_RADec.meta = {}
        reduc_targets_eN_twomass.meta = {}
        reduc_targets_eN_fname = targets_fname + '_north_reduc'
        reduc_targets_eN_RADec_fname = reduc_targets_eN_fname + '_RADec'
        reduc_targets_eN_twomass_fname = reduc_targets_eN_fname + '_2MASS'
        if not args.nowrite and args.writeall:
            #reduc_targets_eN_full.write(reduc_targets_eN_fname, format='ascii.no_header', overwrite=True)
            #reduc_targets_eN_RADec.write(reduc_targets_eN_RADec_fname, format='ascii.no_header', overwrite=True)
            reduc_targets_eN_full.write(reduc_targets_eN_fname, format='ascii.no_header')
            reduc_targets_eN_RADec.write(reduc_targets_eN_RADec_fname, format='ascii.no_header')
            reduc_targets_eN_twomass.tofile(reduc_targets_eN_twomass_fname,sep="\n")
            logging.info('Wrote northern latitude reduced target lists to\n%s,\n%s,\n%s'%(reduc_targets_eN_fname, reduc_targets_eN_RADec_fname, reduc_targets_eN_twomass_fname))
    if reduc_status_eS:
        logging.info('\nReduced list, southern latitudes (%d stars for min daily avail. %d stars per hem.)'%(len(reduc_targets_eS),min_targets_per_day_hemi))
        logging.info(reduc_targets_eS)
        reduc_targets_eS_full = reduc_targets_eS['RA', 'Dec', 'J', 'H', 'K', 'qual', 'idx']
        reduc_targets_eS_RADec = reduc_targets_eS['RA', 'Dec']
        reduc_targets_eS_twomass = reduc_targets_eS['2MASS']
        reduc_targets_eS_RADec.meta = {} 
        reduc_targets_eS_twomass.meta = {} 
        reduc_targets_eS_fname = targets_fname + '_south_reduc'
        reduc_targets_eS_RADec_fname = reduc_targets_eS_fname + '_RADec'
        reduc_targets_eS_twomass_fname = reduc_targets_eS_fname + '_2MASS'
        if not args.nowrite and args.writeall:
            #reduc_targets_eS_full.write(reduc_targets_eS_fname, format='ascii.no_header', overwrite=True)
            #reduc_targets_eS_RADec.write(reduc_targets_eS_RADec_fname, format='ascii.no_header', overwrite=True)
            reduc_targets_eS_full.write(reduc_targets_eS_fname, format='ascii.no_header')
            reduc_targets_eS_RADec.write(reduc_targets_eS_RADec_fname, format='ascii.no_header')
            reduc_targets_eS_twomass.tofile(reduc_targets_eS_twomass_fname,sep="\n")
            logging.info('Wrote southern latitude reduced target lists to\n%s,\n%s,\n%s'%(reduc_targets_eS_fname, reduc_targets_eS_RADec_fname, reduc_targets_eS_twomass_fname))
    if reduc_status_eN and reduc_status_eS:
        reduc_targets_eNS = vstack([reduc_targets_eN, reduc_targets_eS])
    elif reduc_status_eN:
        reduc_targets_eNS = reduc_targets_eN
    elif reduc_status_eS:
        reduc_targets_eNS = reduc_targets_eS
    else:
        logging.warning('\nNo targets survived to reduced list; exiting with no output products')
        sys.exit(1)

    logging.info('\nReduced list, combined northern+southern (%d stars for min daily avail. %d stars per hem.)'%(len(reduc_targets_eNS),min_targets_per_day_hemi))
    logging.info(reduc_targets_eNS)
    reduc_targets_eNS_full = reduc_targets_eNS['RA', 'Dec', 'J', 'H', 'K', 'qual', 'idx']
    reduc_targets_eNS_RADec = reduc_targets_eNS['RA', 'Dec']
    reduc_targets_eNS_twomass = reduc_targets_eNS['2MASS']
    reduc_targets_eNS_apt = reduc_targets_eNS['2MASS', 'RA', 'Dec', 'J', 'H', 'K']
    reduc_targets_eNS_RADec.meta = {}
    reduc_targets_eNS_twomass.meta = {} 
    reduc_targets_eNS_apt.meta = {} 
    reduc_targets_eNS_fname = targets_fname + '_NS_reduc'
    reduc_targets_eNS_RADec_fname = reduc_targets_eNS_fname + '_RADec'
    reduc_targets_eNS_twomass_fname = reduc_targets_eNS_fname + '_2MASS'
    reduc_targets_eNS_apt_fname = reduc_targets_eNS_fname + '_apt.csv'
    if not args.nowrite:
        #reduc_targets_eNS_full.write(reduc_targets_eNS_fname, format='ascii.no_header', overwrite=True)
        #reduc_targets_eNS_RADec.write(reduc_targets_eNS_RADec_fname, format='ascii.no_header', overwrite=True)
        #reduc_targets_eNS_apt.write(reduc_targets_eNS_apt_fname, format='ascii', delimiter=',', overwrite=True)
        reduc_targets_eNS_full.write(reduc_targets_eNS_fname, format='ascii.no_header')
        reduc_targets_eNS_RADec.write(reduc_targets_eNS_RADec_fname, format='ascii.no_header')
        reduc_targets_eNS_apt.write(reduc_targets_eNS_apt_fname, format='ascii', delimiter=',')
        with open(reduc_targets_eNS_twomass_fname, mode='wt') as twomass_list_file:
            twomass_list_file.write('\n'.join(reduc_targets_eNS_twomass))
            twomass_list_file.write('\n')
        logging.info('Wrote combined northern+southern reduced target lists to\n%s,\n%s,\n%s'%(reduc_targets_eNS_fname, reduc_targets_eNS_RADec_fname, reduc_targets_eNS_twomass_fname))
