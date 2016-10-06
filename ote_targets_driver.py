#!/usr/bin/env python2.7

import sys
import os
import datetime
import numpy as np
import subprocess
import multiprocessing
import argparse
import logging

from list_specs import target_lists

def run_ote_target_pipeline(category, target_file_path, new_file_path, log_file_path,
                            all_cores=False, max_final_length=10, min_per_day=3, min_per_day_hemi=3):
#    logging.basicConfig(filename="/dev/null", level=logging.DEBUG)
    logging.basicConfig(level=logging.DEBUG)
#    logging.basicConfig(level=logging.ERROR)
    logger = logging.getLogger()
    fh = logging.FileHandler(log_file_path, mode='w')
    logger.addHandler(fh)

    logger.info('////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////')
    logger.info('////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////')
    logger.info('           Target category {:}\n'.format(category))

    #///////////////////////////////////////////
    # BUILD TARGET LISTS
    #///////////////////////////////////////////
    logger.info('\n////////////////////////////////////////////////////////////////////')
    logger.info('           Build target list')

    if all_cores:
        n_cores = multiprocessing.cpu_count()
    else:
        n_cores = multiprocessing.cpu_count()/2
    owd = os.getcwd()
    os.chdir(os.path.normpath(os.environ['OTE_TARGETS_PIPELINE']))
    output = subprocess.check_output(['python2.7', 'list_builder.py', category,
                                      '--targetfilepath', target_file_path,
                                      '--newfilepath', new_file_path,
                                      '--logfilepath', log_file_path,
                                      '--ncores', str(n_cores)])
    logger.info(output)
    list_name = os.path.join(new_file_path, category, category)
    assert os.path.exists(list_name), 'Expected list_builder product does not exist.'
    fh.close()
    logger.removeHandler(fh)
    fh = logging.FileHandler(log_file_path, mode='a')
    logger.addHandler(fh)

    #///////////////////////////////////////////
    # GET CALENDAR AVAILABILITY TABLE AND PLOT
    #///////////////////////////////////////////
    logger.info('\n////////////////////////////////////////////////////////////////////')
    logger.info('           Availability checker (first pass)')
    output = subprocess.call(['python2.7', 'availability_checker.py', list_name,
                              '--logfilepath', log_file_path, '--lite'])
    avail_name = os.path.join(new_file_path, category, '{:s}_avail.npy'.format(category))
    assert os.path.exists(avail_name), 'Expected availability_checker product does not exist.'

    #///////////////////////////////////////////
    # SORT, SCREEN, AND REDUCE TARGET LISTS
    #///////////////////////////////////////////
    logger.info('\n////////////////////////////////////////////////////////////////////')
    logger.info('           Sort, screen, and reduce')
    output = subprocess.call(['python2.7', 'sort_targets.py', list_name, avail_name,
                              '--logfilepath', log_file_path,
                              '--max_length', str(max_final_length),
                              '--min_per_day', str(min_per_day),
                              '--min_per_day_hemi', str(min_per_day_hemi)])
    reduc_list_name = os.path.join(new_file_path, category, '{:s}_NS_reduc'.format(category))
    reduc_list_2mass_name = os.path.join(new_file_path, category, '{:s}_NS_reduc_2MASS'.format(category))
    reduc_list_csv_name = os.path.join(new_file_path, category, '{:s}_NS_reduc_apt.csv'.format(category))
    assert os.path.exists(reduc_list_name), 'Expected sort_targets product does not exist.'
    assert os.path.exists(reduc_list_2mass_name), 'Expected sort_targets product does not exist.'
    assert os.path.exists(reduc_list_csv_name), 'Expected sort_targets product does not exist.'

    #///////////////////////////////////////////
    # GET CALENDAR AVAILABILITY OF REDUCED LIST
    #///////////////////////////////////////////
    logger.info('\n////////////////////////////////////////////////////////////////////')
    logger.info('           Availability checker (second pass)')
    output = subprocess.call(['python2.7', 'availability_checker.py', list_name, reduc_list_name,
                              '--logfilepath', log_file_path, '--lite'])
    
    #///////////////////////////////////////////
    # WRITE HTML GALLERY
    #///////////////////////////////////////////
    logger.info('\n////////////////////////////////////////////////////////////////////')
    logger.info('           Write HTML gallery')
    cutout_size = 2*float(target_lists[category]['neighbors'][0]['r_arcmin'])/60
    output = subprocess.call(['python2.7', 'write_gallery.py', reduc_list_2mass_name, 
                              '--fov', '{:.3f}'.format(cutout_size),
                              '--logfilepath', log_file_path])

    os.chdir(owd)

if __name__ == "__main__":
    assert "OTE_TARGETS" in os.environ, "Set the OTE_TARGETS shell variable to specify location of new target list products, otherwise provide an explicit path using the --newfilepath command line argument"
    assert "OTE_TARGETS_PIPELINE" in os.environ, "Set the OTE_TARGETS_PIPELINE shell variable to specify the path to the OTE target pipeline"
    parser = argparse.ArgumentParser(description="Top level driver script to produce lists of OTE commissioning target candidates.")
    parser.add_argument("--targetfilepath", type=str,
                        default=os.path.normpath(os.environ["OTE_TARGETS"]),
                        help="Top level directory for OTE target products.")
    parser.add_argument("--newfilepath", type=str,
                        default=os.path.join(os.path.normpath(os.environ["OTE_TARGETS"]),"target_lists_{:s}".format(datetime.datetime.now().strftime("%Y-%m-%d"))),
                        help="Destination directory for new files.")
    parser.add_argument("--allcores", default=False, help="Use all processor cores when running list_builder catalog queries.", action="store_true")
    args = parser.parse_args()
    
    subprocess.call("mkdir -p {:s}".format(args.targetfilepath), shell=True)
    subprocess.call("mkdir -p {:s}".format(args.newfilepath), shell=True)

    log_fname = os.path.join(args.newfilepath, "ote_targets_{:s}.log".format(datetime.datetime.now().strftime("%Y-%m-%d")))
#    logging.basicConfig(filename=log_fname, filemode='w', level=logging.DEBUG)

    #///////////////////////////////////////////
    # EARLY COMMISSIONING
    #///////////////////////////////////////////
    run_ote_target_pipeline('early_comm', args.targetfilepath, args.newfilepath, log_fname, all_cores=args.allcores)

    #///////////////////////////////////////////
    # GLOBAL ALIGNMENT
    #///////////////////////////////////////////
#    run_ote_target_pipeline('global_alignment', args.targetfilepath, args.newfilepath, log_fname, all_cores=args.allcores, min_per_day_hemi=3)

    #///////////////////////////////////////////
    # COARSE PHASING
    #///////////////////////////////////////////
#    run_ote_target_pipeline('coarse_phasing', args.targetfilepath, args.newfilepath, log_fname, all_cores=args.allcores, min_per_day_hemi=3)

    #///////////////////////////////////////////
    # FINE PHASING
    #///////////////////////////////////////////
#    run_ote_target_pipeline('fine_phasing', args.targetfilepath, args.newfilepath, all_cores=args.allcores, min_per_day_hemi=5)

    #///////////////////////////////////////////
    # ROUTINE WFS&C
    #///////////////////////////////////////////
#    run_ote_target_pipeline('routine_wfsc', args.targetfilepath, args.newfilepath, all_cores=args.allcores, 
#                            min_per_day=99999, min_per_day_hemi=99999, max_final_length=99999)
