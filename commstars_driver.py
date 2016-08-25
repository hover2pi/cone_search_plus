import sys
import os
import datetime
import numpy as np
import subprocess
import argparse

from list_specs import target_lists

def run_commstars_pipeline(category, new_file_path, max_final_length=10, min_per_day=3, min_per_day_hemi=3):
    #///////////////////////////////////////////
    # BUILD TARGET LISTS
    #///////////////////////////////////////////
    subprocess.call(['python2', 'list_builder.py', category, '--newfilepath', new_file_path])
    list_name = os.path.join(new_file_path, category, category)
    assert os.path.exists(list_name), 'Expected list_builder product does not exist.'

    #///////////////////////////////////////////
    # GET CALENDAR AVAILABILITY TABLE AND PLOT
    #///////////////////////////////////////////
    subprocess.call(['python2', 'availability_checker.py', list_name, '--lite'])
    avail_name = os.path.join(new_file_path, category, '{:s}_avail.npy'.format(category))
    assert os.path.exists(avail_name), 'Expected availability_checker product does not exist.'

    #///////////////////////////////////////////
    # SORT, SCREEN, AND REDUCE TARGET LISTS
    #///////////////////////////////////////////
    subprocess.call(['python2', 'sort_targets.py', list_name, avail_name,
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
    subprocess.call(['python2', 'availability_checker.py', list_name, reduc_list_name, '--lite'])
    
    #///////////////////////////////////////////
    # WRITE HTML GALLERY
    #///////////////////////////////////////////
    cutout_size = 2*float(target_lists[category]['neighbors'][0]['r_arcmin'])/60
    subprocess.call(['python2', 'write_gallery.py', reduc_list_2mass_name, '{:.3f}'.format(cutout_size)])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Top level driver script to produce lists of OTE commissioning target candidates.")
    parser.add_argument("--newfilepath", type=str, default=os.path.normpath("../target_lists_{:s}".format(datetime.datetime.now().strftime("%Y-%m-%d"))),
                        help="Destination directory for new files.")
    args = parser.parse_args()
    
    subprocess.call("mkdir -p {:s}".format(args.newfilepath), shell=True)

    #///////////////////////////////////////////
    # EARLY COMMISSIONING
    #///////////////////////////////////////////
    run_commstars_pipeline('early_comm', args.newfilepath)

    #///////////////////////////////////////////
    # GLOBAL ALIGNMENT
    #///////////////////////////////////////////
    run_commstars_pipeline('global_alignment', args.newfilepath, min_per_day_hemi=3)

    #///////////////////////////////////////////
    # COARSE PHASING
    #///////////////////////////////////////////
    run_commstars_pipeline('coarse_phasing', args.newfilepath, min_per_day_hemi=3)

    #///////////////////////////////////////////
    # FINE PHASING
    #///////////////////////////////////////////
    run_commstars_pipeline('fine_phasing', args.newfilepath, min_per_day_hemi=5)
