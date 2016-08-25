import sys
import os
import datetime
import numpy as np
import subprocess
import argparse

from list_specs import target_lists

parser = argparse.ArgumentParser(description="Top level driver script to produce lists of OTE commissioning target candidates.")
parser.add_argument("--newfilepath", type=str, default=os.path.normpath("../target_lists_{:s}".format(datetime.datetime.now().strftime("%Y-%m-%d"))),
                    help="Destination directory for new files.")
args = parser.parse_args()

subprocess.call("mkdir -p {:s}".format(args.newfilepath), shell=True)

#///////////////////////////////////////////
# BUILD TARGET LISTS
#///////////////////////////////////////////
subprocess.call(['python2', 'list_builder.py', 'early_comm', '--newfilepath', args.newfilepath])
early_comm_list_name = os.path.join(args.newfilepath, 'early_comm', 'early_comm')
assert os.path.exists(early_comm_list_name), 'Expected list_builder product does not exist.'

subprocess.call(['python2', 'list_builder.py', 'global_alignment', '--newfilepath', args.newfilepath])
global_alignment_list_name = os.path.join(args.newfilepath, 'global_alignment', 'global_alignment')
assert os.path.exists(global_alignment_list_name), 'Expected list_builder product does not exist.'

#///////////////////////////////////////////
# GET CALENDAR AVAILABILITY TABLE AND PLOT
#///////////////////////////////////////////
subprocess.call(['python2', 'availability_checker.py', early_comm_list_name, '--lite'])
early_comm_avail_name = os.path.join(args.newfilepath, 'early_comm', 'early_comm_avail.npy')
assert os.path.exists(early_comm_avail_name), 'Expected availability_checker product does not exist.'

subprocess.call(['python2', 'availability_checker.py', global_alignment_list_name, '--lite'])
global_alignment_avail_name = os.path.join(args.newfilepath, 'global_alignment', 'global_alignment_avail.npy')
assert os.path.exists(global_alignment_avail_name), 'Expected availability_checker product does not exist.'

#///////////////////////////////////////////
# SORT, SCREEN, AND REDUCE TARGET LISTS
#///////////////////////////////////////////
subprocess.call(['python2', 'sort_targets.py', early_comm_list_name, early_comm_avail_name, '--max_length', '10'])
early_comm_reduc_list_name = os.path.join(args.newfilepath, 'early_comm', 'early_comm_reduc')
early_comm_reduc_list_2mass_name = os.path.join(args.newfilepath, 'early_comm', 'early_comm_NS_reduc_2MASS')
early_comm_reduc_list_csv_name = os.path.join(args.newfilepath, 'early_comm', 'early_comm_NS_reduc_apt.csv')
assert os.path.exists(early_comm_reduc_list_name), 'Expected sort_targets product does not exist.'
assert os.path.exists(early_comm_reduc_list_2mass_name), 'Expected sort_targets product does not exist.'
assert os.path.exists(early_comm_reduc_list_csv_name), 'Expected sort_targets product does not exist.'

subprocess.call(['python2', 'sort_targets.py', global_alignment_list_name, global_alignment_avail_name, '--max_length', '10'])
global_alignment_reduc_list_name = os.path.join(args.newfilepath, 'global_alignment', 'global_alignment_reduc')
global_alignment_reduc_list_2mass_name = os.path.join(args.newfilepath, 'global_alignment', 'global_alignment_NS_reduc_2MASS')
global_alignment_reduc_list_csv_name = os.path.join(args.newfilepath, 'global_alignment', 'global_alignment_NS_reduc_apt.csv')
assert os.path.exists(global_alignment_reduc_list_name), 'Expected sort_targets product does not exist.'
assert os.path.exists(global_alignment_reduc_list_2mass_name), 'Expected sort_targets product does not exist.'
assert os.path.exists(global_alignment_reduc_list_csv_name), 'Expected sort_targets product does not exist.'

#///////////////////////////////////////////
# GET CALENDAR AVAILABILITY OF REDUCED LIST
#///////////////////////////////////////////
subprocess.call(['python2', 'availability_checker.py', early_comm_list_name, early_comm_reduc_list_name, '--lite'])

subprocess.call(['python2', 'availability_checker.py', global_alignment_list_name, global_alignment_reduc_list_name, '--lite'])

#///////////////////////////////////////////
# WRITE HTML GALLERY
#///////////////////////////////////////////
early_comm_cutout_size = 2*float(target_lists['early_comm']['neighbors'][0]['r_arcmin'])/60
subprocess.call(['python2', 'write_gallery.py', early_comm_reduc_list_2mass_name, '{:.3f}'.format(early_comm_cutout_size)])

global_alignment_cutout_size = 2*float(target_lists['global_alignment']['neighbors'][0]['r_arcmin'])/60
subprocess.call(['python2', 'write_gallery.py', global_alignment_reduc_list_2mass_name, '{:.3f}'.format(global_alignment_cutout_size)])
