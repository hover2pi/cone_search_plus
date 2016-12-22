#!/usr/bin/env python2.7

'''
Neil Zimmerman
January 2016

SUMMARY

Taking a list of 2MASS IDs as an input argument,
generate an interactive html gallery of sky survey
cutout images from DSS and 2MASS.
You can use the 2MASS ID list written by sort_targets.py as the input file.

EXAMPLE USAGE

$ python2 write_gallery.py ../target_lists/initial_image_mosaic_R45_NS_reduc_2MASS

'''

import sys
import os
from astropy.table import Table
import datetime
import getpass
import argparse
import logging

def gallery(twomass_IDs, html_fname, FoV=0.5):
    """
    Generates a gallery of postage stamps from the 2MASS Point Source Catalog given a list of 2MASS IDS
    """
    head = '''<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
    <html>
    <head>
      <meta content="text/html; charset=ISO-8859-1"
     http-equiv="content-type">
      <title>Comm target gallery</title>
      <!-- include Aladin Lite CSS file in the head section of your page -->
      <link rel="stylesheet" href="http://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.css" />
      <!-- you can skip the following line if your page already integrates the jQuery library -->
      <script type="text/javascript" src="http://code.jquery.com/jquery-1.9.1.min.js" charset="utf-8"></script>
    
      <style>
      div.columns       { width: 700px; height: 350px; margin-left: auto; margin-right: auto;}
      div.left          { float: left; text-align:right }
      div.right         { float: right; text-align:right }
      div.caption       { width: 700px; height: 20px; margin-left: auto; margin-right: auto;}
      div.caption div   { width: 300px; float: left;}
      div.clear         { clear: both; }
      </style>
    
    </head>
    <body>\n'''
    
    closing = '''\n</body>
    </html>'''
    
    fname_split = html_fname.split('/')
    descrip = fname_split[-1]
    
    title = '''<h2><center>%s</center><br></h2>'''%descrip
    
    html_fobj = open(html_fname, "w")
    html_fobj.write(head)
    #html_fobj.write('Hello, world')
    
    now = datetime.datetime.now()
    html_fobj.write('Created by {:s} on {:s}<br>'.format(getpass.getuser(), now.strftime("%Y-%m-%d %H:%M")))
    html_fobj.write('Left side: DSS; Right side: 2MASS. Cutout FoV = {} arcmin<br>'.format(round(FoV*60)))
    
    html_fobj.write(title)

    FoV = str(FoV)
    for ii, star in enumerate(twomass_IDs):
        nm = 'J'+star.replace('+','%2B')
        idx = str(ii+1)
        row = ("<div class=columns>\n<center>\n<h3>",
               idx,
               " <a href='http://simbad.u-strasbg.fr/simbad/sim-id?Ident=2MASS ",
               nm,
               "'>2MASS J",
               star,
               "</a></h3>\n</center>\n<div class=left id='aladin-lite-div",
               idx,
               "-dss' style='width:300px;height:300px;'></div>\n",
               "<script type='text/javascript' src='http://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.js' ",
               "charset='utf-8'></script>\n<script type='text/javascript'>\nvar aladin1 = A.aladin('#aladin-lite-div",
               idx,
               "-dss', {survey: 'P/DSS2/color', fov:",
               FoV,
               " reticleSize: 1, target: '2MASS J",
               star,
               "'});\n</script>\n<div class=right id='aladin-lite-div",
               idx,
               "-2mass' style='width:300px;height:300px;'></div>\n",
               "<script type='text/javascript' src='http://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.js' charset='utf-8'></script>\n",
               "<script type='text/javascript'>\nvar aladin2 = A.aladin('#aladin-lite-div",
               idx,
               "-2mass', {survey: 'P/2MASS/color', fov:",
               FoV,
               " reticleSize: 1, target: '2MASS J",
               star,
               "'});\n</script>\n</div>")
        
        html_fobj.write(''.join(row))
    
    html_fobj.write(closing)
    html_fobj.close()
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Taking a list of 2MASS IDs, generate an interactive html gallery of sky survey cutout images from DSS and 2MASS.")
    parser.add_argument("twomassIDlist_fname", help="file name containing list of 2MASS IDs, one per line")
    parser.add_argument("--fov", help="field of view of sky survey cutout frame in degrees", default='0.5')
    parser.add_argument("--logfilepath", type=str, default=None, help="Log file name")
    args = parser.parse_args()

    if args.logfilepath is None:
        log_fname = os.path.join( os.path.abspath(os.path.join(os.path.dirname(args.twomassIDlist_fname), '..')),
                                  "ote_targets_{:s}.log".format(datetime.datetime.now().strftime("%Y-%m-%d")) )
    else:
        log_fname = args.logfilepath
    logging.basicConfig(filename=log_fname, level=logging.DEBUG, filemode='a')
    logger = logging.getLogger()
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    logger.addHandler(ch)

    #twomass_list_fname = sys.argv[1]
    #twomass_list = Table.read(twomass_list_fname, format='ascii.no_header')
    twomass_list = Table.read(args.twomassIDlist_fname, format='ascii.no_header')
    twomass_IDs = twomass_list['col2']
    N_stars = len(twomass_IDs)
    logging.info('Read in the 2MASS IDs of %d stars.'%N_stars)
   
    #if len(sys.argv) > 2:
    #    FoV = float(sys.argv[2])
    #else:
    #    FoV = 0.5
    
    gallery(twomass_IDs, sys.argv[1], FoV=float(args.fov))
    
    logging.info('Wrote html cutout gallery to %s'%html_fname)
