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

    now = datetime.datetime.now()
    #twomass_list_fname = sys.argv[1]
    #twomass_list = Table.read(twomass_list_fname, format='ascii.no_header')
    twomass_list = Table.read(args.twomassIDlist_fname, format='ascii.no_header')
    twomass_IDs = twomass_list['col2']
    N_stars = len(twomass_IDs)
    logging.info('Read in the 2MASS IDs of %d stars.'%N_stars)
   
    FoV = float(args.fov) 
    #if len(sys.argv) > 2:
    #    FoV = float(sys.argv[2])
    #else:
    #    FoV = 0.5
    
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
    
    html_fname = sys.argv[1] + ".html"
    
    fname_split = sys.argv[1].split('/')
    descrip = fname_split[-1]
    
    title = '''<h2><center>%s</center><br></h2>'''%descrip
    
    html_fobj = open(html_fname, "w")
    html_fobj.write(head)
    #html_fobj.write('Hello, world')
    
    html_fobj.write('Created by {:s} on {:s}<br>'.format(getpass.getuser(), now.strftime("%Y-%m-%d %H:%M")))
    html_fobj.write('Left side: DSS; Right side: 2MASS. Cutout FoV = %d arcmin<br>'%round(FoV*60))
    
    html_fobj.write(title)
    
    for ii, star in enumerate(twomass_IDs):
        row = '''
        <div class=columns>
        <center>
        <h3>%d. <a href="http://simbad.u-strasbg.fr/simbad/sim-id?Ident=2MASS %s">2MASS %s</a></h3>
        </center>
        <div class=left id="aladin-lite-div%d-dss" style="width:300px;height:300px;"></div>
        <script type="text/javascript" src="http://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.js" charset="utf-8"></script>
        <script type="text/javascript">
        var aladin1 = A.aladin('#aladin-lite-div%d-dss', {survey: "P/DSS2/color", fov:%f, reticleSize: 1, target: "2MASS %s"});
        </script>
    
        <div class=right id="aladin-lite-div%d-2mass" style="width:300px;height:300px;"></div>
        <script type="text/javascript" src="http://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.js" charset="utf-8"></script>
        <script type="text/javascript">
        var aladin2 = A.aladin('#aladin-lite-div%d-2mass', {survey: "P/2MASS/color", fov:%f, reticleSize: 1, target: "2MASS %s"});
        </script>
        </div>
        '''%(ii+1, star.replace('+','%2B'), star, ii+1, ii+1, FoV, star, ii+1, ii+1, FoV, star)
        html_fobj.write(row)
    
    html_fobj.write(closing)
    html_fobj.close()
    
    logging.info('Wrote html cutout gallery to %s'%html_fname)
