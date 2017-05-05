#!/usr/bin/env python
from __future__ import print_function, division
import math
import numpy as np
import astropy.table as at
import astropy.units as q
import astropy.coordinates as coords
from astroquery.irsa import Irsa
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
Vizier.ROW_LIMIT = -1

class SourceList(object):
    
    def __init__(self, source, search_radius):
        
        self.target = ''
        self.all_sources = []
        self.sources = []
        self.JH = ''
        self.HK = ''
        
        # Convert RA and Dec into string
        if isinstance(source, (tuple,list)):
            source = [i.value if hasattr(i,'unit') else i for i in source]
            source = ' '.join(map(str,source))
            
        if isinstance(source, (str,coords.builtin_frames.icrs.ICRS)):
            self.target = Simbad.query_region(source, radius=2*q.arcsec)[0]
            
        self.coords = ' '.join(self.target['RA','DEC'])
        
        # Query 2MASS for sources
        if self.target:
            self.all_sources = Vizier.query_region(self.coords, radius=search_radius, catalog=['II/246/out'])[0]
            
        print(len(self.all_sources),'sources found within',str(search_radius))
        
        self.sources = self.all_sources

    def color_cut(self, JH='', HK=''):
        """
        Reduce the list to sources within the given color cuts in J-H and H-K
        """
        # Get total stars
        start = len(self.sources)

        # Add columns for colors and color errors
        self.sources['JHmag'] = self.sources['Jmag']-self.sources['Hmag']
        self.sources['e_JHmag'] = np.sqrt(self.sources['e_Jmag']**2 + self.sources['e_Hmag']**2)
        self.sources['HKmag'] = self.sources['Hmag']-self.sources['Kmag']
        self.sources['e_HKmag'] = np.sqrt(self.sources['e_Hmag']**2 + self.sources['e_Kmag']**2)

        # Apply J-H and H-K color cuts (default values for K0-K5 giants from Straizys+2009)
        if JH:
            self.sources = self.sources[(self.sources['JHmag']+self.sources['e_JHmag']>min(JH))
                                        &(self.sources['JHmag']-self.sources['e_JHmag']<max(JH))]
            start = len(self.sources)
            print('{}/{} stars with {} < J-H > {}'.format(len(self.sources),start,min(JH),max(JH)))
        if HK:
            self.sources = self.sources[(self.sources['HKmag']+self.sources['e_HKmag']>min(JH))
                                        &(self.sources['HKmag']-self.sources['e_HKmag']<max(JH))]
            print('{}/{} stars with {} < H-K > {}'.format(len(self.sources),start,min(HK),max(HK)))
        
        self.JH = JH
        self.HK = HK
    
    def quality_cut(self, quality='A'):
        """
        Reduce the list to sources with any PSC qual flags better than or equal to the value 
        given by **quality** argument, i.e. keep only those with good photometry in JHK.
        For example, if quality='B' then acceptable PSC quality flags include
        'AAA', 'BAA', 'ABA', 'AAB', 'BBA', 'ABB', 'BAB', and 'BBB'
        """
        start = len(self.sources)
        
        keep = np.zeros(start)
        for n in range(start):
            if all([str(self.sources[n]['Qflg'])[i]<=quality for i in [2,3,4]]):
                keep[n] = 1
                
        self.sources = self.sources[np.where(keep)]
        self.quality = quality
        
        print('{}/{} stars with quality flags {} or better'.format(len(self.sources),start,quality))
        
    def proximity_cut(self, radius=8*q.arcsec, delta_K=5.):
        """
        Reduce the list to all sources with no other sources within *radius* and *delta_K*
        """
        start = len(self.sources)
        R = radius.to(q.deg).value
        
        # Iterate through the coordinates and keep only those with no neighbors
        drop = np.zeros(start)
        radec = np.array(self.sources[['RAJ2000','DEJ2000']])
        Kmags = np.array(self.sources['Kmag'])
        Kerrs = np.array(self.sources['e_Kmag'])
        for i,(C,K,E) in enumerate(zip(radec,Kmags,Kerrs)):
            for j,(c,k,e) in enumerate(zip(radec,Kmags,Kerrs)):
                if i!=j:
                    
                    # If the second source is within the radius
                    # and suficiently bright, drop it
                    if distance(C,c)<=R and abs(k-K)>=delta_K:
                        drop[i] = 1
                        break
                        
        self.sources = self.sources[np.where(1-drop)]
        
        print('{}/{} stars with neighbors within {} and {} mags'.format(len(self.sources),start,str(radius),delta_K))
    
    def binary_cut(self):
        """
        Reduce the list to singletons
        """
        pass
    
    def spectral_type_cut(self, SpT_range=('F0','M0')):
        """
        Reduce the list to sources in given spectral type range
        """
        pass
    
    def reset(self):
        self.sources = self.all_sources
        print('Resetting sources to original',len(self.sources))
        
    def show(self):
        cols = ['RAJ2000', 'DEJ2000', 'Jmag', 'e_Jmag', 'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag', 'Qflg', '_2MASS']
        self.sources[cols].pprint(max_width=-1, max_lines=-1)
        
def distance(point_a, point_b):    
    """
    ((long, lat) in deg, (long, lat) in deg, radius in deg) -> True or False

    Calculate the central angle between the two points, as computed by the haversine
    formula

    https://en.wikipedia.org/wiki/Haversine_formula
    """
    long1, lat1 = point_a
    long2, lat2 = point_b

    lambda1, phi1 = long1*math.pi/180.0, lat1*math.pi/180.0
    lambda2, phi2 = long2*math.pi/180.0, lat2*math.pi/180.0

    def haversine(theta):
        return math.sin(theta/2.0)**2

    hav_d_over_r = haversine(phi2-phi1)+math.cos(phi1)*math.cos(phi2)*haversine(lambda2-lambda1)

    central_angle_rad = 2*math.asin(math.sqrt(hav_d_over_r))
    central_angle_deg = central_angle_rad*180.0/math.pi
    
    return central_angle_deg
