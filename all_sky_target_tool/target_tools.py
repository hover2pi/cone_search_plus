#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import sys
import os
import datetime
import glob
import warnings
import numpy as np
import astropy.table as at
import astropy.units as q
import astropy.io.ascii as ii
import astropy.coordinates as coords
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pkg_resources
from astroquery.irsa import Irsa
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from IPython.core.display import display, HTML

warnings.simplefilter('ignore', UserWarning)

def vet_list(candidates, write_to=''):
    """
    Vets a list of candidates by checking catalog magnitudes, binarity, spectral type, distance, and angular diameter.
    
    Parameters
    ----------
    candidates: astropy.table.Table
        The table of candidates to vet
    write_to: str (optional)
        The file to write the resulting table to
        
    """
    # Iterate through list
    IDs, idx, SPT, PLX, W1list, W2list, kw1, w1w2 = [], [], [], [], [], [], [], []
    for n,row in enumerate(candidates):
        r, d = row['RA'], row['Dec']
        radec = coords.ICRS(ra=r*q.deg, dec=d*q.deg)
        
        rj, rh, rk = row['J'], row['H'], row['K']
        source = Vizier.query_region(radec, radius=10*q.arcsec, catalog=['II/246/out'])
        simb = Simbad.query_region(radec, radius=10*q.marcsec)
        try:
            name = simb[0]['MAIN_ID'].decode("utf-8")
        except:
            name = source[0]['_2MASS'][0].decode("utf-8")

        # Checking WDS, Speckle, Fourth Catalog of Interferometric Measurements of Binary Stars
        wds = Vizier.query_region(radec, radius=2.0*q.arcsec, catalog=['B/wds/wds','J/AJ/143/124','J/AJ/119/3084/table2'])
        if not wds:

            # Check WISE colors
            wise = Vizier.query_region(radec, radius=2.0*q.arcsec, catalog=['II/328/allwise'])
            if wise:
                j, h, k, w1, w2 = wise[0][['Jmag', 'Hmag', 'Kmag', 'W1mag', 'W2mag']][0]
                kw1.append(k-w1)
                w1w2.append(w1-w2)

                # Make sure it's a good 2MASS/WISE cross-match by checking JHK mags
                if  (j>rj-0.1) and (j<rj+0.1) \
                and (h>rh-0.1) and (h<rh+0.1) \
                and (k>rk-0.1) and (k<rk+0.1):
                    pass
                else:
                    print('Bad cross-match?',[j,h,k],[rj,rh,rk])


                #if (k-w1<=K-W1+delta) and (k-w1>=K-W1-delta) \
                #and (w1-w2<=W1-W2+delta) and (w1-w2>=W1-W2-delta):
                if True:
                    W1list.append(w1)
                    W2list.append(w2)

                    # Checking [Skiff, SEGUE, Michigan Spectral Survey]
                    catalogs = ['B/mk/mktypes', 'J/ApJ/816/80/table3', 'III/214/vol5', 'III/133/vol4', 
                                'III/76A/ls', 'III/230/catalog', 'J/ApJS/151/387/library', 'V/137D/XHIP']
                    MK = Vizier.query_region(radec, radius=8.0*q.arcsec, catalog=catalogs)
                    spts = []
                    if MK:
                        for cat in MK:
                            spts += [i.decode("utf-8") for i in cat['SpType']]
                    SPT.append(', '.join(list(set(spts))))

                    # Check Gaia for parallax
                    plx = Vizier.query_region(radec, radius=2.0*q.arcsec, 
                                              catalog=['I/337/tgas','I/311/hip2'])
                    if plx:
                        PLX.append(round(float(plx[0]['Plx']),3))
                    else:
                        PLX.append(np.nan)

                    idx.append(n)
                    IDs.append(name)
                else:
                    print(k-w1,K-W1,w1-w2,W1-W2)
                    print('Bad colors! Skipping...')
            else:
                print('Not in WISE! Skipping...')
        else:
            print('Binary at {}, {}. Skipping...'.format(r,d))

    # Add the catalog values to the final list
    final = candidates[idx]
    final['RA'].unit = final['Dec'].unit = q.deg
    final['SpT'] = SPT
    final.rename_column('arcmin', 'sep')
    final['sep'].unit = q.arcmin
    final['d'] = pi2pc(PLX)
    spts = [i.split(',')[0] for i in SPT]
    final['theta_D'] = theta_D(final['d'], spt=spts)
    final['W1'] = W1list
    final['W2'] = W2list
    final['2MASS'] = IDs
    
    simb = []
    for l in final['2MASS']:
        try:
            nm = str(Simbad.query_object("2MASS J{}".format(l))['MAIN_ID'][0]).replace("b'",'').replace("'",'')
        except:
            nm = '-'
        simb.append(nm)
    final.add_column(at.Column(simb, name='SIMBAD'), index=0)
    
    pfinal = final[[k for k in final.colnames if k not in ['idx','delta_RA','delta_Dec','J-H','H-K']]]

    # Print the table
    print('\nDone!',len(final),'candidate{} found.\n'.format('' if len(final)==1 else 's'))
    pfinal.pprint(max_width=120)
    
    # SIMBAD links
    html = []
    for row in pfinal:
        try:
            nm = row['SIMBAD'] if row['SIMBAD']!='-' else row['2MASS']
            link = 'http://simbad.u-strasbg.fr/simbad/sim-id?Ident=2MASS%20J{}'.format(row['2MASS'])
            html.append("<a href='{}' target='_'>{}</a>".format(link,nm))
        except:
            html.append("-")
    display(HTML(', '.join(html)))
    
    # Write the final list to file
    if final and write_to:
        final.meta = None
        final.write(write_to, format='ascii.fixed_width')
        print('\nTable written to',write_to)
    else:
        print('\nNo file written.')
        
    return pfinal

def polynomial_fit(n, m, degree=1, sig='', x='x', y='y', title='', c='k', ls='--', lw=2, legend=True, ax='',
                      output_data=False, dictionary=True, plot_rms=True, plot=True, verbose=False):
    """
    Fits a polynomial to the given data
    
    Parameters
    ----------
    n: array-like
        The x-axis values to fit
    m: array-like
        The y-axis values to fit
    degree: int
        The degree o the polynomial
    sig: array-like
        The uncertainties in the y-axis values
    x: str
        The x-axis parameter name
    y: str
        The y-axis parameter name
    title: str
        The title of the plot
    """
    p, residuals, rank, singular_values, rcond = np.polyfit(np.asarray(n), np.asarray(m), degree,
                                                 w=1./np.array([i if i else 1 for i in sig]) if sig != '' else None, full=True)
    f = np.poly1d(p)
    w = np.linspace(min(n), max(n), 50)
    rms = np.sqrt(sum((m - f(n)) ** 2) / len(n))
    data = [[y, (min(n), max(n)), rms] + list(reversed(p))]

    # Plotting
    if plot:
        ax.plot(w, f(w), color=c, ls=ls, lw=lw, label='${}$'.format(poly_print(p, x=x, y=y)) if legend else '', zorder=10)
        if plot_rms:
            ax.fill_between(w, f(w) - rms, f(w) + rms, color=c, alpha=0.1, zorder=-1)
    
    # Make a dictionary
    D = {'yparam': y, 'xparam': x, 'rms': round(rms, 3), 'min': round(min(n), 1), 'max': round(max(n), 1)}
    D.update({'c{}'.format(str(o)): v for o, v in enumerate(list(reversed(p)))})

    if verbose:
        print_data = np.asarray([[y, r'{:.1f} < {} < {:.1f}'.format(min(n), x, max(n)), '{:.3f}'.format(rms)] 
                   + ['{:.3e}'.format(v) for v in list(reversed(p))]])

        at.Table(print_data, names=['P(x)', 'x', 'rms'] + [r'$c_{}$'.format(str(i)) for i in range(len(p))]).pprint()
        print('\n')

    return D if dictionary else data

def polynomial_eval(values, coeffs, plot=False, color='g', ls='-', lw=2):
    '''
    Evaluates *values* given the list of ascending polynomial *coeffs*.

    Parameters
    ----------
    values: int, list, tuple, array
      The value or values to to evaluated
    coeffs: list, tuple, array
      The sequence of ascending polynomial coefficients beginning with zeroth order
    plot: bool (optional)
      Plot the results in the given color
    color: str
      The color of the line or fill color of the point to plot
    ls: str
      The linestyle of the line to draw
    lw: int
      The linewidth of the line to draw

    Returns
    -------
    out: float, list
      The evaluated results

    '''

    def poly_eval(val):
        return sum([c * (val ** o) for o, c in enumerate(coeffs)])

    if isinstance(coeffs, dict):
        coeffs = [coeffs[j] for j in sorted([i for i in coeffs.keys() if i.startswith('c')])]

    if isinstance(values, (int, float)):
        out = poly_eval(values)
        if plot:
            plt.errorbar([values], [out], marker='*', markersize=18, color=color, markeredgecolor='k',
                         markeredgewidth=2, zorder=10)

    elif isinstance(values, (tuple, list, np.ndarray)):
        out = [poly_eval(v) for v in values]
        if plot:
            plt.plot(values, out, color=color, lw=lw, ls=ls)

    else:
        out = None
        print("Input values must be an integer, float, or sequence of integers or floats!")

    return out

def specType(SpT, types=[i for i in 'OBAFGKMLTY']):
    """
    Converts between float and letter/number spectral types (e.g. 14.5 => 'B4.5' and 'A3' => 23).
    
    Parameters
    ----------
    SpT: float, str
        Float spectral type or letter/number spectral type between O0.0 and Y9.9
    types: list
        The MK spectral type letters to include, e.g. ['M','L','T','Y']
      
    Returns
    -------
    list, str
        The converted spectral type string or (spectral type, luminosity class) numbers
    """
    try:
        # String input
        if isinstance(SpT, str) and SpT[0] in types and SpT!='':
            MK, LC = SpT[0], 'V'
            suf = SpT[1:].replace('n','').replace('e','').replace('w','')\
                         .replace('m','').replace('a','').replace('Fe','')\
                         .replace('-1','').replace(':','').replace('?','')\
                         .replace('-V','').replace('p','')

            for cl in ['III','V','IV']:
                try:
                    idx = suf.find(cl)
                    val = float(suf[:idx].split('/')[0])
                    LC = suf[idx:].split('/')[0].split(',')[0]
                    break
                except:
                    try:
                        val = float(suf)
                    except:
                        continue

            return [types.index(MK)*10+val-(4 if MK in ['M','L','T','Y'] else 0), LC]

        # Numerical input
        elif isinstance(SpT, float) or isinstance(SpT, int) and 0.0 <= SpT < len(types)*10:
            letter = ''.join(types)[int(SpT // 10)]
            number = int(SpT % 10) if SpT%10==int(SpT%10) else SpT%10
            return '{}{}'.format(letter, number)

        # Bogus input
        else:
            print('Spectral type',SpT,'must be a float between 0 and',len(types)*10,'or a string of class',types)
            return [np.nan, '']
        
    except:
        return [np.nan, '']

def proximity_plot(x, y, target, radius, xlabel='x', ylabel='y'):
    """
    Create a scatter plot of sources relative to a target source
    
    Parameters
    ----------
    x: array-like
        The x values
    y: array-like
        The y values
    target: array-like
        The (x,y) values of the target
    radius: float
        The radius from the target to draw
    xlabel: str
        The x axis label
    ylabel: str
        The y axis label
    """
    # Get the target coordinates
    X, Y = target
    
    # Draw the figure
    plt.figure()
    plt.scatter([X], [Y], color='r', marker='*', s=200)
    plt.scatter(x, y, color='b')
    
    # Make the circle
    circle1=plt.Circle(target, radius, color='r', alpha=0.1)
    plt.gcf().gca().add_artist(circle1)
    
    # Fix the labels and limits
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.ylim(min(plt.ylim()[0],Y-radius), max(plt.ylim()[1],Y+radius))
    plt.xlim(min(plt.xlim()[0],X-radius), max(plt.xlim()[1],X+radius))

def pi2pc(parallax, parallax_unc=0):
    """
    Calculate the distance in parsecs given the parallax in mas
    
    Parameters
    ----------
    parallax: float, array-like
        The parallax(s) of the target(s) in [mas]
    parallax_unc: float, array-like
        The uncertainties (same shape and units as the parallax)
    
    Returns
    -------
    astropy.units.quantity.Quantity
        The distance(s) in parsecs
    """
    # Add units and calculate distance
    pi = parallax * q.arcsec / 1000.
    d = (1 * q.pc * q.arcsec) / pi

    # Uncertainties if provided
    if parallax_unc:
        sig_pi = parallax_unc * q.arcsec / 1000.
        sig_d = sig_pi * q.pc * q.arcsec / pi ** 2
    
        return d.round(2), sig_d.round(2)
    
    else:
        return d.round(2)
    
def theta_D(distance, radius='', spt=''):
    """
    Calculate the angular diameter of a source given its distance 
    and radius or spectral type
    
    Parameters
    ----------
    distance: astropy.units.quantity.Quantity
        The distance to the target
    radius: astropy.units.quantity.Quantity (optional)
        The radius of the target
    spt: str
        The spectral type of the target (optional)
        
    Returns
    -------
    float
        The angular diameter of the target
    """
    try:
        if radius=='':
            radius = spt2radius(spt)[0]

        # Decompose the R/D value to unitless quantity
        RD = (radius/distance).decompose()

        return (2*np.arctan(RD)).to(q.marcsec).round(4)
    except:
        return np.nan

def spt2radius(SpT, plot=False):
    """
    Get radius estimate from spectral type
    
    Parameters
    ----------
    SpT: str, array-like
        The spectral type(s) to look up
    
    Returns
    -------
    astropy.units.quantity.Quantity
        The radius given the spectral type(s)
    """
    # Ingest list of sources with spectral type and radius measurements
    # Read in spt/radii tables
    tables = []
    for t in glob.glob(pkg_resources.resource_filename('all_sky_target_tool', 'data/catalogs/radii/*.txt')):
        tab = ii.read(t, header_start=0, delimiter='|')[['spt','radius']]
        tables.append(tab)
    R = at.vstack(tables)

    # Convert spt to integer
    spt = [i.split(',')[0].split('-')[0] for i in R['spt']]
    spt, lc = np.asarray([specType(sp) for sp in spt]).T
    spt = list(map(float,spt))
    R['SpT'], R['lum_class'] = spt, lc

    # Plot it
    dicts = {}
    for cl,c,m,o in zip(['III','IV','V'],['r','m','b'],['s','d','o'],[2,3,5]):
        x, y = np.asarray(R[R['lum_class']==cl]['SpT']), np.asarray(R[R['lum_class']==cl]['radius'])
        x, y = np.asarray(list(map(float,x))), np.asarray(list(map(float,y)))

        # Color cuts
        if cl=='III':
            x, y = x[x>40], y[x>40]
            x, y = x[y>3], y[y>3]
            
            x, y = x[~np.logical_and(x<50,y>50)], y[~np.logical_and(x<50,y>50)]
        elif cl=='V':
            x, y = x[y<13], y[y<13]            
            
        if plot:
            ax = plt.gca()
            ax.scatter(x, y, label=cl, marker=m, color=c)
            ax.set_xlabel('SpT')
            ax.set_ylabel('R_sun')

        # Fit polynomial to data points
        dicts[cl] = polynomial_fit(x, y, sig='', degree=o, x='SpT', y='R_sun', c=c, ax=ax if plot else False, legend=False, plot=plot)

    if plot:
        # Change numbers to spectral types
        ax.set_xlim(-1,66)
        ax.set_xticks([0,10,20,30,40,50,55,65])
        ax.set_xticklabels(['O','B','A','F','G','K','M','L'])
        plt.legend(loc=2)

    # Evaluate the polynomial
    radius, unc = [], []
    for s in SpT:
        try:
            sp, lc = specType(s)
            if s and sp>=dicts[lc]['min'] and sp<=dicts[lc]['max']:
                r = polynomial_eval(sp, dicts[lc], plot=plot)
                rms = dicts[lc]['rms']
            else:
                r, rms = np.nan, np.nan
        except KeyError:
            r, rms = np.nan, np.nan

        radius.append(r)
        unc.append(rms)

    radius = np.asarray(radius)
    unc = np.asarray(unc)
    
    if isinstance(SpT, str):
        radius = radius[0]
        unc = unc[0]
    
    return [radius*q.R_sun, unc*q.R_sun]