# -*- coding: utf-8 -*-
from __future__ import print_function
import glob
import matplotlib
import matplotlib.animation as animation
from matplotlib import pyplot as plt
import numpy as np
import datetime
import os.path
from astropy.io import ascii
import ephem
import pdb

def read_commstars(filename):
    """Read in a text table with columns for RA, Dec, J, H, K mag,
    quality flags (qual), and an index (discarded), and add
    'el' and 'eb' columns for ecliptic longitude and latitude
    respectively."""

    table = ascii.read(filename, names=['RA', 'Dec', 'J', 'H', 'K', 'qual', 'idx'])
    table_ecl = np.empty(len(table), dtype=[
        ('RA', '<f8'),
        ('Dec', '<f8'),
        ('el', '<f8'),
        ('eb', '<f8'),
        ('J', '<f8'),
        ('H', '<f8'),
        ('K', '<f8'),
        ('qual', 'S3'),
    ])
    for idx, row in enumerate(table):
        l, b = equatorial_deg_to_ecliptic_deg(row['RA'], row['Dec'])
        new_row = (row['RA'], row['Dec'], l, b, row['J'], row['H'], row['K'], row['qual'])
        table_ecl[idx] = new_row

    return table_ecl

def read_jaystars(filename):
    try:
        table = ascii.read(filename, names=['RA', 'Dec', 'K',])
    except ascii.InconsistentTableError:
        # if it's not one format of Jay table, it's the other
        table = ascii.read(filename, names=['RA', 'Dec', 'K', 'dx1', 'dy1', 'mm1'])
        
    table_ecl = np.empty(len(table), dtype=[
        ('RA', '<f8'),
        ('Dec', '<f8'),
        ('el', '<f8'),
        ('eb', '<f8'),
        ('K', '<f8'),
    ])
    for idx, row in enumerate(table):
        l, b = equatorial_deg_to_ecliptic_deg(row['RA'], row['Dec'])
        new_row = (row['RA'], row['Dec'], l, b, row['K'])
        table_ecl[idx] = new_row

    return table_ecl

def equatorial_deg_to_ecliptic_deg(ra, dec):
    """Convert RA and declination as decimal degrees to
    ecliptic longitude and latitude in decimal degrees"""
    eq = ephem.Equatorial(ra * ephem.degree, dec * ephem.degree)
    ec = ephem.Ecliptic(eq)
    return ec.lon / ephem.degree, ec.lat / ephem.degree

def ecliptic_deg_to_equatorial_deg(lon, lat):
    """Convert ecliptic longitude and latitude as decimal degrees to
    equatorial RA and declination as decimal degrees"""
    ec = ephem.Ecliptic(lon * ephem.degree, lat * ephem.degree)
    eq = ephem.Equatorial(ec)
    return eq.ra / ephem.degree, eq.dec / ephem.degree

def equatorial_plane_rad():
    """Get ecliptic longitude and latitude coordinates in radians
    for the celestial equator (declination = 0 deg)"""
    equatorial_plane_ra = np.linspace(0, 360, 200)
    equatorial_plane_dec = np.zeros_like(equatorial_plane_ra)
    equatorial_plane_l, equatorial_plane_b = np.empty_like(equatorial_plane_ra), np.empty_like(equatorial_plane_dec)
    for idx, (ra, dec) in enumerate(zip(equatorial_plane_ra, equatorial_plane_dec)):
        l, b = equatorial_deg_to_ecliptic_deg(ra, dec)
        equatorial_plane_l[idx] = l
        equatorial_plane_b[idx] = b
    return np.deg2rad(equatorial_plane_l), np.deg2rad(equatorial_plane_b)

def galactic_plane_rad():
    """Get ecliptic longitude and latitude coordinates in radians
    for the galactic plane"""
    galactic_l = np.linspace(0, 2 * np.pi, 200)
    galactic_b = np.zeros_like(galactic_l)
    ecliptic_l, ecliptic_b = np.empty_like(galactic_l), np.empty_like(galactic_b)
    for idx, (g_l, g_b) in enumerate(zip(galactic_l, galactic_b)):
        gal = ephem.Galactic(g_l, g_b)
        ec = ephem.Ecliptic(gal)
        ecliptic_l[idx] = ec.lon
        ecliptic_b[idx] = ec.lat
    return ecliptic_l, ecliptic_b

def dest_from_start_and_bearing_rad(lon, lat, bearing, angular_dist):
    """
    Find a destination's latitude and longitude given a
    starting (longitude, latitude) pair, a bearing in
    [0, 2pi], and an angular distance. (All arguments in radians.)
    
    From the explanation at
    http://www.movable-type.co.uk/scripts/latlong.html#destPoint
    """
    lat2 = np.arcsin(np.sin(lat) * np.cos(angular_dist) + np.cos(lat) * np.sin(angular_dist) * np.cos(bearing))
    lon2 = lon + np.arctan2(np.sin(bearing) * np.sin(angular_dist) * np.cos(lat), np.cos(angular_dist) - np.sin(lat) * np.sin(lat2))
    return lon2, lat2

def small_circle_rad(lon, lat, angular_radius, npoints=100):
    """Compute longitude and latitude coordinates for `npoints`
    points in a circle around longitude `lon` and latitude `lat` in radians"""
    lons = np.ones(npoints) * lon
    lats = np.ones(npoints) * lat
    angular_radii = np.ones(npoints) * angular_radius
    bearings = np.linspace(0, 2 * np.pi, npoints)
    lonpts, latpts = dest_from_start_and_bearing_rad(lons, lats, bearings, angular_radii)
    return lonpts, latpts

def _plot_deg(func, lon, lat, **kwargs):
    """Take longitude degrees [0, 360] to radians [-pi, pi], latitude degrees to radians, and plot"""
    if not np.isscalar(lon):
        lon = lon.copy()
        lon %= 360.0
        lon[lon > 180] -= 360
        lon_rad, lat_rad = nanify(np.deg2rad(lon), np.deg2rad(lat))
    else:
        lon %= 360.0
        lon = lon - 360 if lon > 180 else lon
        lon_rad, lat_rad = np.deg2rad(lon), np.deg2rad(lat)
    func(lon_rad, lat_rad, **kwargs)

def _plot_rad(func, lon, lat, **kwargs):
    """Wrap radians at longitude pi and plot"""
    if not np.isscalar(lon):
        lon = lon.copy()
        lon %= 2 * np.pi
        lon[lon > np.pi] -= 2 * np.pi
        lon, lat = nanify(lon, lat)
    else:
        lon %= 2 * np.pi
        lon = lon - 2 * np.pi if lon > np.pi else lon
    func(lon, lat, **kwargs)

def plot_deg(ax, lon, lat, **kwargs):
    _plot_deg(ax.plot, lon, lat, **kwargs)
def scatter_deg(ax, lon, lat, **kwargs):
    _plot_deg(ax.scatter, lon, lat, **kwargs)

def plot_rad(ax, lon, lat, **kwargs):
    _plot_rad(ax.plot, lon, lat, **kwargs)
def scatter_rad(ax, lon, lat, **kwargs):
    _plot_rad(ax.scatter, lon, lat, **kwargs)

def nanify(lon, lat):
    """Take longitude and latitude in radians, replace
    certain `lon` values where coordinates wrap around
    with `np.nan` to create discontinuities
    (and prevent weird wrapping behavior)"""
    wjump = np.where(((lon[:-1] >  np.pi/2) & (lon[1:] < -np.pi/2)) | 
                     ((lon[:-1] < -np.pi/2) & (lon[1:] >  np.pi/2)))
    lat = lat.copy()
    lat[wjump] = np.nan

    return lon, lat

sun_angle_from, sun_angle_to = 85, 85 + 50  # can point up to 5 deg towards sun, up to 50 deg away from sun

def check_separation(point_a, point_b, radius_degrees):
    """
    ((long, lat) in deg, (long, lat) in deg, radius in deg) -> True or False

    If the central angle between the two points, as computed by the haversine
    formula, is less than or equal to the radius in degrees, True is returned.

    https://en.wikipedia.org/wiki/Haversine_formula
    """
    import math

#     long1, lat1 = point_a
#     long2, lat2 = point_b

#     lambda1, phi1 = long1 * math.pi / 180.0, lat1 * math.pi / 180.0
#     lambda2, phi2 = long2 * math.pi / 180.0, lat2 * math.pi / 180.0

    def haversine(theta):
        return math.sin(theta / 2.0) ** 2

    hav_d_over_r = haversine(phi2 - phi1) + \
        math.cos(phi1) * math.cos(phi2) * haversine(lambda2 - lambda1)

    central_angle_rad = 2 * math.asin(math.sqrt(hav_d_over_r))
    central_angle_deg = central_angle_rad * 180.0 / math.pi

    if central_angle_deg <= radius_degrees:
        return True
    else:
        return False


def angular_separation_rad(lambda1, phi1, lambda2, phi2):
    """
    Compute the angular separation from (lambda1, phi1)
    in radians to (lambda2, phi2) in radians.
    
    lambda1, phi1 : arrays or scalars
    lambda2, phi2 : scalars
    """
    def haversine(theta):
        return np.sin(theta / 2.0) ** 2
    hav_d_over_r = haversine(phi2 - phi1) + \
        np.cos(phi1) * np.cos(phi2) * haversine(lambda2 - lambda1)
    central_angle_rad = 2 * np.arcsin(np.sqrt(hav_d_over_r))
    return central_angle_rad
    
def field_of_regard_filter(catalog, sun):
    lambdas, phis = np.deg2rad(catalog['el']), np.deg2rad(catalog['eb'])
    separations = angular_separation_rad(lambdas, phis, sun.lon, sun.lat)
    sun_angle_from_rad, sun_angle_to_rad = np.deg2rad(sun_angle_from), np.deg2rad(sun_angle_to)
    rows_with_separation_in_range = (separations > sun_angle_from_rad) & (separations < sun_angle_to_rad)

    return rows_with_separation_in_range

def plot_sky(sky_ax, availability_ax, frame_num, catalog, sun_instances, availability_flags):
    """Plot the sky in ecliptic coordinates for a given catalog and date"""
    sky_ax.set_title('Ecliptic Coordinates')
    sky_ax.grid()
    scatter_deg(sky_ax, catalog['el'], catalog['eb'], marker='.', alpha=0.5)
    # catalog_in_field = catalog[field_of_regard_filter(catalog, sun_ec)]
    catalog_in_field = catalog[availability_flags[:, frame_num]]
    scatter_deg(sky_ax, catalog_in_field['el'], catalog_in_field['eb'], marker='.', color='red', alpha=0.5)

    sun_ec = sun_instances[frame_num]
    # plot sun
    plot_deg(sky_ax, sun_ec.lon / ephem.degree, sun_ec.lat / ephem.degree, markersize=15, marker='o', color='yellow')
    plot_deg(sky_ax, sun_ec.lon / ephem.degree + 180, 0, markersize=15, marker='o', color='red')

    # add equatorial plane
    equatorial_plane_l, equatorial_plane_b = equatorial_plane_rad()
    plot_rad(sky_ax, equatorial_plane_l, equatorial_plane_b, label='Equator', lw=2, alpha=0.5)

    # add galactic plane
    gal_l, gal_b = galactic_plane_rad()
    plot_rad(sky_ax, gal_l, gal_b, lw=2, label='Galactic', alpha=0.5)
    
    # plot field of regard
    cl, cb = small_circle_rad(sun_ec.lon, sun_ec.lat, np.deg2rad(sun_angle_from))
    plot_rad(sky_ax, cl, cb, lw=4, label=u"{}º to {}º from sun".format(sun_angle_from, sun_angle_to), color='orange')
    cl, cb = small_circle_rad(sun_ec.lon, sun_ec.lat, np.deg2rad(sun_angle_to))
    plot_rad(sky_ax, cl, cb, lw=4, color='orange')
    cl, cb = small_circle_rad(sun_ec.lon, sun_ec.lat, np.deg2rad((sun_angle_to + sun_angle_from) / 2.0))
    plot_rad(sky_ax, cl, cb, lw=4, ls='--', color='orange')
    # sky_ax.legend()
    
    # plot stars available vs. day of year
    days = np.arange(availability_flags.shape[1])
    star_totals = availability_flags.sum(axis=0)
    availability_ax.plot(days, star_totals)
    availability_ax.axvline(x=frame_num)

def plot_sky2(sky_ax, availability_ax, frame_num, catalog, reduc_catalog, sun_instances, availability_flags):
    """Plot the sky in ecliptic coordinates for a given catalog and date"""
    sky_ax.set_title('Ecliptic Coordinates')
    sky_ax.grid()
    scatter_deg(sky_ax, catalog['el'], catalog['eb'], marker='.', alpha=0.5)
    scatter_deg(sky_ax, reduc_catalog['el'], reduc_catalog['eb'], marker='*', color='gray', edgecolor='gray', alpha=0.5, s=40)
    # catalog_in_field = catalog[field_of_regard_filter(catalog, sun_ec)]
    catalog_in_field = reduc_catalog[availability_flags[:, frame_num]]
    # scatter_deg(sky_ax, catalog_in_field['el'], catalog_in_field['eb'], marker='.', color='red', alpha=0.5)
    scatter_deg(sky_ax, catalog_in_field['el'], catalog_in_field['eb'], marker='*', color='red', edgecolor='red', s=40)

    sun_ec = sun_instances[frame_num]
    # plot sun
    plot_deg(sky_ax, sun_ec.lon / ephem.degree, sun_ec.lat / ephem.degree, markersize=15, marker='o', color='yellow')
    plot_deg(sky_ax, sun_ec.lon / ephem.degree + 180, 0, markersize=15, marker='o', color='red')

    # add equatorial plane
    equatorial_plane_l, equatorial_plane_b = equatorial_plane_rad()
    plot_rad(sky_ax, equatorial_plane_l, equatorial_plane_b, label='Equator', lw=2, alpha=0.5)

    # add galactic plane
    gal_l, gal_b = galactic_plane_rad()
    plot_rad(sky_ax, gal_l, gal_b, lw=2, label='Galactic', alpha=0.5)
    
    # plot field of regard
    cl, cb = small_circle_rad(sun_ec.lon, sun_ec.lat, np.deg2rad(sun_angle_from))
    plot_rad(sky_ax, cl, cb, lw=4, label=u"{}º to {}º from sun".format(sun_angle_from, sun_angle_to), color='orange')
    cl, cb = small_circle_rad(sun_ec.lon, sun_ec.lat, np.deg2rad(sun_angle_to))
    plot_rad(sky_ax, cl, cb, lw=4, color='orange')
    cl, cb = small_circle_rad(sun_ec.lon, sun_ec.lat, np.deg2rad((sun_angle_to + sun_angle_from) / 2.0))
    plot_rad(sky_ax, cl, cb, lw=4, ls='--', color='orange')
    # sky_ax.legend()
    
    # plot stars available vs. day of year
    days = np.arange(availability_flags.shape[1])
    # star_totals = availability_flags.sum(axis=0)
    #availability_ax.plot(days, star_totals)
    #availability_ax.imshow(availability_flags, interpolation='nearest', cmap='gray_r', aspect='auto')
    availability_ax.imshow(availability_flags.astype(float)*0.5, vmax=1, interpolation='nearest', cmap='gray_r', aspect='auto')
    availability_ax.axvline(x=frame_num)
    N_stars = availability_flags.shape[0]
    plt.yticks(np.arange(N_stars), (np.arange(N_stars)+1)[::-1])

def format_date(date_or_datetime):
    return '{}/{:02}/{:02}'.format(date_or_datetime.year, date_or_datetime.month, date_or_datetime.day)

def precompute_availability(catalog, start_date, n_days):
    sun_instances = []
    availability_flags = np.zeros((len(catalog), n_days), dtype=np.bool)
    for i in range(n_days):
        sun = ephem.Sun()
        date_string = format_date(start_date + datetime.timedelta(i))
        sun.compute(date_string)
#        print(date_string, "Sun RA = ", sun.ra, " Dec = ", sun.dec)
        sun_ec = ephem.Ecliptic(sun)
#        print(date_string, "Sun l = ", sun_ec.lon, " b = ", sun_ec.lat)
        bitmask = field_of_regard_filter(catalog, sun_ec)
        availability_flags[:,i] = bitmask
        sun_instances.append(sun_ec)

    return sun_instances, availability_flags

#def analyze_catalog(catalog_path, kind=None):
def analyze_catalog(catalog_path, reduc_catalog_path, kind='jay'):
    _, catalog_name = os.path.split(catalog_path)
    catalog_name, _ = os.path.splitext(catalog_name)
    if kind == 'jay':
        catalog = read_jaystars(catalog_path)
    else:
        catalog = read_commstars(catalog_path)

    if reduc_catalog_path is not None:
        _, reduc_catalog_name = os.path.split(reduc_catalog_path)
        reduc_catalog_name, _ = os.path.splitext(reduc_catalog_name)
        if kind == 'jay':
            reduc_catalog = read_jaystars(reduc_catalog_path)
        else:
            reduc_catalog = read_commstars(reduc_catalog_path)
    else:
        reduc_catalog = None

#    pdb.set_trace()
    n_days = int(1.5 * 365)
    launch_date = datetime.datetime(2018, 10, 1)  # placeholder

    # First light will occur at about 28 days after launch, initiating
    # wavefront sensing and control activities to align the mirror segments.
    # Instrument checkout will start 37 days after launch, well before the
    # final L2 orbit insertion is obtained after 106 days. Hereafter the
    # full commissioning starts and the observatory will be ready for normal
    # science operations approximately 6 months after launch.
    commissioning_begins = datetime.datetime(2018, 10, 1) + datetime.timedelta(days=28)
    suns, availability = precompute_availability(
        catalog,
        commissioning_begins,
        n_days=n_days
    )
    avail_fname = "%s_avail.npy" % catalog_name
    np.save(avail_fname, availability)
    
    fig = plt.figure(figsize=(14, 8))
    sky_ax = plt.subplot2grid(
        (3,2), (0, 0),
        colspan=2, rowspan=2,
        projection='mollweide'
    )
    availability_ax = plt.subplot2grid(
        (3,2), (2, 0),
        colspan=2, rowspan=1
    )

    def update_for_day(frame_num):
        # current_date = commissioning_begins + datetime.timedelta(days=frame_num)
        sky_ax.clear()
        availability_ax.clear()
        plot_sky(sky_ax, availability_ax, frame_num, catalog, suns, availability)
        availability_ax.set_xlabel('Day since {}'.format(commissioning_begins))
        availability_ax.set_xlim(0, n_days)
        availability_ax.set_ylabel('N targets')
        sky_ax.set_title('{} ({} + {})'.format(
            commissioning_begins + datetime.timedelta(frame_num),
            commissioning_begins,
            datetime.timedelta(frame_num),
        ))

    psf_ani = animation.FuncAnimation(fig, update_for_day, n_days,
                                      interval=100, blit=True)
    psf_ani.save(
        '{}.mp4'.format(catalog_name),
        writer='ffmpeg', bitrate=2000
#        writer='ffmpeg', bitrate=200
#        '{}.gif'.format(catalog_name),
#        writer='imagemagick'
    )
    plt.clf()

    if reduc_catalog is not None:
        print(reduc_catalog)
        suns, availability = precompute_availability(
            reduc_catalog,
            commissioning_begins,
            n_days=n_days
        )
        fig2 = plt.figure(figsize=(14, 8))
        sky_ax = plt.subplot2grid(
            (3,2), (0, 0),
            colspan=2, rowspan=2,
            projection='mollweide'
        )
        availability_ax = plt.subplot2grid(
            (3,2), (2, 0),
            colspan=2, rowspan=1
        )
        def update_for_day2(frame_num):
            # current_date = commissioning_begins + datetime.timedelta(days=frame_num)
            sky_ax.clear()
            availability_ax.clear()
            plot_sky2(sky_ax, availability_ax, frame_num, catalog, reduc_catalog, suns, availability)
            availability_ax.set_xlabel('Day since {}'.format(commissioning_begins))
            availability_ax.set_xlim(0, n_days)
            # availability_ax.set_ylabel('N targets')
            sky_ax.set_title('{} ({} + {})'.format(
                commissioning_begins + datetime.timedelta(frame_num),
                commissioning_begins,
                datetime.timedelta(frame_num),
            ))
        reduc_psf_ani = animation.FuncAnimation(fig2, update_for_day2, n_days,
                                                interval=100, blit=True)
        reduc_psf_ani.save(
            '{}.mp4'.format(reduc_catalog_name),
            writer='ffmpeg', bitrate=2000
        )
        plt.clf()

if __name__ == "__main__":
    import sys
#    analyze_catalog(sys.argv[-1],'not')
    if len(sys.argv) == 2:
        analyze_catalog(sys.argv[1], None, 'not')
    elif len(sys.argv) > 2:
        analyze_catalog(sys.argv[1], sys.argv[2], 'not')
