#!/usr/bin/env python
from __future__ import print_function, division
import sys
import os
import datetime
import multiprocessing
import subprocess
import shutil
import numpy as np
import requests
from astropy.io import ascii
from astropy.table import Table
from scipy.interpolate import interp1d
from pprint import pformat
from os.path import split, join, exists

# DONE: warn when neighbor criteria push below completeness limit in K for 2MASS
# N/A : flag stars at poles in RA/Dec as potentially having undetected neighbors
#        -> per casual tests and Jay's email, handles singularities at poles
# TODO: generate sky plots of results

from list_specs import target_lists, jay_lists

# Skip anything that touches the filesystem (for debugging)
PRETEND = False
# Query the ecliptic poles only (faster, saves in target_lists_cvz_only)
CVZ_ONLY = True
# Number of workers (32 for telserv1)
N_PROCESSES = 32
# Number of star entries per chunk of list (chunks are inputs to cone
# searches)
# 1000 star chunks leads to 300 - 1500 MB neighbor lists
CHUNK_SIZE = 1000
# Switch off multiprocessing for better tracebacks in debugging
MULTIPROCESS_CHUNKS = True
# Keep intermediate files for debugging
KEEP_INTERMEDIATES = False
# This is maybe the wrong way to handle this cutoff, but here's the
# reasoning. 2MASS is complete to K < 14.3 in "unconfused regions" of
# the sky. Elsewhere in the 2MASS PSC manual, it says the limits can
# be affected by up to 1 mag in regions near the galactic plane,
# so out of an abundance of caution I picked 13.3 as the cutoff where
# 2MASS PSC is supplemented with additional cone searches in GSC-II
#
# See http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec2_2.html
# for details
TWOMASS_COMPLETENESS_K = 13.3
# For debugging, feed a list of stars with supposedly no neighbors
# back into an RDLIST+ search and ensure no neighbors are found
# (doesn't work when the spec says no *brighter* neighbors, since
# that needs more than just query_2mass to apply)
DOUBLE_CHECK_NEIGHBORS = False

RADIUS_DEGREES_FORMAT = '{:1.5f}'
SERVICE_URL_TEMPLATE = "http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?RA={ra:03.9}&DEC={dec:02.8}&SR={rdeg}&FORMAT=CSV&CAT=GSC23"

# These globals are set in the if __name__ == "__main__" block
_pool = None
_manager = None

def _log(*args):
    """
    Helper to log things with the process name (MainProcess, PoolWorker-N, etc.)
    """
    name = multiprocessing.current_process().name
    ts = datetime.datetime.now().isoformat()
    print(
        '[{}:{}]'.format(name, ts) + 
        ' ' + 
        ' '.join(map(str, args))
    )

def check_separation(point_a, point_b, radius_degrees):
    """
    ((long, lat) in deg, (long, lat) in deg, radius in deg) -> True or False

    If the central angle between the two points, as computed by the haversine
    formula, is less than or equal to the radius in degrees, True is returned.

    https://en.wikipedia.org/wiki/Haversine_formula
    """
    import math

    long1, lat1 = point_a
    long2, lat2 = point_b

    lambda1, phi1 = long1 * math.pi / 180.0, lat1 * math.pi / 180.0
    lambda2, phi2 = long2 * math.pi / 180.0, lat2 * math.pi / 180.0

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

def run_command(command_args, input_from=None, output_to=None):
    """
    Take an argument list starting with a path to an executable
    and run, redirecting output to files specified by a path name
    (`input_from` or `output_to`)
    """
    command_string = ' '.join(command_args) + ' '
    if input_from is not None:
        command_string += "< {} ".format(input_from)
    if output_to is not None:
        command_string += "> {}".format(output_to)
    _log(command_string)

    if not PRETEND:
        input_file = open(input_from) if input_from else None
        output_file = open(output_to, 'w') if output_to else None
        try:
            return_code = subprocess.check_call(
                command_args,
                stdin=input_file,
                stdout=output_file,
                stderr=subprocess.STDOUT
            )
        except subprocess.CalledProcessError:
            if input_from is not None:
                input_file.close()
            if output_to is not None:
                output_file.flush()
                output_file.close()
                if not KEEP_INTERMEDIATES:
                    os.remove(output_to)
            raise
        else:
            if input_from is not None:
                input_file.close()
            if output_to is not None:
                output_file.flush()
                output_file.close()

def load_base_stars(filename):
    """A slow / memory-inefficient algorithm to build a map from
    U00000XYZ indices back to the original lines in the star list

    Returns a set `indices` and a dict `lookup` (which is maybe
    inefficient because dict keys are a set, but I think I use the
    set intersection operators somewhere)
    """
    indices = set()
    lookup = {}
    with open(filename, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            idx = line.split()[-1]  # Keep the identifier Jay's code spits out
            lookup[idx] = line
            indices.add(idx)
    return indices, lookup

def chunk_list(base_list_path, chunk_size=CHUNK_SIZE):
    """Takes a path to a list, makes new files with the ".chunk_00000"
    suffix for every `CHUNK_SIZE` lines from the input
    """
    chunk_template = '.chunk_{:05}'
    if PRETEND and not exists(base_list_path):
        return [base_list_path + chunk_template.format(i+1) for i in range(10)]
    line_count = 0
    chunk_count = 1
    if not PRETEND:
        chunk = open(base_list_path + chunk_template.format(1), 'w')
    else:
        chunk = None
    chunk_paths = [base_list_path + chunk_template.format(1)]
    # _log('first chunk = {}'.format(chunk_paths[0]))
    with open(base_list_path, 'r') as base_list:
        for idx, line in enumerate(base_list):
            if not PRETEND:
                chunk.write(line)
            
            if (idx + 1) % chunk_size == 0:
                # write/close chunk, open new file for next
                if not PRETEND:
                    chunk.close()
                chunk_count += 1
                if not PRETEND:
                    chunk = open(base_list_path + chunk_template.format(chunk_count), 'w')
                chunk_paths.append(base_list_path + chunk_template.format(chunk_count))
                # _log('next chunk = {}'.format(chunk_paths[-1]))
    if not PRETEND:
        chunk.close()
    return chunk_paths

def compute_base_list(k_min, k_max):
    """Simply query for stars in some magnitude range, optionally
    limiting to the CVZ. The resulting file path is returned, and will
    be of the form `cache/base_<k_min>_k_<k_max>`"""
    # generate filename
    k_min_string = '{:1.1f}'.format(k_min)
    k_max_string = '{:1.1f}'.format(k_max)
    base_list_name = 'base_{}_k_{}'.format(k_min_string, k_max_string)
    if CVZ_ONLY:
        base_list_name += "_cvz"
    base_list_path = join('cache', base_list_name)
    # build args list
    args = [
        './query_2mass',
        'MK_MIN={}'.format(k_min_string),
        'MK_MAX={}'.format(k_max_string)
    ]
    # test existence
    if exists(base_list_path):
        _log("{} exists".format(base_list_path))
        return base_list_path, count_non_comment_lines(base_list_path)
    if not CVZ_ONLY:
        # call subprocess
        # write to scratch file
        run_command(args, output_to=base_list_path)
    else:
        north_args = args + ['EB_MIN={}'.format(85), 'EB_MAX={}'.format(90)]
        north_path = base_list_path + '.north'
        run_command(north_args, output_to=north_path)

        south_args = args + ['EB_MAX={}'.format(-85), 'EB_MIN={}'.format(-90)]
        south_path = base_list_path + '.south'
        run_command(south_args, output_to=south_path)

        run_command(['cat', north_path, south_path], output_to=base_list_path)

    return base_list_path, count_non_comment_lines(base_list_path)

def prune_stars_with_neighbors(base_list_path, neighbor_list_path, output_list_path, only_reject_brighter_neighbors):
    """Takes a base list (or chunk of a base list) and a corresponding
    neighbor list. Produces an output list at `output_list_path` with
    the stars from `base_list_path` that have *no* neighbors
    in the neighbor list"""
    if exists(output_list_path):
        _log("{} exists".format(output_list_path))
        return output_list_path, count_non_comment_lines(output_list_path)

    if not PRETEND:
        base_indices, base_lookup = load_base_stars(base_list_path)
        with open(neighbor_list_path, 'r') as f:
            comment_lines = 0
            for idx, line in enumerate(f):
                if line[0] == '#':
                    comment_lines += 1
                    continue
                if len(line.strip()) == 0:
                    continue
                parts = line.split()
                if parts[0] == 'END':
                    continue
                base_idx = parts[-1]
                if only_reject_brighter_neighbors:
                    try:
                        neighbor_mag = float(parts[4])
                        base_mag = float(parts[13])
                        if neighbor_mag < base_mag:  # smaller is brighter
                            base_indices.discard(base_idx)
                    except:
                        print(parts)
                        raise
                else:
                    base_indices.discard(base_idx)
        _log("pruning base list {} with neighbor list {}".format(base_list_path, neighbor_list_path))
        _log("base stars: {} entries".format(len(base_lookup.keys())))
        _log("neighbor stars: {} entries".format(
           idx + 1 - comment_lines
        ))
        _log("keeping {} entries from base list".format(len(base_indices)))
        with open(output_list_path, 'w') as f:
            for base_idx in base_indices:
                f.write(base_lookup[base_idx])
    if only_reject_brighter_neighbors:
        _log("create {} with brighter-neighbor-having stars removed".format(output_list_path))
    else:
        _log("create {} with neighbor-having stars removed".format(output_list_path))
    return output_list_path, len(base_indices)

def gsc_prune_stars_with_neighbors(base_list_path, output_list_path, delta_k, radius_arcmin, only_reject_brighter_neighbors):
    radius_degrees_string = RADIUS_DEGREES_FORMAT.format(radius_arcmin / 60.)  # radius to degrees
    base_list = ascii.read(base_list_path, names=['RA', 'Dec', 'J', 'H', 'K', 'qual', 'idx'])
    # output_list = Table(dtype=base_list.dtype)
    # The above doesn't work to make an empty table (it says it has no columns)
    # so be annoyingly thorough in specifying it...
    output_list = Table(names=['RA', 'Dec', 'J', 'H', 'K', 'qual', 'idx'], dtype=['<f8', '<f8', '<f8', '<f8', '<f8', 'S3', 'S9'])
    count_zero_neighbors, count_2mass_artifacts, count_missing_gscmag, count_cant_interpolate, count_approx_k_too_bright = 0, 0, 0, 0, 0
    bj_rf = np.array([0.0, 0.5, 1.0, 2.0, 2.5, 3.0])  # B_J is JpgMag, R_F is FpgMag, bj_rf = JpgMag - FpgMag
    c_bk = np.array([1.5, 1.5, 2.3, 4.8, 7.0, 7.0])
    c_bk_interpolator = interp1d(bj_rf, c_bk)

    for row in base_list:
        query_url = SERVICE_URL_TEMPLATE.format(
            ra=row['RA'],
            dec=row['Dec'],
            rdeg=radius_degrees_string
        )
        resp = requests.get(query_url)
        # make an astropy table out of the votable
        table = ascii.read(resp.text)
        # neighbors will be those stars without KMag (because if they had KMags, they'd be in 2MASS)
        neighbors = table[table['KMag'] == 99.99]
        # find any non-99 kmag values. there should not be any, if the query_2mass tool worked right
        if len(table[(table['KMag'] != 99.99) & (table['class'] == 0)]) != 1:
            # found a 2MASS neighbor that wasn't found in the query_2mass step
            # it's safest to just discard this star; we don't need to transform
            # GSC magnitudes to KMag or anything
            count_2mass_artifacts += 1
            continue
        if len(neighbors) == 0:
            # if no neighbors at *all*, go right to appending to output list
            output_list.add_row(row)
            count_zero_neighbors += 1
            continue
        if np.any(neighbors['FpgMag'] == 99.99) or np.any(neighbors['JpgMag'] == 99.99):
            # can't compute transformed magnitudes for one or more neighbors
            # so we must discard this star
            count_missing_gscmag += 1
            continue
        else:
            # interpolate a K mag using color relation in Anderson 2009
            try:
                c_bk_approx = c_bk_interpolator(neighbors['JpgMag'] - neighbors['FpgMag']) # (B_J - R_F)
            except ValueError:
                _log("Out of range (JpgMag - FpgMag) value processing neighbors of {}".format(base_list_path))
                _log(neighbors)
                count_cant_interpolate += 1
                continue
            # c_bk = b_j - k_2mass -> k_2mass (approx) = b_j - c_bk
            k_approx = neighbors['JpgMag'] - c_bk_approx
            # test length of set of matching rows with abs(computed_neighbor_k_mag - base_k_mag) < delta_k
            neighbor_delta_ks = np.abs(k_approx - row['K'])
            if np.any(neighbor_delta_ks >= delta_k):
                count_approx_k_too_bright += 1
                continue
            else:
                output_list.add_row(row)
    output_list.write(output_list_path, format='ascii.no_header')
    _log("Kept", len(output_list), "of", len(base_list))
    _log("Entries kept because zero neighbors were found:", count_zero_neighbors)
    _log("Entries discrarded because of 2MASS artifacts: ", count_2mass_artifacts)
    _log("Entries discarded because Fpg or Jpg magnitudes weren't available for neighbors:", count_missing_gscmag)
    _log("Entries discarded because can't approximate K mag (input out of range):", count_cant_interpolate)
    _log("Entries discarded because approximate K mag (from Jpg - Fpg color) too bright:", count_approx_k_too_bright)
    return output_list_path, len(output_list)

def find_neighbors(base_list_path, delta_k, radius_arcmin):
    """Query 2MASS for all the neighbors within radius_arcmin arcminutes
    and (if not None) delta_k magnitudes"""
    _, base_list_name = split(base_list_path)

    radius_degrees_string = RADIUS_DEGREES_FORMAT.format(radius_arcmin / 60.)  # radius to degrees
    args = [
        './query_2mass',
        'RDLIST+',
        'DR_MIN=0.001',
        'DR_MAX={}'.format(radius_degrees_string)
    ]
    output_list_path = base_list_path + '.neighbors_r_{}'.format(radius_degrees_string)

    if delta_k is not None:
        delta_k_string = '{:1.1f}'.format(delta_k)
        output_list_path += '_dk_' + delta_k_string
        args.append('DK_MAX={}'.format(delta_k_string))

    run_command(args, input_from=base_list_path, output_to=output_list_path)

    return output_list_path

def _process_chunk_for_neighbors(chunk_path, delta_k, radius_arcmin, output_list_path, output_list_lock, only_reject_brighter_neighbors, must_check_gsc):
    """Used in `find_stars_without_neighbors` to neighbor-search and
    then write out those stars with no neighbors to `output_list_path`

    Uses a multiprocessing Lock (`output_list_lock`) to ensure
    the multiple jobs writing to the output list are well behaved"""
    _log("neighbor searching {}".format(chunk_path))
    _log("args:", (chunk_path, delta_k, radius_arcmin, output_list_path, output_list_lock, only_reject_brighter_neighbors, must_check_gsc))
    starting_total = count_non_comment_lines(chunk_path)
    # run cone search
    chunk_neighbors_path = find_neighbors(chunk_path, delta_k, radius_arcmin)
    # filter base list chunk
    pruned_chunk_path, n_kept_stars = prune_stars_with_neighbors(chunk_path, chunk_neighbors_path, chunk_path + '.pruned', only_reject_brighter_neighbors)
    _log("Chunk pruned with dk < {} r < {}: {} of {} remaining".format(delta_k, radius_arcmin, n_kept_stars, starting_total))
    # pedantically double-check this output
    if DOUBLE_CHECK_NEIGHBORS and delta_k is not None:
        radius_degrees_string = RADIUS_DEGREES_FORMAT.format(radius_arcmin / 60.)  # radius to degrees
        commandstring = "./query_2mass RDLIST+ DR_MIN=0.001 DR_MAX={} DK_MAX={} < {}".format(radius_degrees_string, delta_k, pruned_chunk_path)
        _log(commandstring)
        output = subprocess.check_output(commandstring, shell=True)
        if output != '#NARGs:  4\n# OPEN BASE.DATA for map...\n# OPEN BASE.DATA for work...\n  \n END OF INPUT FILE...\n  \n':
            raise ValueError("Failed double-check of {} with dr {} and dk {}. Output was\n\n".format(pruned_chunk_path, radius_degrees_string, delta_k, output))
        else:
            _log("Checked output from the above command, consistent with all 2MASS-neighbor-having-stars removed")
    # check GSC if necessary
    if must_check_gsc and n_kept_stars > 0:
        pruned_chunk_path, n_kept_stars = gsc_prune_stars_with_neighbors(pruned_chunk_path, pruned_chunk_path + '.gsc_pruned', delta_k, radius_arcmin, only_reject_brighter_neighbors)
        _log("Chunk pruned using GSC with dk < {} r < {}: {} of {} remaining".format(delta_k, radius_arcmin, n_kept_stars, starting_total))

    # append kept stars to final list
    _log("Append remaining stars from {} to {}".format(pruned_chunk_path, output_list_path))
    if PRETEND:
        return  # bail before we touch any files
    if n_kept_stars > 0:
        # don't bother waiting on the lock if there's nothing to write
        if output_list_lock is not None:
            with output_list_lock:
                _log("Acquired lock on {}".format(output_list_path))
                with open(output_list_path, 'a') as output_list:
                    for line in open(pruned_chunk_path):
                        output_list.write(line)
        else:
            with open(output_list_path, 'a') as output_list:
                for line in open(pruned_chunk_path):
                    output_list.write(line)
    # remove chunk and neighbors for chunk
    _log("remove {}".format(pruned_chunk_path))
    _log("remove {}".format(chunk_neighbors_path))
    _log("remove {}".format(chunk_path))
    if not KEEP_INTERMEDIATES:
        os.remove(pruned_chunk_path)
        os.remove(chunk_neighbors_path)
        os.remove(chunk_path)
    return n_kept_stars

def find_stars_without_neighbors(base_list_path, delta_k, radius_arcmin, only_reject_brighter_neighbors, must_check_gsc):
    """Given a base list and constraints on `delta_k`, `radius_arcmin`,
    and whether only brighter neighbors disqualify a target, execute
    a parallelized search for neighbors in chunks from the base list.

    Uses a multiprocessing Pool (`_pool`) and apply_async to farm out
    units of work."""
    # generate filename incorporating base list name
    radius_degrees_string = RADIUS_DEGREES_FORMAT.format(radius_arcmin / 60.)  # radius to degrees
    _, base_list_name = split(base_list_path)
    if delta_k is not None:
        delta_k_string = '{:1.1f}'.format(delta_k)
        neighbor_list_name = 'neighbors_r_{}_dk_{}_for_{}'.format(radius_degrees_string, delta_k_string, base_list_name)
    else:
        delta_k_string = ''
        neighbor_list_name = 'neighbors_r_{}_for_{}'.format(radius_degrees_string, base_list_name)
    neighbor_list_path = join('cache', neighbor_list_name)
    # _log("base name for neighbor search chunks {}".format(neighbor_list_name))
    
    _, base_list_name = split(base_list_path)
    neighbor_list_name, _ = neighbor_list_name.split('_for_')
    output_list_name = '{}_without_{}'.format(base_list_name, neighbor_list_name)
    output_list_path = join('cache', output_list_name)
    if exists(output_list_path):
        return output_list_path, count_non_comment_lines(output_list_path)

    if not PRETEND:
        output_list = open(output_list_path, 'w')
        header = """# base list: {}
# rejecting stars with neighbors
# ... of delta K mag <= {}
# ... within radius <= {} degrees
# ... and only if they're brighter? {}
"""
        output_list.write(header.format(base_list_path, delta_k_string, radius_degrees_string, only_reject_brighter_neighbors))
        output_list.close()
    else:
        output_list = None

    # use a lock to coordinate access to output_list_path
    output_list_lock = _manager.Lock()

    # chunk the base list
    chunk_paths = chunk_list(base_list_path)
    # _log("chunk_paths = {}".format(chunk_paths))
    # neighbor search all the base chunks
    if MULTIPROCESS_CHUNKS:
        results = []
        for chunk_path in chunk_paths:
            _log("Spawning for {}".format(chunk_path))
            res = _pool.apply_async(
                _process_chunk_for_neighbors,
                (
                    chunk_path,
                    delta_k,
                    radius_arcmin,
                    output_list_path,
                    output_list_lock,
                    only_reject_brighter_neighbors,
                    must_check_gsc
                )
            )
            results.append(res)

        n_stars_kept = 0
        for r in results:
            n_stars_kept += r.get()
    else:
        n_stars_kept = 0
        for chunk_path in chunk_paths:
            n_stars_kept += _process_chunk_for_neighbors(
                chunk_path,
                delta_k,
                radius_arcmin,
                output_list_path,
                None,
                only_reject_brighter_neighbors,
                must_check_gsc
            )
    if not PRETEND:
        output_list.close()
    return output_list_path, n_stars_kept

def count_non_comment_lines(path):
    """Counts lines in file that have data (non-empty, non-#-prefixed-comment)"""
    non_comment_lines = 0
    if not PRETEND:
        with open(path, 'r') as f:
            for line in f:
                if len(line) > 0 and line[0] == '#':
                    continue
                elif len(line.strip()) == 0:
                    continue
                else:
                    non_comment_lines += 1
    return non_comment_lines

def remove_non_AAA_sources(path):
    """Remove all sources with PSC qual flags other than AAA
    (in other words, keep only those with good photometry in J/H/K)"""
    outpath = path + "_good"
    n_pruned = 0
    n_base_stars = 0
    if not PRETEND:
        with open(outpath, 'w') as fout:
            with open(path, 'r') as fin:
                for line in fin:
                    if len(line) > 0 and line[0] == '#':
                        continue
                    else:
                        n_base_stars += 1

                    if 'AAA' not in line:
                        n_pruned += 1
                        continue
                    fout.write(line)
        _log("Pruned {} sources with non-AAA quality flags".format(n_pruned))
    n_kept = n_base_stars - n_pruned
    return outpath, n_kept

def compute_list(name, spec):
    """Given a data structure corresponding to target criteria
    (i.e. a dict from `list_specs`), compute the base list and apply
    target criteria that can be handled by this script
    (e.g. neighbor star exclusion, but not enforcing
    guide star availability)"""
    report = [name,]

    def _report(message):
        report.append(message)

    _log("Computing {} from {}".format(name, pformat(spec)))
    k_min, k_max = spec['k_mag']
    base_list_path, n_base_sources = compute_base_list(k_min, k_max)
    msg = "Base list {} < K {} has {} sources".format(k_min, k_max, n_base_sources)
    _log(msg)
    _report(msg)
    base_list_path, n_base_minus_extended = remove_non_AAA_sources(base_list_path)
    msg = "After removing non-AAA sources: {}".format(n_base_minus_extended)
    _log(msg)
    _report(msg)
    neighbor_criteria = spec.get('neighbors', [])

    # neighbor criteria are applied iteratively, so base stars pruned by one set of criteria
    # won't be used as the base for a subsequent cone-search
    # the intermediate_list_path holds the current base list
    intermediate_list_path = base_list_path
    # If one of the criteria puts us over the 2MASS completeness limit
    # we need to go check GSC-II later:
    must_check_gsc = False
    if len(neighbor_criteria) > 0:
        for nc in neighbor_criteria:
            if k_max + nc['delta_k'] > TWOMASS_COMPLETENESS_K:
                warning = "Warning: For {} < K < {} and deltaK < {}, max neighbor Kmag we care about is {}, but 2MASS is only complete to K={}".format(
                    k_min, k_max, nc['delta_k'], k_max + nc['delta_k'], TWOMASS_COMPLETENESS_K
                )
                _log(warning)
                _report(warning)
                must_check_gsc = True
            _log("for {} prune neighbors deltaK < {} mag; r < {} arcmin".format(intermediate_list_path, nc['delta_k'], nc['r_arcmin']))
            intermediate_list_path, n_kept = find_stars_without_neighbors(intermediate_list_path, nc['delta_k'], nc['r_arcmin'], only_reject_brighter_neighbors=False, must_check_gsc=must_check_gsc)
            _log("Got {}".format(intermediate_list_path))
            msg = "After excluding stars with neighbors in deltaK < {} mag; r < {} arcmin: {} sources".format(nc['delta_k'], nc['r_arcmin'], n_kept)
            _log(msg)
            _report(msg)

    brighter_neighbors_r_arcmin = spec.get('no_brighter_neighbors_r_arcmin')
    if brighter_neighbors_r_arcmin is not None:
        # can't supply DK_MAX here, since we'll miss some bright neighbors unless we set DK_MAX high enough that it's useless
        _log("for {} prune brighter neighbors in r {}".format(intermediate_list_path, brighter_neighbors_r_arcmin))
        pruned_list_path, n_kept = find_stars_without_neighbors(intermediate_list_path, None, brighter_neighbors_r_arcmin, only_reject_brighter_neighbors=True, must_check_gsc=False)
        _log("Got {}".format(pruned_list_path))
        msg = "After excluding stars with brighter neighbors in r < {}: {} sources".format(brighter_neighbors_r_arcmin, n_kept)
        _log(msg)
        _report(msg)
    else:
        pruned_list_path = intermediate_list_path

    if CVZ_ONLY:
        dest_path = join('target_lists_cvz_only', name)
        _report("CVZ only")
    else:
        dest_path = join('target_lists', name)
    _log("cp {} {}".format(pruned_list_path, dest_path))
    if not PRETEND:
        shutil.copy(pruned_list_path, dest_path)

    with open(dest_path + '.report', 'w') as f:
        f.write('\n'.join(report))
        f.write("\n")
        _log("Wrote report to {}".format(dest_path + ".report"))
    return dest_path

if __name__ == "__main__":
    # Set up these shared/global variables
    _pool = multiprocessing.Pool(N_PROCESSES)
    _manager = multiprocessing.Manager()
    # Make sure destination directories exist
    subprocess.call('mkdir -p ./cache ./target_lists ./target_lists_cvz_only', shell=True)

    # Without writing a full dependency solver, this should be enough to ensure
    # that target lists sharing the same base mag criteria don't clobber
    # each other when run with multiprocessing
    # early_commissioning, global_alignment - both 4.5-5.5
    compute_list('early_commissioning', target_lists['early_commissioning'])
    # coarse_phasing, fine_phasing_routine_maintenance - both 8.5-9.5
    compute_list('coarse_phasing', target_lists['coarse_phasing'])


    # Now that intermediate lists are in place for those that conflict,
    # do the rest all at once
    rest_of_the_target_lists = {
        'initial_image_mosaic': target_lists['initial_image_mosaic'],
        'global_alignment': target_lists['global_alignment'],
        'fine_phasing_routine_maintenance': target_lists['fine_phasing_routine_maintenance'],
        # 'mimf_miri': target_lists['mimf_miri'],
    }

    for name, spec in rest_of_the_target_lists.items():
        compute_list(name, spec)
    _pool.close()
    _pool.join()
