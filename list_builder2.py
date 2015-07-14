#!/usr/bin/env python
from __future__ import print_function, division
import sys
import os
import datetime
import multiprocessing
import subprocess
import shutil
from pprint import pformat
from os.path import split, join, exists

# TODO: warn when neighbor criteria push below completeness limit in K for 2MASS
# TODO: flag stars at poles in RA/Dec as potentially having undetected neighbors
# ----> per casual tests and Jay's email, handles singularities at poles
# TODO: generate sky plots of results

from list_specs import target_lists

PRETEND = False
CHUNK_SIZE = 1000 # 1000 star chunks leads to 300 - 1500 MB neighbor lists
TWOMASS_COMPLETENESS_K = 15.5
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
            input_file.close()
            output_file.flush()
            output_file.close()
            os.remove(output_to)
            raise
        else:
            input_file.close()
            output_file.flush()
            output_file.close()

def load_base_stars(filename):
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
    chunk_template = '.chunk_{:03}'
    if PRETEND and not exists(base_list_path):
        return [base_list_path + chunk_template.format(i+1) for i in range(10)]
    line_count = 0
    chunk_count = 1
    if not PRETEND:
        chunk = open(base_list_path + chunk_template.format(1), 'w')
    else:
        chunk = None
    chunk_paths = [base_list_path + chunk_template.format(1)]
    _log('first chunk = {}'.format(chunk_paths[0]))
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
                _log('next chunk = {}'.format(chunk_paths[-1]))
    if not PRETEND:
        chunk.close()
    return chunk_paths

def compute_base_list(k_min, k_max):
    # generate filename
    k_min_string = '{:1.1f}'.format(k_min)
    k_max_string = '{:1.1f}'.format(k_max)
    base_list_name = 'base_{}_k_{}'.format(k_min_string, k_max_string)
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
        return base_list_path
    # call subprocess
    # write to scratch file
    run_command(args, output_to=base_list_path)

    # return file path
    return base_list_path

def prune_stars_with_neighbors(base_list_path, neighbor_list_path, output_list_path, only_reject_brighter_neighbors):
    if exists(output_list_path):
        _log("{} exists".format(output_list_path))
        return output_list_path

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
        _log("# base stars: {} entries".format(len(base_lookup.keys())))
        _log("# neighbor stars: {} entries".format(
           idx + 1 - comment_lines
        ))
        _log("# keeping {} entries from base list".format(len(base_indices)))
        with open(output_list_path, 'w') as f:
            for base_idx in base_indices:
                f.write(base_lookup[base_idx])
    if only_reject_brighter_neighbors:
        _log("create {} with brighter-neighbor-having stars removed".format(output_list_path))
    else:
        _log("create {} with neighbor-having stars removed".format(output_list_path))
    return output_list_path

def find_neighbors(base_list_path, delta_k, radius_arcmin):
    _, base_list_name = split(base_list_path)

    radius_degrees_string = '{:1.3f}'.format(radius_arcmin / 60.)  # radius to degrees
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

def _process_chunk_for_neighbors(chunk_path, delta_k, radius_arcmin, output_list_path, output_list_lock, only_reject_brighter_neighbors):
    _log("neighbor searching {}".format(chunk_path))
    # run cone search
    chunk_neighbors_path = find_neighbors(chunk_path, delta_k, radius_arcmin)
    # filter base list chunk
    pruned_chunk_path = prune_stars_with_neighbors(chunk_path, chunk_neighbors_path, chunk_neighbors_path + '.pruned', only_reject_brighter_neighbors)

    # append kept stars to final list
    _log("Append remaining stars from {} to {}".format(pruned_chunk_path, output_list_path))
    if PRETEND:
        return  # bail before we touch any files
    with output_list_lock:
        _log("Acquired lock on {}".format(output_list_path))
        with open(output_list_path, 'a') as output_list:
            for line in open(pruned_chunk_path):
                output_list.write(line)
    # remove chunk and neighbors for chunk
    _log("remove {}".format(pruned_chunk_path))
    # os.remove(pruned_chunk_path)
    _log("remove {}".format(chunk_neighbors_path))
    # os.remove(chunk_neighbors_path)
    _log("remove {}".format(chunk_path))
    # os.remove(chunk_path)

def find_stars_without_neighbors(base_list_path, delta_k, radius_arcmin, only_reject_brighter_neighbors):
    # generate filename incorporating base list name
    radius_degrees_string = '{:1.3f}'.format(radius_arcmin / 60.)  # radius to degrees
    _, base_list_name = split(base_list_path)
    if delta_k is not None:
        delta_k_string = '{:1.1f}'.format(delta_k)
        neighbor_list_name = 'neighbors_r_{}_dk_{}_for_{}'.format(radius_degrees_string, delta_k_string, base_list_name)
    else:
        delta_k_string = ''
        neighbor_list_name = 'neighbors_r_{}_for_{}'.format(radius_degrees_string, base_list_name)
    neighbor_list_path = join('cache', neighbor_list_name)
    _log("base name for neighbor search chunks {}".format(neighbor_list_name))
    
    _, base_list_name = split(base_list_path)
    neighbor_list_name, _ = neighbor_list_name.split('_for_')
    output_list_name = '{}_without_{}'.format(base_list_name, neighbor_list_name)
    output_list_path = join('cache', output_list_name)
    if exists(output_list_path):
        return output_list_path

    if not PRETEND:
        output_list = open(output_list_path, 'w')
        header = """
# base list: {}
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
    _log("chunk_paths = {}".format(chunk_paths))
    # neighbor search all the base chunks
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
                only_reject_brighter_neighbors
            )
        )
        results.append(res)

    for r in results:
        r.get()  # n.b. even though the helper returns None, want to see any exceptions propagate up

    # for idx, chunk_path in enumerate(chunk_paths):
    #     _log("neighbor searching {}".format(chunk_path))
    #     # run cone search
    #     chunk_neighbors_path = find_neighbors(chunk_path, delta_k, radius_arcmin)
    #     # filter base list chunk
    #     pruned_chunk_path = prune_stars_with_neighbors(chunk_path, chunk_neighbors_path, chunk_neighbors_path + '.pruned', delta_k)
    #
    #     # append kept stars to final list
    #     _log("Append remaining stars from {} to {}".format(pruned_chunk_path, output_list_path))
    #     if PRETEND:
    #         continue  # bail before we touch any files
    #     for line in open(pruned_chunk_path):
    #         output_list.write(line)
    #     # remove chunk and neighbors for chunk
    #     _log("remove {}".format(pruned_chunk_path))
    #     os.remove(pruned_chunk_path)
    #     _log("remove {}".format(chunk_neighbors_path))
    #     os.remove(chunk_neighbors_path)
    #     _log("remove {}".format(chunk_path))
    #     os.remove(chunk_path)

    if not PRETEND:
        output_list.close()
    return output_list_path

def compute_list(name, spec):
    _log("Computing {} from {}".format(name, pformat(spec)))
    k_min, k_max = spec['k_mag']
    base_list_path = compute_base_list(k_min, k_max)
    neighbor_criteria = spec.get('neighbors', [])

    # neighbor criteria are applied iteratively, so base stars pruned by one set of criteria
    # won't be used as the base for a subsequent cone-search
    # the intermediate_list_path holds the current base list
    intermediate_list_path = base_list_path
    if len(neighbor_criteria) > 0:
        for nc in neighbor_criteria:
            if k_max + nc['delta_k'] > TWOMASS_COMPLETENESS_K:
                raise RuntimeError("For {} < K < {} and deltaK <= {}, max neighbor Kmag we care about is {}, but 2MASS is only complete to K={}".format(
                    k_min, k_max, k_max + nc['delta_k'], 
                ))
            _log("for {} prune neighbors deltaK <= {} mag; r <= {} arcmin".format(intermediate_list_path, nc['delta_k'], nc['r_arcmin']))
            intermediate_list_path = find_stars_without_neighbors(intermediate_list_path, nc['delta_k'], nc['r_arcmin'], only_reject_brighter_neighbors=False)

    brighter_neighbors_r_arcmin = spec.get('no_brighter_neighbors_r_arcmin')
    if brighter_neighbors_r_arcmin is not None:
        # can't supply DK_MAX here, since we'll miss some bright neighbors unless we set DK_MAX high enough that it's useless
        _log("for {} prune brighter neighbors in r {}".format(intermediate_list_path, brighter_neighbors_r_arcmin))
        pruned_list_path = prune_neighbors(intermediate_list_path, None, brighter_neighbors_r_arcmin, only_reject_brighter_neighbors=True)
    else:
        pruned_list_path = intermediate_list_path

    dest_path = join('target_lists', name)
    _log("cp {} {}".format(pruned_list_path, dest_path))
    if not PRETEND:
        shutil.copy(pruned_list_path, join('target_lists', name))

    return dest_path

def _test_log():
    _log("test message")
    return 'foo'

if __name__ == "__main__":
    # Set up these shared/global variables
    _pool = multiprocessing.Pool(8)  # 8 processors in my Mac Pro
    _manager = multiprocessing.Manager()
    # Make sure destination directories exist
    subprocess.call('mkdir -p ./cache ./target_lists', shell=True)
    compute_list('early_commissioning', target_lists['early_commissioning'])
    res = _pool.apply_async(_test_log)
    res.get()
    _pool.close()
    _pool.join()
    