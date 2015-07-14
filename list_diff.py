#!/usr/bin/env python
from __future__ import division, print_function
from itertools import chain
import decimal

def load_list(list_path, ra_col, dec_col, precision):
    list_dict = {}
    with open(list_path) as f:
        for line in f:
            if line[0] == '#':
                continue
            parts = line.split()
            if len(parts) == 0:
                continue
            ra1, ra2 = parts[ra_col - 1].split('.')
            ra2 = ra2[:precision]
            ra = ra1 + '.' + ra2
            dec1, dec2 = parts[dec_col - 1].split('.')
            dec2 = dec2[:precision]
            dec = dec1 + '.' + dec2
            list_dict[(decimal.Decimal(ra), decimal.Decimal(dec))] = line
    return list_dict

def compare_lists(left, right):
    left_entries = set(left.iterkeys())
    right_entries = set(right.iterkeys())
    left_only = left_entries - right_entries
    right_only = right_entries - left_entries

    print("left: {} entries".format(len(left_entries)))
    if len(left_only) > 0:
        print("These {} stars only appeared in the left list:".format(len(left_only)))
        for key in left_only:
            print("\t" + left[key], end='')
    print("right: {} entries".format(len(right_entries)))
    if len(right_only) > 0:
        print("These {} stars only appeared in the right list:".format(len(right_only)))
        for key in right:
            print("\t" + right[key], end='')
    if len(left_only) == 0 and len(right_only) == 0:
        print("The lists are identical")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Diff two target lists in memory')
    parser.add_argument('left_list', type=str, help='left side of comparison')
    parser.add_argument('right_list', type=str, help='right side of comparison')
    parser.add_argument('-lr', '--left-ra-col', type=int, required=True, help="Number of column containing RA for the left list (1-indexed)")
    parser.add_argument('-ld', '--left-dec-col', type=int, required=True, help="Number of column containing Dec for the left list (1-indexed)")
    parser.add_argument('-rr', '--right-ra-col', type=int, required=True, help="Number of column containing RA for the right list (1-indexed)")
    parser.add_argument('-rd', '--right-dec-col', type=int, required=True, help="Number of column containing Dec for the right list (1-indexed)")
    parser.add_argument('-p', '--precision', type=int, required=True, help="How many places after the decimal to keep for comparison")
    args = parser.parse_args()
    if args.left_ra_col == args.left_dec_col or args.right_ra_col == args.right_dec_col:
        raise RuntimeError("Same column number supplied for both RA and Dec")
    left = load_list(args.left_list, args.left_ra_col, args.left_dec_col, args.precision)
    right = load_list(args.right_list, args.right_ra_col, args.right_dec_col, args.precision)
    compare_lists(left, right)
