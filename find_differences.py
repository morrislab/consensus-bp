# Find all breakpoints that are in one release but not a second.
from __future__ import print_function, division
import csv
import glob
import os
import sys
from collections import defaultdict

def sort_pos(positions):
  for chrom in positions.keys():
    # True > False, so "P.postype == 'start'" will place starts after ends if
    # they've both at same coordinate.
    # Order by postype: ends, then starts, then undirecteds
    positions[chrom].sort()

def chrom_key(chrom):
  if chrom.isdigit():
    return int(chrom)
  elif chrom == 'X':
    return 100
  elif chrom == 'Y':
    return 101
  else:
    raise Exception('Unknown chrom: %s' % chrom)

def write_bp(bp, outfn):
  with open(outfn, 'w') as outf:
    print('chrom', 'pos', sep='\t', file=outf)
    for chrom in sorted(bp.keys(), key = chrom_key):
      for B in bp[chrom]:
        print(chrom, B, sep='\t', file=outf)

class BpComparer(object):
  def _load(self, bpfn):
    positions = defaultdict(list)

    with open(bpfn) as F:
      reader = csv.DictReader(F, delimiter='\t')
      for row in reader:
        chrom = row['chrom'].upper()
        positions[chrom].append(int(row['pos']))

    return positions

  def _subtract(self, bp_old, bp_new, threshold):
    '''Perform bp_new - bp_old.'''
    matches = []
    avail_old = set(bp_old)
    avail_new = set(bp_new)

    while len(avail_old) > 0 and len(avail_new) > 0:
      closest_dist = float('inf')
      match = None
      for bp_old in avail_old:
        for bp_new in avail_new:
          dist = abs(bp_old - bp_new)
          if dist < closest_dist and dist <= threshold:
            closest_dist = dist
            match = (bp_old, bp_new)

      if match is not None:
        matches.append(match)
        matched_old, matched_new = match
        avail_old.remove(matched_old)
        avail_new.remove(matched_new)
      else:
        break

    return list(avail_new)
    
  def compare(self, bpfn_old, bpfn_new):
    '''Perform new - old.'''
    bp_old = self._load(bpfn_old)
    bp_new = self._load(bpfn_new)
    all_chroms = set(bp_old.keys()) | set(bp_new.keys())
    threshold = 1e4
    remaining_new = {}

    for chrom in all_chroms:
      remaining_new[chrom] = self._subtract(bp_old[chrom], bp_new[chrom], threshold)

    sort_pos(remaining_new)
    return remaining_new

def list_bpfiles(bpdir):
  bpfiles = glob.glob('%s/*-*.txt' % bpdir)
  bpfiles = { os.path.basename(F).split('.')[0]: F for F in bpfiles }
  return bpfiles

def compare(bpfiles_old, bpfiles_new, outdir):
  common = set(bpfiles_old.keys()) & set(bpfiles_new.keys())
  for dataset in common:
    added_bp = BpComparer().compare(bpfiles_old[dataset], bpfiles_new[dataset])
    write_bp(added_bp, os.path.join(outdir, '%s.txt' % dataset))

def main():
  bpfiles_old = list_bpfiles(sys.argv[1])
  bpfiles_new = list_bpfiles(sys.argv[2])
  outdir = sys.argv[3]
  compare(bpfiles_old, bpfiles_new, outdir)

main()
