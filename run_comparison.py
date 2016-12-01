from __future__ import print_function, division
import argparse
import glob
import os
from collections import defaultdict
import glob
import math
import sys

def load_blacklist(blacklistfn):
  with open(blacklistfn) as blist:
    return set([l.strip() for l in blist.readlines()])

def generate_method_combos(methods):
  methods = sorted(set(methods))
  N = len(methods)
  for idx in range(1, 2**N):
    active_mask = bin(idx)[2:] # Remove '0b' prefix

    active_mask = active_mask.zfill(N)
    active_methods = [methods[I] for I, C in enumerate(active_mask) if C == '1']
    yield (active_mask, active_methods)

def generate_command(methods, guid, window_size, centromere_fn, sv_dir, out_dir, support_masks):
  cnv_calls = ' '.join(['%s=%s/%s_segments.txt' % (method, method, guid) for method in sorted(methods)])
  assert len(support_masks) > 0

  cmd = 'python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/make_consensus_breakpoints.py '
  cmd += ' --required-methods %s' % (','.join(methods))
  cmd += ' --support-masks %s' % ','.join(support_masks)
  cmd += ' --num-needed-methods %s' % len(methods)
  cmd += ' --dataset-name %s' % guid

  sv_path = glob.glob(os.path.join(sv_dir, '%s.*.sv.vcf.gz' % guid))
  assert len(sv_path) <= 1
  if len(sv_path) == 1:
    cmd += ' --sv-filename %s' % sv_path[0]
  cmd += ' --window-size %s' % window_size
  cmd += ' --centromere-filename %s' % centromere_fn
  cmd += ' --consensus-bps %s/%s.txt' % (out_dir, guid)
  cmd += ' --bp-details %s/%s.json' % (out_dir, guid)
  cmd += ' --verbose'
  cmd += ' %s' % cnv_calls
  cmd += ' >%s/%s.stdout' % (out_dir, guid)
  cmd += ' 2>%s/%s.stderr' % (out_dir, guid)
  return cmd

def print_safely(S):
  # Account for EOF generated on STDOUT when output piped to "head" or
  # whatever.
  try:
    if S is None:
      return
    print(S)
  except IOError, e:
    assert e.errno == 32
    sys.exit()

def run_custom(guid, methods_for_guid, window_size, centromere_fn, sv_dir, outdir, should_use):
  '''Combine any3 with any2, in which the any2 are seleected from a certain set of "reliable" methods.'''
  if not os.path.exists(outdir):
    os.makedirs(outdir)

  masks = set()

  for active_mask, active_methods in generate_method_combos(methods_for_guid):
    num_active = sum([int(x) for x in active_mask])
    if should_use(num_active, active_methods):
      masks.add(active_mask)

  cmd = generate_command(
    methods_for_guid,
    guid,
    window_size,
    centromere_fn,
    sv_dir,
    outdir,
    masks,
  )
  print_safely(cmd)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument('--methods', dest='methods',
    help='Methods over which we compute possible combinations')
  parser.add_argument('--window-size', dest='window_size', type=int, default=5000,
    help='Window within which breakpoints must be placed to be considered equivalent')
  parser.add_argument('--blacklist', dest='blacklist',
    help='List of blacklisted tumor samples')
  parser.add_argument('--centromere-filename', dest='centromere_fn', required=True,
      help='File containing centromeres (e.g., http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz)')
  parser.add_argument('--num-needed-methods', dest='num_needed_methods', type=int,
    help='Number of available (optional or required) methods necessary to establish consensus')
  parser.add_argument('sv_dir', help='Directory containing SVs in VCF format')
  parser.add_argument('out_dir', help='Output directory')
  parser.add_argument('methods', nargs='+', help='Methods whose CNV calls you wish to use')
  args = parser.parse_args()

  segfiles_by_guid = defaultdict(list)
  blacklist = load_blacklist(args.blacklist)

  for method in args.methods:
    if method.endswith('/'):
      method = method[:-1]

    segfiles = glob.glob('%s/*_segments.txt' % method)
    for segfile in segfiles:
      guid = segfile.split('/')[1].split('_')[0]
      segfiles_by_guid[guid].append(method)

  for guid, methods_for_guid in segfiles_by_guid.items():
    if guid in blacklist:
      continue
    if args.num_needed_methods is not None and len(methods_for_guid) < args.num_needed_methods:
      continue

    must_exclude = set(('mustonen095', 'dkfz'))
    should_use = lambda num_active, active_methods: \
      (num_active >= 3) or \
      (num_active == 2 and len(set(active_methods) & must_exclude) == 0)
    run_custom(
      guid,
      methods_for_guid,
      args.window_size,
      args.centromere_fn,
      args.sv_dir,
      args.out_dir,
      should_use
    )

main()
