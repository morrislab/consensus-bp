import argparse
import glob
import os
from collections import defaultdict
import glob
import math

def load_blacklist(blacklistfn):
  with open(blacklistfn) as blist:
    return set([l.strip() for l in blist.readlines()])

def generate_method_combos(methods):
  methods = sorted(set(methods))
  N = len(methods)
  for idx in range(1, 2**N):
    active_mask = bin(idx)[2:] # Remove '0b' prefix

    # Skip cases of only one method, since they won't produce any breakpoints.
    # Our consensus BP algorithm places intervals around each breakpoint, then
    # finds intersections between those intervals. A lone method won't have any
    # intersections with just itself, and so won't produce any BPs except for
    # SVs and centromere/telomere boundaries.
    if sum([int(x) for x in active_mask]) == 1:
      continue
    if sum([int(x) for x in active_mask]) != 4:
      continue

    active_mask = active_mask.zfill(N)
    active_methods = [methods[I] for I, C in enumerate(active_mask) if C == '1']
    yield (active_mask, active_methods)

def generate_command(methods, guid, window_size, centromere_fn, sv_dir, out_dir, support_threshold=None):
  cnv_calls = ' '.join(['%s=%s/%s_segments.txt' % (method, method, guid) for method in methods])

  cmd = 'python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/make_consensus_breakpoints.py '
  cmd += ' --required-methods %s' % (','.join(methods))
  if support_threshold is None:
    cmd += ' --support-threshold %s' % len(methods)
    cmd += ' --num-needed-methods %s' % len(methods)
  else:
    cmd += ' --support-threshold %s' % support_threshold
    cmd += ' --num-needed-methods %s' % support_threshold
  cmd += ' --dataset-name %s' % guid

  sv_path = glob.glob(os.path.join(sv_dir, '%s.*.sv.vcf.gz' % guid))
  if len(sv_path) != 1:
    return None
  cmd += ' --window-size %s' % window_size
  cmd += ' --sv-filename %s' % sv_path[0]
  cmd += ' --centromere-filename %s' % centromere_fn
  cmd += ' --consensus-bps %s/%s.txt' % (out_dir, guid)
  cmd += ' --bp-details %s/%s.json' % (out_dir, guid)
  cmd += ' --verbose'
  cmd += ' %s' % cnv_calls
  cmd += ' >%s/%s.stdout' % (out_dir, guid)
  cmd += ' 2>%s/%s.stderr' % (out_dir, guid)
  return cmd

def run_specific_methods(guid, methods_for_guid, window_size, centromere_fn, sv_dir, base_outdir):
  for active_mask, active_methods in generate_method_combos(methods_for_guid):
    outdir = os.path.join(base_outdir, 'methods.%s' % active_mask)
    if not os.path.exists(outdir):
      os.makedirs(outdir)
    cmd = generate_command(
      active_methods,
      guid,
      window_size,
      centromere_fn,
      sv_dir,
      outdir,
    )
    if cmd is None:
      continue
    try:
      print(cmd)
    except IOError, e:
      assert e.errno == 32
      return

def run_any_combo(guid, methods_for_guid, window_size, centromere_fn, sv_dir, base_outdir):
  for support_threshold in range(2, len(methods_for_guid) + 1):
    outdir = os.path.join(base_outdir, 'methods.any%s' % support_threshold)
    if not os.path.exists(outdir):
      os.makedirs(outdir)
    cmd = generate_command(
      methods_for_guid,
      guid,
      window_size,
      centromere_fn,
      sv_dir,
      outdir,
      support_threshold = support_threshold
    )
    if cmd is None:
      continue
    try:
      print(cmd)
    except IOError, e:
      assert e.errno == 32
      return

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
  parser.add_argument('sv_dir', help='Directory containing SVs in VCF format')
  parser.add_argument('out_dir', help='Output directory')
  parser.add_argument('method', nargs='+', help='Methods whose CNV calls you wish to use')
  args = parser.parse_args()

  methods = args.methods.split(',')

  segfiles_by_guid = defaultdict(list)
  blacklist = load_blacklist(args.blacklist)

  for method in methods:
    if method.endswith('/'):
      method = method[:-1]

    segfiles = glob.glob('%s/*_segments.txt' % method)
    for segfile in segfiles:
      guid = segfile.split('/')[1].split('_')[0]
      segfiles_by_guid[guid].append(method)

  for guid, methods_for_guid in segfiles_by_guid.items():
    if guid in blacklist:
      continue
    if set(methods_for_guid) != set(methods):
      continue

    run_specific_methods(guid, methods_for_guid, args.window_size, args.centromere_fn, args.sv_dir, args.out_dir)
    #run_any_combo(guid, methods_for_guid, args.window_size, args.centromere_fn, args.sv_dir, args.out_dir)

main()
