from __future__ import print_function, division
import argparse
import glob
import os
from collections import defaultdict
import glob
import math
import sys
import csv

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

def generate_command(methods, guid, window_size, centromere_fn, sv_dir, out_dir, support_masks, include_only_chrom=None):
  cnv_calls = ' '.join(['%s=%s/%s_segments.txt' % (method, method, guid) for method in sorted(methods)])
  assert len(support_masks) > 0

  cmd = 'python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/make_consensus_breakpoints.py '
  cmd += ' --required-methods %s' % (','.join(methods))
  cmd += ' --support-masks %s' % ','.join(support_masks)
  cmd += ' --num-needed-methods %s' % len(methods)
  cmd += ' --dataset-name %s' % guid

  sv_path = glob.glob(os.path.join(sv_dir, '%s.*.sv.vcf.gz' % guid))
  if len(sv_path) != 1:
    return None
  cmd += ' --window-size %s' % window_size
  cmd += ' --sv-filename %s' % sv_path[0]
  cmd += ' --centromere-filename %s' % centromere_fn
  cmd += ' --consensus-bps %s/%s.txt' % (out_dir, guid)
  cmd += ' --bp-details %s/%s.json' % (out_dir, guid)
  if include_only_chrom is not None:
    cmd += ' --include-only-chrom %s' % include_only_chrom
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

def nCr(n,r):
  f = math.factorial
  return int(f(n) / f(r) / f(n-r))

def run_individual_methods(guid, methods_for_guid, window_size, centromere_fn, sv_dir, base_outdir, out_prefix='', include_only_chrom=None):
  for active_mask, active_methods in generate_method_combos(methods_for_guid):
    if len(active_methods) != 1:
      continue
    outdir = os.path.join(base_outdir, 'methods.%s%s' % (out_prefix, active_methods[0]))
    if not os.path.exists(outdir):
      os.makedirs(outdir)
    cmd = generate_command(
      methods_for_guid,
      guid,
      window_size,
      centromere_fn,
      sv_dir,
      outdir,
      (active_mask,),
      include_only_chrom
    )
    print_safely(cmd)

def run_any_N(guid, methods_for_guid, window_size, centromere_fn, sv_dir, base_outdir, N, exclude=None):
  if exclude is None:
    exclude = set()
  methods_for_guid = set(methods_for_guid) - exclude
  assert N <= len(methods_for_guid)

  outdir = os.path.join(base_outdir, 'methods.any%s' % N)
  if len(exclude) > 0:
    outdir += '_exclude_' + ','.join(sorted(exclude))
  if not os.path.exists(outdir):
    os.makedirs(outdir)

  masks = set()

  for active_mask, active_methods in generate_method_combos(methods_for_guid):
    num_active = sum([int(x) for x in active_mask])
    if num_active < N:
      continue
    masks.add(active_mask)
  M = len(methods_for_guid)
  expected = sum([nCr(M, i) for i in range(N, M + 1)])
  assert len(masks) == expected

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

def run_custom(guid, methods_for_guid, window_size, centromere_fn, sv_dir, base_outdir, should_use, run_name, include_only_chrom=None):
  outdir = os.path.join(base_outdir, 'methods.%s' % run_name)
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
    include_only_chrom = include_only_chrom
  )
  print_safely(cmd)

def parse_sex(sexfn):
  sex = {}
  with open(sexfn) as F:
    reader = csv.DictReader(F, delimiter='\t')
    for row in reader:
      assert row['pred_gender'] in ('male', 'female')
      sex[row['tumourid']] = row['pred_gender']
  return sex

def run_autosome_strats(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir):
  run_individual_methods(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir)

  for exclude in (tuple(), ('dkfz',), ('peifer',), ('dkfz', 'peifer')):
    for N in range(2, len(methods_for_guid) + 1):
      if N > len(methods_for_guid) - len(exclude):
        continue
      run_any_N(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir, N, exclude=set(exclude))

  for must_include_one in (
    set(('broad', 'peifer')),
    set(('broad', 'vanloo_wedge_segs')),
    set(('peifer', 'vanloo_wedge_segs')),
    set(('broad', 'peifer', 'vanloo_wedge_segs')),
  ):
    assert must_include_one.issubset(methods_for_guid)
    should_use = lambda num_active, active_methods: \
      (num_active >= 3) or \
      (num_active == 2 and len(set(active_methods) & must_include_one) > 0)
    run_custom(
      guid,
      methods_for_guid,
      window_size,
      centromere_fn,
      sv_dir,
      out_dir,
      should_use,
      'any3_any2_include_%s' % ','.join(must_include_one)
    )

  for must_exclude in (
    set(('mustonen095',)),
    set(('mustonen095', 'dkfz')),
    set(('mustonen095', 'dkfz', 'peifer'))
  ):
    assert must_exclude.issubset(methods_for_guid)
    should_use = lambda num_active, active_methods: \
      (num_active >= 3) or \
      (num_active == 2 and len(set(active_methods) & must_exclude) == 0)
    run_custom(
      guid,
      methods_for_guid,
      window_size,
      centromere_fn,
      sv_dir,
      out_dir,
      should_use,
      'any3_any2_exclude_%s' % ','.join(must_exclude)
    )

  # Run any3_any2_conservative
  conservative_methods = set(('broad', 'peifer', 'vanloo_wedge_segs'))
  should_use = lambda num_active, active_methods: \
    (num_active >= 3) or \
    (num_active == 2 and set(active_methods).issubset(conservative_methods))
  run_custom(
    guid,
    methods_for_guid,
    window_size,
    centromere_fn,
    sv_dir,
    out_dir,
    should_use,
    'any3_any2_conservative'
  )

def run_X_strats(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir, sample_sex):
  prefix = 'X_%s' % sample_sex

  should_use = lambda num_active, active_methods: \
    (num_active >= 3) or \
    (num_active == 2 and len(set(('dkfz', 'mustonen095')) & set(active_methods)) == 0)
  run_custom(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir, should_use, '%s_any3_any2_except_dkfz,mustonen095' % prefix, include_only_chrom='X')

  should_use = lambda num_active, active_methods: \
    (num_active >= 3)
  run_custom(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir, should_use, '%s_any3' % prefix, include_only_chrom='X')

  should_use = lambda num_active, active_methods: \
    (num_active >= 2)
  run_custom(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir, should_use, '%s_any2' % prefix, include_only_chrom='X')

  run_individual_methods(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir, out_prefix='X_%s_' % sample_sex, include_only_chrom='X')

def run_Y_strats(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir):
  assert set(methods_for_guid) == set(('dkfz', 'jabba'))
  should_use = lambda num_active, active_methods: \
    (num_active >= 2)
  run_custom(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir, should_use, 'Y_dkfz,jabba', include_only_chrom='Y')

  run_individual_methods(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir, out_prefix='Y_', include_only_chrom='Y')

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
  parser.add_argument('--sex-filename', dest='sex_fn', required=True,
      help='File containing inferred sex of each sample from Stefan')
  parser.add_argument('sv_dir', help='Directory containing SVs in VCF format')
  parser.add_argument('out_dir', help='Output directory')
  parser.add_argument('methods', nargs='+', help='Methods whose CNV calls you wish to use')
  args = parser.parse_args()

  methods = args.methods.split(',')

  segfiles_by_guid = defaultdict(list)
  blacklist = load_blacklist(args.blacklist)
  sex = parse_sex(args.sex_fn)

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
    #run_autosome_strats(guid, methods_for_guid, args.window_size, args.centromere_fn, args.sv_dir, args.out_dir)

    sample_sex = sex[guid]
    if sample_sex == 'male':
      # Exclude Battenberg's calls on males, as its calls are poor quality
      # because of haploid state and consequent lack of heterozygous SNPs.
      methods_for_guid = [M for M in methods_for_guid if not M.startswith('vanloo_wedge')]

      # Peifer doesn't usually provide calls for X. This confirms so.
      assert 'peifer' not in methods
      run_X_strats(guid, methods_for_guid, args.window_size, args.centromere_fn, args.sv_dir, args.out_dir, sample_sex)

      # Only Jabba and DKFZ provide Y calls.
      methods_for_guid = list(set(methods_for_guid) & set(('dkfz', 'jabba')))
      run_Y_strats(guid, methods_for_guid, args.window_size, args.centromere_fn, args.sv_dir, args.out_dir)
    elif sample_sex == 'female':
      run_X_strats(guid, methods_for_guid, args.window_size, args.centromere_fn, args.sv_dir, args.out_dir, sample_sex)
    else:
      raise Exception('Unknown sex: %s' % sample_sex)

main()
