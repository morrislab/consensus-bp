from __future__ import print_function, division
import argparse
import glob
import os
from collections import defaultdict
import math
import sys
import csv

def load_sample_list(samplistfn):
  with open(samplistfn) as slist:
    header = next(slist)
    return set([l.strip().split('\t')[0] for l in slist])

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

  cmd = 'python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/make_consensus_breakpoints.py '
  cmd += ' --required-methods %s' % (','.join(methods))
  if len(support_masks) > 0:
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

def run_custom(guid, methods_for_guid, window_size, centromere_fn, sv_dir, base_outdir, should_use, run_name, include_only_chrom=None):
  outdir = os.path.join(base_outdir, 'bp.%s' % run_name)
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

def _process_any3_any2_conservative(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir, run_name, include_only_chrom):
  must_exclude = set(('mustonen095', 'dkfz'))
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
    run_name = run_name,
    include_only_chrom = include_only_chrom,
  )

def process_Y(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir):
  should_use = lambda num_active, active_methods: \
    (num_active >= 2)
  run_custom(
    guid,
    methods_for_guid,
    window_size,
    centromere_fn,
    sv_dir,
    out_dir,
    should_use,
    run_name = 'Y',
    include_only_chrom = 'Y',
  )

def process(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir, sex, mode):
  patient_sex = sex[guid]
  assert patient_sex in ('male', 'female')
  assert mode in ('autosomes', 'X_female', 'X_male', 'Y')

  if mode == 'autosomes':
    _process_any3_any2_conservative(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir, run_name = mode, include_only_chrom=None)
  elif mode == 'X_female' and patient_sex == 'female':
    _process_any3_any2_conservative(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir, run_name = mode, include_only_chrom='X')
  elif mode == 'X_male' and patient_sex == 'male':
    _process_any3_any2_conservative(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir, run_name = mode, include_only_chrom='X')
  elif mode == 'Y' and patient_sex == 'male':
    process_Y(guid, methods_for_guid, window_size, centromere_fn, sv_dir, out_dir)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument('--methods', dest='methods',
    help='Methods over which we compute possible combinations')
  parser.add_argument('--window-size', dest='window_size', type=int, default=5000,
    help='Window within which breakpoints must be placed to be considered equivalent')
  parser.add_argument('--whitelist', dest='whitelist',
    help='List of whitelisted tumor samples')
  parser.add_argument('--greylist', dest='greylist',
    help='List of greylisted tumor samples')
  parser.add_argument('--centromere-filename', dest='centromere_fn', required=True,
    help='File containing centromeres (e.g., http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz)')
  parser.add_argument('--sex-filename', dest='sex_fn', required=True,
    help='File containing inferred sex of each sample from Stefan')
  parser.add_argument('--mode', dest='mode', choices=('autosomes', 'X_female', 'X_male', 'Y'), required=True,
    help='What type of consensus to produce, since different inputs and parameters are typically used for each.')
  parser.add_argument('sv_dir', help='Directory containing SVs in VCF format')
  parser.add_argument('out_dir', help='Output directory')
  parser.add_argument('methods', nargs='+', help='Methods whose CNV calls you wish to use')
  args = parser.parse_args()

  segfiles_by_guid = defaultdict(list)
  valid_samples = load_sample_list(args.whitelist) | load_sample_list(args.greylist)
  assert len(valid_samples) == 2778
  sex = parse_sex(args.sex_fn)

  for method in args.methods:
    if method.endswith('/'):
      method = method[:-1]

    segfiles = glob.glob('%s/*_segments.txt' % method)
    for segfile in segfiles:
      guid = segfile.split('/')[1].split('_')[0]
      segfiles_by_guid[guid].append(method)

  for guid, methods_for_guid in segfiles_by_guid.items():
    if guid not in valid_samples:
      continue
    process(guid, methods_for_guid, args.window_size, args.centromere_fn, args.sv_dir, args.out_dir, sex, args.mode)

if __name__ == '__main__':
  main()
