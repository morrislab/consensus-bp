from __future__ import print_function
import sys
import glob
import os

def parse_list(listfn):
  with open(listfn) as F:
    lines = F.readlines()
  lines = lines[1:]
  samples = [L.strip().split('\t')[0] for L in lines]
  return set(samples)

def fetch_svs(svdir):
  svfn = [os.path.basename(fn).split('.')[0] for fn in glob.glob('%s/*.vcf.gz' % svdir)]
  return set(svfn)

def write_set(S, fn):
  with open(fn, 'w') as F:
    if len(S) == 0:
      return
    print('\n'.join(sorted(S)), file=F)

def parse_released_bp(rlfn):
  samples = set()
  with open(rlfn) as F:
    for line in F:
      fields = line.strip().split('\t')
      sampid = fields[0]
      samples.add(sampid)
  return samples

def main():
  blacklistfn = sys.argv[1]
  greylistfn = sys.argv[2]
  whitelistfn = sys.argv[3]
  svdir = sys.argv[4]
  released_bps_fn = sys.argv[5]
  method_dirs = sys.argv[6:]

  bl = parse_list(blacklistfn)
  gl = parse_list(greylistfn)
  wl = parse_list(whitelistfn)
  desired = wl | gl
  assert len(desired) == 2778
  assert len(desired & bl) == 0

  released = parse_released_bp(released_bps_fn)
  remaining = desired - released
  write_set(remaining, 'missing.consensus.txt')

  svs = fetch_svs(svdir)
  print('Missing SVs:', len(desired - svs))
  print('SVs not on whitelist or greylist:', len(svs - desired))

  samples = {}

  for method_dir in method_dirs:
    method_name = os.path.basename(method_dir)
    if method_name.endswith('/'):
      method_name = method_name[:-1]
    samples[method_name] = glob.glob('%s/*_segments.txt' % method_dir)
    samples[method_name] = set([os.path.basename(M).split('_')[0] for M in samples[method_name]])

  for method in sorted(samples.keys()):
    write_set((remaining & svs) - samples[method], 'missing.%s.has_sv.txt' % method)
    write_set((remaining - svs) - samples[method], 'missing.%s.no_sv.txt' % method)

main()
