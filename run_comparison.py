import argparse
import glob
import os
from collections import defaultdict
import glob

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument('--required-methods', dest='required_methods',
    help='Methods that must be present to include entire dataset and each segment within as part of consensus')
  parser.add_argument('--optional-methods', dest='optional_methods',
    help='Methods that will be incorporated if available')
  parser.add_argument('--num-needed-methods', dest='num_needed_methods', type=int, default=3,
    help='Number of available (optional or required) methods necessary to establish consensus')
  parser.add_argument('--window-size', dest='window_size', type=int, default=5000,
    help='Window within which breakpoints must be placed to be considered equivalent')
  parser.add_argument('--undirected', dest='directed', action='store_false',
    help='Whether we should treat breakpoints as directed (i.e., "start" distinct from "end")')
  parser.add_argument('sv_dir', help='Directory containing SVs in VCF format')
  parser.add_argument('method', nargs='+', help='Methods whose CNV calls you wish to use')
  args = parser.parse_args()

  segfiles_by_guid = defaultdict(list)

  for method in args.method:
    if method.endswith('/'):
      method = method[:-1]

    segfiles = glob.glob('%s/*_segments.txt' % method)
    for segfile in segfiles:
      guid = segfile.split('/')[1].split('_')[0]
      segfiles_by_guid[guid].append(method)

  cnt = 0
  for guid, methods_for_guid in segfiles_by_guid.items():
    cnt += 1
    if cnt > 5 and False:
      break

    if len(methods_for_guid) <= args.num_needed_methods:
      continue
    cnv_calls = ' '.join(['%s=%s/%s_segments.txt' % (method, method, guid) for method in methods_for_guid])

    cmd = 'python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/compare_breakpoints.py '
    if args.required_methods:
      cmd += ' --required-methods %s' % args.required_methods
    if args.optional_methods:
      cmd += ' --optional-methods %s' % args.optional_methods
    if not args.directed:
      cmd += ' --undirected'
    cmd += '  --num-needed-methods %s' % args.num_needed_methods
    cmd += '  --window-size %s' % args.window_size

    sv_path = glob.glob(os.path.join(args.sv_dir, '%s.*.sv.vcf.gz' % guid))
    if len(sv_path) != 1:
      continue
    cmd += ' %s %s %s' % (guid, sv_path[0], cnv_calls)

    try:
      print(cmd)
    except IOError, e:
      assert e.errno == 32
      return

main()
