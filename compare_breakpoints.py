from __future__ import print_function
import argparse
import csv
from collections import defaultdict, namedtuple
import json
import sys

class CnvFileParser(object):
  def load(self, filename):
    cnvs = []

    with open(filename) as cnvf:
      reader = csv.DictReader(cnvf, delimiter='\t')
      for row in reader:
        try:
          cnv = {
            'chrom': row['chromosome'].upper(),
            'start': int(row['start']),
            # mustonen, at least, produces fully closed intervals -- i.e.,
            # [start, end]. (We know this because some intervals will start and
            # end at the same coordinate.) We want half-open -- i.e., [start,
            # end) because this simplifies the algorithm.
            #
            # DKFZ VCFs sometimes appear half-open and sometimes appear
            # half-closed, depending on the particular interval. After running
            # David's script to create consensus CNV calls, I see intervals both
            # of the type [a, a] (fully closed) and ([a, b], [b, c]) (half open).
            # For now, I treat them as half-open by subtracting one in the
            # conversion script to make them fully closed, then adding one to the
            # end coordinate here (as with all other datasets) to make them
            # half-open again.
            #
            # Note we're assuming that all methods produce fully-closed
            # intervals, which we convert to half-open. This assumption may be
            # wrong.
            'end': int(row['end']) + 1,
            'cell_prev': float(row['clonal_frequency']),
          }
        except ValueError, e:
          print(e, row, filename, file=sys.stderr)
          continue

        total, minor, major = row['copy_number'], row['minor_cn'], row['major_cn']
        # On occasion, ABSOLUTE will report major but not minor.
        if major == 'NA' or minor == 'NA':
          continue
        else:
          # Convert to float first to allow for CN values with decimal point such as "1.0".
          minor = float(minor)
          major = float(major)

        assert minor.is_integer() and major.is_integer()
        # Ensure minor <= major.
        if minor > major:
          minor, major = major, minor
        cnv['minor'] = int(minor)
        cnv['major'] = int(major)

        # 'mustonen*' methods encode X as chr23.
        if cnv['chrom'] == '23':
          cnv['chrom'] = 'X'

        try:
          assert cnv['start'] < cnv['end'], '%s is not < %s in %s' % (cnv['start'], cnv['end'], filename)
        except AssertionError, e:
          print(e, file=sys.stderr)
          continue
        if not (0.0 <= cnv['cell_prev'] <= 1.0):
          print('Cellular prevalence is %s in %s' % (cnv['cell_prev'], filename), file=sys.stderr)
          if cnv['cell_prev'] < 0 or cnv['cell_prev'] >= 1.05:
            continue
          cnv['cell_prev'] = 1.0
        cnvs.append(cnv)

    return cnvs

class CnvOrganizer(object):
  def __init__(self, cnvs):
    cnvs = self._organize_cnvs(cnvs)
    cnvs = self._eliminate_duplicate_breakpoints(cnvs)
    self._ensure_no_overlap(cnvs)
    self._cnvs = cnvs

  def _organize_cnvs(self, cnvs):
    organized = defaultdict(list)

    for cnv in cnvs:
      chrom = cnv['chrom']
      del cnv['chrom']
      organized[chrom].append(cnv)

    for chrom, chrom_cnvs in organized.items():
      chrom_cnvs.sort(key = lambda C: C['start'])

    return organized

  def _eliminate_duplicate_breakpoints(self, cnvs):
    filtered = defaultdict(list)

    for chrom, chrom_cnvs in cnvs.items():
      for cnv in chrom_cnvs:
        if len(filtered[chrom]) > 0 and filtered[chrom][-1]['start'] == cnv['start'] and filtered[chrom][-1]['end'] == cnv['end']:
          continue
        filtered[chrom].append(cnv)

    return filtered

  def _ensure_no_overlap(self, cnvs):
    for chrom, chrom_cnvs in cnvs.items():
      for idx in range(len(chrom_cnvs) - 1):
        current, next = chrom_cnvs[idx], chrom_cnvs[idx + 1]
        assert current['start'] < current['end'] <= next['start'] < next['end']

  def organize(self):
    return self._cnvs

class BreakpointScorer(object):
  def __init__(self, cncalls, all_methods):
    self._cncalls = self._combine_cnvs(cncalls)
    self._all_methods = all_methods

  def _combine_cnvs(self, cncalls):
    all_chroms = set([C for L in cncalls.values() for C in L.keys()])
    all_cnvs = defaultdict(list)

    for chrom in all_chroms:
      chrom_cnvs = []
      for method in cncalls.keys():
        for cnv in cncalls[method][chrom]:
          cnv = dict(cnv) # Copy object before modifying
          cnv['method'] = method
          all_cnvs[chrom].append(cnv)

    return all_cnvs

  def score(self, window=5000, directed=True):
    # method_scores[method][group][score]
    # method \in { broad, dkfz, mustonen095, peifer, vanloo_wedge }
    # group  \in { broad, dkfz, mustonen095, peifer, vanloo_wedge }
    method_scores = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    Position = namedtuple('Position', ('pos', 'postype', 'method'))

    for chrom, chrom_cnvs in self._cncalls.items():
      positions = {
        'start': [Position(pos=C['start'], postype='start', method=C['method']) for C in chrom_cnvs],
        'end': [Position(pos=C['end'], postype='end', method=C['method']) for C in chrom_cnvs],
      }
      positions['all'] = positions['start'] + positions['end']
      for postype, pos in positions.items():
        positions[postype] = sorted(pos, key = lambda P: P.pos)

      if directed:
        keys = ('start', 'end')
      else:
        keys = ('all',)

      support = defaultdict(set)

      for postype in keys:
        # Ensure no duplicates.
        assert len(positions[postype]) == len(set(positions[postype]))
        open_pos = set()

        # Determine how many other breakpoints are close to current one.
        for pos in positions[postype]:
          open_pos.add(pos)
          for P in open_pos:
            assert pos.pos >= P.pos
          open_pos = set([P for P in open_pos if pos.pos - P.pos <= window])
          #print(json.dumps((pos, list(open_pos))))
          for P in open_pos:
            support[P].add(pos.method)

      for pos, supporting_methods in support.items():
        S = len(supporting_methods)
        assert pos.method in supporting_methods

        if S == 2 and len(self._all_methods) == 5:
          group = supporting_methods - set([pos.method])
        elif S == 4 and len(self._all_methods) == 5:
          group = self._all_methods - supporting_methods
        else:
          group = set(['indeterminate'])

        assert len(group) == 1
        group = group.pop()
        method_scores[pos.method][group][S] += 1

    return method_scores

def load_cn_calls(cnv_files):
  cn_calls = {}
  for dataset in cnv_files:
    method_name, cnv_fn = dataset.split('=', 1)
    assert method_name not in cn_calls
    cn_calls[method_name] = CnvFileParser().load(cnv_fn)
  method_names = set(cn_calls.keys())

  return (cn_calls, method_names)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
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
  parser.add_argument('dataset_name', help='Dataset name')
  parser.add_argument('cnv_files', nargs='+', help='CNV files')
  args = parser.parse_args()

  dataset_name = args.dataset_name

  if args.required_methods is not None:
    required_methods = set(args.required_methods.split(','))
  else:
    required_methods = set()
  if args.optional_methods is not None:
    optional_methods = set(args.optional_methods.split(','))
  else:
    optional_methods = set()

  cn_calls, avail_methods = load_cn_calls(args.cnv_files)
  consensus_methods = (required_methods | optional_methods) & avail_methods
  # Ensure we have no extraneous methods.
  assert consensus_methods == avail_methods

  assert avail_methods.issuperset(required_methods) and len(consensus_methods) >= args.num_needed_methods

  for method, cnvs in cn_calls.items():
    cn_calls[method] = CnvOrganizer(cnvs).organize()

  bs = BreakpointScorer(cn_calls, consensus_methods)
  scores = bs.score(window=args.window_size, directed=args.directed)

  print(json.dumps({
    'dataset': args.dataset_name,
    'scores': scores,
  }))
  return

  import numpy as np
  for M, MS in scores.items():
    scores[M] = np.mean([T for S, C in MS.items() for T in C*[S]])
  print(json.dumps({
    'dataset': args.dataset_name,
    'scores': scores,
  }))

if __name__ == '__main__':
  main()
