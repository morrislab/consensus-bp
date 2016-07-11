from __future__ import print_function
import argparse
import csv
from collections import defaultdict, namedtuple
import json
import sys
import gzip
import re

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

Position = namedtuple('Position', ('chrom', 'pos', 'postype', 'method'))
StructVar = namedtuple('StructVar', ('chrom', 'pos', 'svclass'))

class BreakpointScorer(object):
  def __init__(self, cncalls, sv_filename, all_methods, window=5000, directed=True):
    self._cncalls = self._combine_cnvs(cncalls)
    self._all_methods = all_methods
    self._sv = self._load_sv(sv_filename)
    self._sv_class_counts = self._calc_sv_class_counts()

    self._window = window
    self._directed = directed
    self._bp_mutual_scores, self._bp_mutual_scores_by_method = self._calc_breakpoint_mutual_scores()

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

  def _load_sv(self, sv_filename):
    sv = defaultdict(list)

    with gzip.open(sv_filename) as svf:
      for line in svf:
        line = line.strip()
        if line.startswith('#'):
          continue
        fields = line.split('\t')
        chrom, pos, filter, info = fields[0].upper(), int(fields[1]), fields[6], fields[7]
        assert fields[6] == 'PASS'
        svclass = re.compile('SVCLASS=([^;]+)').findall(info)[0]
        sv[chrom].append(StructVar(chrom=chrom, pos=pos, svclass=svclass))

    return sv

  def _calc_sv_class_counts(self):
    sv_class_counts = defaultdict(int)
    for chrom in self._sv.keys():
      for sv in self._sv[chrom]:
        sv_class_counts[sv.svclass] += 1
    return sv_class_counts

  def _extract_pos(self, chrom, cnvs):
    positions = {
      'start': [Position(chrom=chrom, pos=C['start'], postype='start', method=C['method']) for C in cnvs],
      'end': [Position(chrom=chrom, pos=C['end'], postype='end', method=C['method']) for C in cnvs],
    }
    positions['all'] = positions['start'] + positions['end']
    return positions

  def score_breakpoints(self):
    return self._bp_mutual_scores_by_method

  def _calc_breakpoint_mutual_scores(self):
    # method_scores[method][group][score]
    # method \in { broad, dkfz, mustonen095, peifer, vanloo_wedge }
    # group  \in { broad, dkfz, mustonen095, peifer, vanloo_wedge, indeterminate }
    bp_scores_by_method = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    bp_scores = {}

    for chrom, chrom_cnvs in self._cncalls.items():
      positions = self._extract_pos(chrom, chrom_cnvs)
      for postype, pos in positions.items():
        positions[postype] = sorted(pos, key = lambda P: P.pos)

      if self._directed:
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
          open_pos = set([P for P in open_pos if pos.pos - P.pos <= self._window])
          #print(json.dumps((pos, list(open_pos))))
          for P in open_pos:
            support[P].add(pos.method)

      for pos, supporting_methods in support.items():
        S = len(supporting_methods)
        assert pos.method in supporting_methods

        assert pos not in bp_scores
        bp_scores[pos] = S

        if S == 2 and len(self._all_methods) == 5:
          group = supporting_methods - set([pos.method])
        elif S == 4 and len(self._all_methods) == 5:
          group = self._all_methods - supporting_methods
        else:
          group = set(['indeterminate'])

        assert len(group) == 1
        group = group.pop()
        bp_scores_by_method[pos.method][group][S] += 1

    return (bp_scores, bp_scores_by_method)

  def _find_closest_sv(self, bp, chrom, window):
    closest_dist = float('inf')
    closest_sv = None

    for sv in self._sv[chrom]:
      dist = abs(bp.pos - sv.pos)
      if dist < closest_dist and dist <= window:
        closest_sv = sv
        closest_dist = dist

    return closest_sv

  def score_sv(self, required_score=0):
    bp_sv_scores = defaultdict(lambda: defaultdict(int))
    sv_bp_scores = defaultdict(lambda: defaultdict(lambda:  {True: 0, False: 0}))

    for chrom, chrom_cnvs in self._cncalls.items():
      cnv_bp = self._extract_pos(chrom, chrom_cnvs)
      for bp in cnv_bp['all']:
        if self._bp_mutual_scores[bp] < required_score:
          continue
        closest_sv = self._find_closest_sv(bp, chrom, self._window)
        if closest_sv is not None:
          bp_sv_scores[bp.method][closest_sv.svclass] += 1
        else:
          bp_sv_scores[bp.method][None] += 1

    for chrom in self._sv.keys():
      cnv_bp = self._extract_pos(chrom, self._cncalls[chrom])
      for sv in self._sv[chrom]:
        supported_methods = set([
          bp.method for bp in cnv_bp['all'] \
          if abs(bp.pos - sv.pos) <= self._window \
          and self._bp_mutual_scores[bp] >= required_score
        ])
        for method in self._all_methods:
          sv_bp_scores[method][sv.svclass][method in supported_methods] += 1

    return (bp_sv_scores, sv_bp_scores)

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
  parser.add_argument('sv_filename', help='Consensus structural variants filename (VCF format)')
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

  bs = BreakpointScorer(cn_calls, args.sv_filename, consensus_methods, window=args.window_size, directed=args.directed)
  bp_mutual_scores = bs.score_breakpoints()
  bp_sv_scores, sv_bp_scores = bs.score_sv(required_score=4)

  print(json.dumps({
    'dataset': args.dataset_name,
    'bp_mutual_scores': bp_mutual_scores,
    # Number of CNV breakpoints each SV is supported by.
    'sv_bp_scores': sv_bp_scores,
    # Number of SVs each CNV breakpoint is supported by.
    'bp_sv_scores': bp_sv_scores,
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
