from __future__ import print_function
import argparse
import csv
import gzip
import sys
import json
from collections import defaultdict, namedtuple

Position = namedtuple('Position', ('chrom', 'pos', 'postype', 'method'))
StructVar = namedtuple('StructVar', ('chrom', 'pos', 'svclass'))
MissingExemplarInterval = namedtuple('MissingExemplarInterval', ('chrom', 'start_idx', 'end_idx', 'expected_postype'))

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

class StructVarParser(object):
  def __init__(self, sv_filename):
    self._sv_filename = sv_filename

  def parse(self):
    sv = defaultdict(list)

    with gzip.open(self._sv_filename) as svf:
      for line in svf:
        line = line.strip()
        if line.startswith('#'):
          continue
        fields = line.split('\t')
        chrom, pos, filter, info_raw = fields[0].upper(), int(fields[1]), fields[6], fields[7]
        assert fields[6] == 'PASS'

        info = {}
        for tokens in info_raw.split(';'):
          if '=' in tokens:
            K, V = tokens.split('=', 1)
          else:
            K, V = tokens, None
          info[K] = V

        svclass = info['SVCLASS']
        if 'BKDIST' in info:
          bkdist = int(info['BKDIST'])
        else:
          # -1 signifies translocation, so we use -2 to indicate missing data.
          bkdist = -2

        min_sv_size = 10000
        if bkdist < min_sv_size:
          continue

        sv[chrom].append(StructVar(chrom=chrom, pos=pos, svclass=svclass))

    return sv

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

class ConsensusMaker(object):
  def __init__(self, cncalls, window, support_threshold):
    self._window = window
    self._support_threshold = support_threshold

    cncalls = self._combine_cnvs(cncalls)
    self._directed_positions = self._extract_pos(cncalls)
    self._directed_position_indices = self._index_pos(self._directed_positions)

  def _index_pos(self, positions):
    indices = {}
    total_pos = 0

    for chrom in positions.keys():
      pos = positions[chrom]
      indices.update({ P: idx for (idx, P) in enumerate(pos) })
      total_pos += len(pos)

    # Ensure no duplicate positions.
    assert len(indices) == total_pos
    return indices

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

  def _extract_pos(self, cnvs):
    positions = {}
    total_cnvs = 0

    for chrom in cnvs.keys():
      total_cnvs += len(cnvs[chrom])
      pos = [
        Position(chrom=chrom, pos=C[postype], postype=postype, method=C['method'])
        for C in cnvs[chrom]
        for postype in ('start', 'end')
      ]
      # True > False, so "P.postype == 'end'" will place ends after starts if
      # they've both at same coordinate.
      positions[chrom] = sorted(pos, key = lambda P: (P.pos, P.postype == 'end', P.method))

    assert sum([len(P) for P in positions.values()]) == 2*total_cnvs
    return positions

  def _find_clusters(self, chrom, postype, support_threshold, start_idx=None, end_idx=None):
    positions = self._directed_positions[chrom]

    if start_idx is None:
      start_idx = 0
    if end_idx is None:
      end_idx = len(positions)
    assert 0 <= start_idx < end_idx <= len(positions)

    # Ensure no duplicates.
    assert len(positions) == len(set(positions))
    prev_pos = []

    idx = start_idx
    while idx < end_idx:
      pos = positions[idx]
      if pos.postype != postype:
        idx += 1
        continue

      prev_pos.append(pos)
      for P in prev_pos:
        assert pos.pos >= P.pos
      prev_pos = [P for P in prev_pos if pos.pos - P.pos <= self._window]
      supporting_methods = set([P.method for P in prev_pos])

      if len(supporting_methods) < support_threshold:
        idx += 1
        continue

      # If we reach this point, we know we've found a cluster member. It and
      # all breakpoints preceding it within the window will be in prev_pos.
      # Now, we must find the following positions, which we store in next_pos.
      nidx = idx + 1
      remaining_window = pos.pos - prev_pos[0].pos
      assert 0 <= remaining_window <= self._window
      next_pos = []
      while nidx < end_idx and positions[nidx].pos - pos.pos <= remaining_window:
        assert positions[nidx].pos >= pos.pos
        if positions[nidx].postype == postype:
          next_pos.append(positions[nidx])
        nidx += 1

      cluster = prev_pos + next_pos
      cluster = sorted(cluster, key = lambda P: (P.pos, P.method))
      # Ensure no duplicates.
      assert len(cluster) == len(set(cluster))
      yield cluster
      idx = nidx
      prev_pos = []

  def _get_median(self, L):
    # Assumes L is sorted already, since sorting criterion you desire may vary
    # from list to list.
    # If len(L) is even, this will return the lower element (e.g., if len(L) = 4,
    # this returns L[1] rather than L[2]).
    idx = int((len(L) - 0.5) / 2)
    return L[idx]

  def _find_method_medians(self, cluster):
    partitioned = defaultdict(list)
    for pos in cluster:
      partitioned[pos.method].append(pos)

    medians = set()
    for method in partitioned.keys():
      medians.add(self._get_median(partitioned[method]))

    return sorted(medians, key = lambda P: (P.pos, P.method))

  def _balance_segments(self, exemplars):
    # Fix positions so we have no unopened ends or unclosed starts.
    bad_exemplars = set()

    for chrom in exemplars.keys():
      first_exemplar, last_exemplar = exemplars[chrom][0], exemplars[chrom][-1]
      if first_exemplar.postype != 'start':
        bad_exemplars.add(MissingExemplarInterval(
          chrom = chrom,
          start_idx = None,
          end_idx = self._directed_position_indices[first_exemplar],
          expected_postype = 'start'
        ))

      if last_exemplar.postype != 'end':
        bad_exemplars.add(MissingExemplarInterval(
          chrom = chrom,
          start_idx = self._directed_position_indices[last_exemplar] + 1,
          end_idx = None,
          expected_postype = 'end'
        ))

      last_postype = first_exemplar.postype
      last_idx = self._directed_position_indices[first_exemplar]
      for exemplar in exemplars[chrom][1:-1]:
        expected_postype = (last_postype == 'start' and 'end' or 'start')
        if exemplar.postype != expected_postype:
          bad_exemplars.add(MissingExemplarInterval(
            chrom = chrom,
            start_idx = last_idx + 1,
            end_idx = self._directed_position_indices[exemplar],
            expected_postype = expected_postype
          ))
        last_postype = exemplar.postype
        last_idx = self._directed_position_indices[exemplar]

    for bad_exemplar in bad_exemplars:
      for reduced_support_threshold in reversed(range(1, self._support_threshold + 1)):
        balancing_clusters = list(self._find_clusters(
          bad_exemplar.chrom,
          bad_exemplar.expected_postype,
          reduced_support_threshold,
          bad_exemplar.start_idx,
          bad_exemplar.end_idx
        ))
        print(json.dumps([bad_exemplar, reduced_support_threshold, balancing_clusters]))

  def make_consensus(self):
    self._exemplars = defaultdict(list)
    self._directed_clusters = {}

    for chrom in self._directed_positions.keys():
      for postype in ('start', 'end'):
        for cluster in self._find_clusters(chrom, postype, self._support_threshold):
          method_medians = self._find_method_medians(cluster)
          exemplar = self._get_median(method_medians)
          self._exemplars[chrom].append(exemplar)
          self._directed_clusters[exemplar] = cluster

      # The list currently consists of all the start points, followed by all
      # the endpoints. We want to sort by position so we can figure out when we
      # have unclosed starts and unopened ends.
      self._exemplars[chrom] = sorted(self._exemplars[chrom], key = lambda P: P.pos)

      # Remove chromosomes with no breakpoints.
      if len(self._exemplars[chrom]) == 0:
        del self._exemplars[chrom]

    self._balance_segments(self._exemplars)

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
  parser.add_argument('--support-threshold', dest='support_threshold', type=int, default=4,
    help='Number of methods that must support a cluster to place a consensus breakpoint within it')
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

  cm = ConsensusMaker(cn_calls, args.window_size, args.support_threshold)
  cm.make_consensus()

if __name__ == '__main__':
  main()
