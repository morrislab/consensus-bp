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

        min_sv_size = 0
        if bkdist < min_sv_size:
          continue

        sv[chrom].append(StructVar(chrom=chrom, pos=pos, svclass=svclass))

    return sv

# To get centromeric regions:
#   curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip -c | grep acen
class CentromereParser(object):
  def load(self, filename):
    regions = defaultdict(list)

    with gzip.open(filename) as cf:
      for line in cf:
        fields = line.strip().split()
        chrom, start, end, arm, region_type = fields[0].upper(), int(fields[1]), int(fields[2]), fields[3].lower(), fields[4].lower()

        assert chrom.startswith('CHR')
        chrom = chrom[3:]

        if region_type != 'acen':
          continue
        regions[chrom].append((start, end))

    for chrom in regions.keys():
      assert len(regions[chrom]) == 2
      regions[chrom].sort()
      assert regions[chrom][0][1] == regions[chrom][1][0]
      regions[chrom] = (regions[chrom][0][0], regions[chrom][1][1])
      assert regions[chrom][0] < regions[chrom][1]

    return dict(regions)

class CentromereAndTelomereBreaker(object):
  def __init__(self, threshold = 1e6):
    self._threshold = threshold
    self._chr_lens = {
      '1': 248956422,
      '2': 242193529,
      '3': 198295559,
      '4': 190214555,
      '5': 181538259,
      '6': 170805979,
      '7': 159345973,
      '8': 145138636,
      '9': 138394717,
      '10': 133797422,
      '11': 135086622,
      '12': 133275309,
      '13': 114364328,
      '14': 107043718,
      '15': 101991189,
      '16': 90338345,
      '17': 83257441,
      '18': 80373285,
      '19': 58617616,
      '20': 64444167,
      '21': 46709983,
      '22': 50818468,
      'X': 156040895,
      'Y': 57227415
    }

  def add_breakpoints(self, positions, centromeres):
    # Exclude chrY, since males lack it.
    chroms = set(self._chr_lens.keys()) - set(['Y'])

    for chrom in chroms:
      points = {
        'chrom_start': 1,
        'chrom_end': self._chr_lens[chrom],
        'centromere_start': centromeres[chrom][0],
        'centromere_end': centromeres[chrom][1]
      }
      presence = {k: None for k in points.keys()}

      for pos in positions[chrom]:
        for ptype, point in points.items():
          if presence[ptype] is not None:
            continue
          if abs(point - pos.pos) <= self._threshold:
            presence[ptype] = pos

      for ptype, pos in presence.items():
        if pos is not None:
          log('Already have chr%s.%s at %s' % (chrom, ptype, pos))
          continue
        log('Adding %s.%s at %s' % (chrom, ptype, points[ptype]))
        positions[chrom].append(Position(
          chrom = chrom,
          pos = points[ptype],
          postype = 'undirected',
          method = ptype,
        ))

    sort_pos(positions)

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

class BreakpointClusterer(object):
  def __init__(self, cncalls, window, support_threshold):
    self._window = window
    self._support_threshold = support_threshold

    cncalls = self._combine_cnvs(cncalls)
    self.directed_positions = self._extract_pos(cncalls)
    self._directed_position_indices = self._index_pos(self.directed_positions)

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
      positions[chrom] = [
        Position(chrom=chrom, pos=C[postype], postype=postype, method=C['method'])
        for C in cnvs[chrom]
        for postype in ('start', 'end')
      ]
    sort_pos(positions)

    assert sum([len(P) for P in positions.values()]) == 2*total_cnvs
    return positions

  def _should_use_bp(self, pos, postype):
    if pos.postype != postype:
      return False
    if pos in self.used_bps:
      return False
    return True

  def _find_clusters(self, chrom, postype, support_threshold, start_idx=None, end_idx=None):
    positions = self.directed_positions[chrom]

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
      if not self._should_use_bp(pos, postype):
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
        if self._should_use_bp(positions[nidx], postype):
          next_pos.append(positions[nidx])
        nidx += 1

      cluster = prev_pos + next_pos
      cluster = sorted(cluster, key = lambda P: (P.pos, P.method))
      # Ensure no duplicates.
      assert len(cluster) == len(set(cluster))
      yield cluster
      idx = self._directed_position_indices[cluster[0]]
      prev_pos = []

  def _get_median(self, L):
    # Assumes L is sorted already, since sorting criterion you desire may vary
    # from list to list.
    # If len(L) is even, this will return the lower element (e.g., if len(L) = 4,
    # this returns L[1] rather than L[2]).
    idx = int((len(L) - 0.5) / 2)
    return L[idx]

  def _find_method_representatives(self, cluster):
    partitioned = defaultdict(list)
    for pos in cluster:
      partitioned[pos.method].append(pos)

    representatives = set()
    for method in partitioned.keys():
      representatives.add(partitioned[method][0])

    return sorted(representatives, key = lambda P: (P.pos, P.method))

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
      balancer_found = False
      for reduced_support_threshold in reversed(range(1, self._support_threshold + 1)):
        for balancer in self._find_clusters(
          bad_exemplar.chrom,
          bad_exemplar.expected_postype,
          reduced_support_threshold,
          bad_exemplar.start_idx,
          bad_exemplar.end_idx
        ):
          yield balancer
          balancer_found = True
        if balancer_found:
          break
      #assert balancer_found is True, ('No balancing cluster found for %s' % str(bad_exemplar))

  def _draw_exemplar_from_cluster(self, cluster):
    method_repr = self._find_method_representatives(cluster)
    for member in method_repr:
      self.used_bps[member] = method_repr
    exemplar = self._get_median(method_repr)
    chrom = exemplar.chrom
    self._exemplars[chrom].append(exemplar)

  def _remove_redundant(self):
    for chrom in self._exemplars.keys():
      exemplars = self._exemplars[chrom]

      idx = len(exemplars) - 1
      while idx > 0:
        while exemplars[idx].pos == exemplars[idx - 1].pos:
          # Remove element at idx.
          log('Removing redundant %s because %s' % (exemplars[idx], exemplars[idx - 1]))
          exemplars = exemplars[:idx] + exemplars[(idx + 1):]
          idx -= 1
        idx -= 1

      self._exemplars[chrom] = exemplars

  def _check_sanity(self):
    for chrom in self._exemplars.keys():
      exemplars = self._exemplars[chrom]
      #assert exemplars[0].postype == 'start', ('Exemplar %s on chrom %s is not start' % (exemplars[0], chrom))
      #assert exemplars[-1].postype == 'end', ('Exemplar %s on chrom %s is not end' % (exemplars[-1], chrom))
      for idx in range(len(exemplars) - 1):
        assert exemplars[idx].pos <= exemplars[idx + 1].pos
      # Ideally, I could also check to see whether the breakpoints in between
      # the first and last alternate between opening and closing, but they
      # won't -- as we may add multiple new breakpoints in the case of a tie
      # when balancing, we won't preserve this alternating property.

  def cluster(self):
    self._exemplars = defaultdict(list)
    self.used_bps = {}

    for chrom in self.directed_positions.keys():
      for postype in ('start', 'end'):
        for cluster in self._find_clusters(chrom, postype, self._support_threshold):
          self._draw_exemplar_from_cluster(cluster)

      # The list currently consists of all the start points, followed by all
      # the endpoints. We want to sort by position so we can figure out when we
      # have unclosed starts and unopened ends.
      self._exemplars[chrom] = sorted(self._exemplars[chrom], key = lambda P: P.pos)

      # Remove chromosomes with no breakpoints.
      if len(self._exemplars[chrom]) == 0:
        del self._exemplars[chrom]

    for cluster in self._balance_segments(self._exemplars):
      self._draw_exemplar_from_cluster(cluster)

    sort_pos(self._exemplars)
    # Remove breakpoints at same position -- e.g., because an end and then a
    # start occur on top of each other.
    self._exemplars = BreakpointFilter().remove_adjacent(self._exemplars, threshold=0)
    self._check_sanity()
    return self._exemplars

class StructVarIntegrator(object):
  def __init__(self, sv_filename):
    self._sv = StructVarParser(sv_filename).parse()
    self.matched_sv_to_bp = {}
    self.unmatched_sv = []

  def _find_closest_exemplar(self, sv, exemplars, window):
    closest_exemplar = None
    closest_dist = float('inf')
    for exemplar in exemplars:
      assert exemplar.chrom == sv.chrom
      dist = abs(sv.pos - exemplar.pos)
      if dist < closest_dist and dist <= window:
        closest_dist = dist
        closest_exemplar = exemplar
    return closest_exemplar

  def _match_svs_to_exemplar(self, svs, exemplars, window):
    matches = []
    avail_sv = set(svs)
    avail_exemplars = set(exemplars)

    while len(avail_sv) > 0 and len(avail_exemplars) > 0:
      closest_dist = float('inf')
      match = None

      for sv in avail_sv:
        for exemplar in avail_exemplars:
          dist = abs(sv.pos - exemplar.pos)
          if dist < closest_dist and dist < window:
            closest_dist = dist
            match = (sv, exemplar)

      if match is not None:
        matches.append(match)
        matched_sv, matched_exemplar = match
        avail_sv.remove(matched_sv)
        avail_exemplars.remove(matched_exemplar)
      else:
        break

    return (matches, avail_sv)

  def integrate(self, exemplars, window):
    for chrom in self._sv.keys():
      num_initial_exemplars = len(exemplars[chrom])
      chrom_exemplars = set(exemplars[chrom])
      matches, unmatched_svs = self._match_svs_to_exemplar(self._sv[chrom], chrom_exemplars, window)
      num_moved = 0
      num_added = 0

      for matched_sv, matched_exemplar in matches:
        moved = Position(
          chrom = chrom,
          pos = matched_sv.pos,
          postype = matched_exemplar.postype,
          method = 'sv_%s' % matched_exemplar.method
        )
        chrom_exemplars.remove(matched_exemplar)
        chrom_exemplars.add(moved)
        self.matched_sv_to_bp[moved] = matched_exemplar
        num_moved += 1

      for unmatched_sv in unmatched_svs:
        sv_bp = Position(chrom = chrom, pos = unmatched_sv.pos, postype = 'undirected', method = 'sv')
        if sv_bp in chrom_exemplars:
          # We may have multiple SV breakpoints at the same coordinates (e.g.,
          # a DUP and t2tINV). If this occurs, proceed no further so that our
          # count of the SVs added remains accurate.
          continue
        self.unmatched_sv.append(sv_bp)
        chrom_exemplars.add(sv_bp)
        num_added += 1

      exemplars[chrom] = list(chrom_exemplars)
      assert len(exemplars[chrom]) == num_initial_exemplars + num_added, ('Exemplars = %s, num_initial_exemplars = %s, num_added = %s' % (len(exemplars[chrom]), num_initial_exemplars, num_added))

    sort_pos(exemplars)

class OutputWriter(object):
  def _create_associate(self, bp):
    return {
      'method': bp.method,
      'postype': bp.postype,
      'chrom': bp.chrom,
      'pos': bp.pos,
    }

  def _create_posmap(self, bps, used_bps):
    posmap = defaultdict(lambda: defaultdict(list))
    self._bp_to_entry_map = {}

    for chrom in bps.keys():
      for bp in bps[chrom]:
        entry = {
          'postype': bp.postype,
          'pos': bp.pos,
          'method': bp.method,
          'associates': []
        }

        if bp in used_bps:
          for associate in used_bps[bp]:
            assert associate.chrom == bp.chrom
            assert associate.postype == bp.postype
            entry['associates'].append(self._create_associate(associate))
        else:
          entry['associates'].append(self._create_associate(bp))

        self._bp_to_entry_map[bp] = entry
        posmap[bp.method][chrom].append(entry)

    return posmap

  def _add_svs_to_posmap(self, posmap, used_bps, matched_sv_to_bp, unmatched_sv):
    for usv in unmatched_sv:
      entry = {
        'postype': 'undirected',
        'pos': usv.pos,
        'method': usv.method,
        'associates': []
      }
      posmap['sv'][usv.chrom].append(entry)
      self._bp_to_entry_map[usv] = entry

    for msv, mbp in matched_sv_to_bp.items():
      matched_bps = used_bps[mbp]
      entry = {
        'postype': 'undirected',
        'pos': msv.pos,
        'method': msv.method,
        'associates': [self._create_associate(bp) for bp in matched_bps]
      }
      posmap['sv'][msv.chrom].append(entry)
      self._bp_to_entry_map[msv] = entry

      for bp in matched_bps:
        entry = self._bp_to_entry_map[bp]
        entry['associates'].append({
          'method': 'sv',
          'postype': 'undirected',
          'chrom': msv.chrom,
          'pos': msv.pos,
          })

  def _add_exemplars_to_posmap(self, posmap, exemplars, used_bps, matched_sv_to_bp):
    for chrom in exemplars.keys():
      for exemplar in exemplars[chrom]:
        entry = {
          'postype': exemplar.postype,
          'pos': exemplar.pos,
          'method': exemplar.method,
          'associates': []
        }
        exemplar_associate = self._create_associate(exemplar)
        # Override this, since 'method' field no longer corresponds to
        # something on plot.
        exemplar_associate['method'] = 'consensus'

        if exemplar.method.startswith('sv_'):
          before_move = matched_sv_to_bp[exemplar]
          before_bp = used_bps[before_move]
          after_associates = [self._create_associate(bp) for bp in before_bp] + [exemplar_associate]
          for bp in before_bp:
            self._bp_to_entry_map[bp]['associates'] = after_associates
          entry['associates'] = after_associates
        elif exemplar.method == 'sv':
          entry['associates'].append(self._create_associate(exemplar))
          sv_entry = self._bp_to_entry_map[exemplar]
          sv_entry['associates'].append(exemplar_associate)
        else:
          if exemplar in used_bps:
            associate_bps = used_bps[exemplar]
          else:
            associate_bps = []
          associates = [self._create_associate(bp) for bp in associate_bps] + [exemplar_associate]
          entry['associates'] = associates
          for bp in associate_bps:
            self._bp_to_entry_map[bp]['associates'] = associates

        posmap['consensus'][chrom].append(entry)

  def write_details(self, exemplars, bps, used_bps, matched_sv_to_bp, unmatched_sv, methods, outfn):
    posmap = self._create_posmap(bps, used_bps)
    self._add_svs_to_posmap(posmap, used_bps, matched_sv_to_bp, unmatched_sv)
    self._add_exemplars_to_posmap(posmap, exemplars, used_bps, matched_sv_to_bp)
    with open(outfn, 'w') as outf:
      json.dump({
        'methods': list(methods),
        'bp': posmap,
      }, outf)

  def write_consensus(self, exemplars, outfn):
    with open(outfn, 'w') as consensus_bpf:
      print('chrom', 'pos', sep='\t', file=consensus_bpf)
      for chrom in sorted(exemplars.keys(), key = chrom_key):
        for exemplar in exemplars[chrom]:
          print(chrom, exemplar.pos, sep='\t', file=consensus_bpf)

class BreakpointFilter(object):
  def remove_adjacent(self, breakpoints, threshold):
    filtered = defaultdict(list)

    for chrom in breakpoints.keys():
      points = breakpoints[chrom]

      idx = len(points) - 1
      while idx > 0:
        while points[idx].pos - points[idx - 1].pos <= threshold:
          # Remove element at idx.
          log('Removing %s because of preceding %s' % (points[idx], points[idx - 1]))
          points = points[:idx] + points[(idx + 1):]
          idx -= 1
        idx -= 1

      filtered[chrom] = points

    return filtered

  def remove_sex(self, breakpoints):
    filtered = defaultdict(list)

    for chrom in breakpoints.keys():
      if chrom in ('X', 'Y'):
        continue
      filtered[chrom] = breakpoints[chrom]

    return filtered

def load_cn_calls(cnv_files):
  cn_calls = {}
  for dataset in cnv_files:
    method_name, cnv_fn = dataset.split('=', 1)
    assert method_name not in cn_calls
    cn_calls[method_name] = CnvFileParser().load(cnv_fn)
  method_names = set(cn_calls.keys())

  return (cn_calls, method_names)

def chrom_key(chrom):
  if chrom.isdigit():
    return int(chrom)
  elif chrom == 'X':
    return 100
  elif chrom == 'Y':
    return 101
  else:
    raise Exception('Unknown chrom: %s' % chrom)

def sort_pos(positions):
  for chrom in positions.keys():
    # True > False, so "P.postype == 'start'" will place starts after ends if
    # they've both at same coordinate.
    # Order by postype: ends, then starts, then undirecteds
    positions[chrom].sort(key = lambda P: (P.pos, P.postype == 'undirected', P.postype == 'start', P.method))

def log(*msgs):
  if log.verbose:
    print(*msgs, file=sys.stderr)
log.verbose = False

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
  parser.add_argument('--dataset-name', dest='dataset_name', required=True,
    help='Dataset name')
  parser.add_argument('--sv-filename', dest='sv_fn', required=True,
    help='Consensus structural variants filename (VCF format)')
  parser.add_argument('--centromere-filename', dest='centromere_fn', required=True,
      help='File containing centromeres (e.g., http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz)')
  parser.add_argument('--consensus-bps', dest='consensus_bp_fn', required=True,
    help='Path to consensus BPs output')
  parser.add_argument('--bp-details', dest='bp_details_fn', required=True,
    help='Path to BP details output')
  parser.add_argument('--verbose', dest='verbose', action='store_true')
  parser.add_argument('cnv_files', nargs='+', help='CNV files')
  args = parser.parse_args()

  log.verbose = args.verbose

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

  bc = BreakpointClusterer(cn_calls, args.window_size, args.support_threshold)
  exemplars = bc.cluster()

  svi = StructVarIntegrator(args.sv_fn)
  svi.integrate(exemplars, args.window_size)

  centromere_and_telomere_threshold = 1e6
  adjacent_threshold = 1e4
  centromeres = CentromereParser().load(args.centromere_fn)
  CentromereAndTelomereBreaker(centromere_and_telomere_threshold).add_breakpoints(exemplars, centromeres)
  filtered = BreakpointFilter().remove_adjacent(exemplars, adjacent_threshold)
  filtered = BreakpointFilter().remove_sex(filtered)

  ow = OutputWriter()
  ow.write_consensus(filtered, args.consensus_bp_fn)
  ow.write_details(filtered, bc.directed_positions, bc.used_bps, svi.matched_sv_to_bp, svi.unmatched_sv, consensus_methods, args.bp_details_fn)

if __name__ == '__main__':
  main()
