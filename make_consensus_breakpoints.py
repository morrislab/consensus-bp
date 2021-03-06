from __future__ import print_function
import argparse
import csv
import gzip
import sys
import json
from collections import defaultdict, namedtuple

StructVar = namedtuple('StructVar', ('chrom', 'pos', 'svclass'))

# Note that `associates` are distinct from `breakpoints.
#   * Breakpoints: Position objects that lie in the interval.
#   * Associates: Position objects that may not be in the interval, but are
#     associated with it in some way. I use this to track all the breakpoints
#     that originally lay in the different intervals that were subsequently
#     intersected. It can thus provide insight into what different breakpoints
#     gave rise to the intervals whose intersection yielded a consensus
#     breakpoint. This is especially useful in cases when the breakpoints that
#     generated an intersection all lie outside the interval.
Interval = namedtuple('Interval', ('start', 'end', 'breakpoints', 'associates', 'method'))

class Position(object):
  # Position was originally a namedtuple, but we want it to be mutable.
  def __init__(self, chrom, pos, postype, method):
    self.chrom = chrom
    self.pos = pos
    self.postype = postype
    self.method = method
    self._interval = None

  @property
  def interval(self):
    return self._interval

  @interval.setter
  def interval(self, val):
    self._interval = val

  def __str__(self):
    return 'chr%s(%s, %s, %s)' % (self.chrom, self.pos, self.postype, self.method)

  def __repr__(self):
    return str(self)


# Taken from https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes.
CHROM_LENS = {
  '1':249250621,
  '2':243199373,
  '3':198022430,
  '4':191154276,
  '5':180915260,
  '6':171115067,
  '7':159138663,
  '8':146364022,
  '9':141213431,
  '10':135534747,
  '11':135006516,
  '12':133851895,
  '13':115169878,
  '14':107349540,
  '15':102531392,
  '16':90354753,
  '17':81195210,
  '18':78077248,
  '19':59128983,
  '20':63025520,
  '21':48129895,
  '22':51304566,
  'X':155270560,
  'Y':59373566,
}

class CnvFileParser(object):
  def load(self, filename):
    cnvs = []

    with open(filename) as cnvf:
      reader = csv.DictReader(cnvf, delimiter='\t')
      for row in reader:
        try:
          cnv = {
            'chrom': row['chromosome'].upper(),
            'start': int(float(row['start'])),
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
            'end': int(float(row['end'])) + 1,
          }
        except ValueError, e:
          print(e, row, filename, file=sys.stderr)
          continue

        # 'mustonen*' methods encode X as chr23.
        if cnv['chrom'] == '23':
          cnv['chrom'] = 'X'
        elif cnv['chrom'] == '24':
          cnv['chrom'] = 'Y'
        elif cnv['chrom'] == 'M':
          continue

        # DKFZ intervals are left-closed, right-open, as are Jabba intervals.
        if 'dkfz/' in filename or 'jabba/' in filename:
          cnv['end'] -= 1

        try:
          assert cnv['start'] < cnv['end'], '%s is not < %s in %s' % (cnv['start'], cnv['end'], filename)
          assert cnv['start'] >= 1, ('Start position %s is < 1 in %s ' % (cnv['start'], filename))
          assert cnv['end'] <= CHROM_LENS[cnv['chrom']], ('End position %s exceeds length of chromosome %s (%s) in %s' % (cnv['end'], cnv['chrom'], CHROM_LENS[cnv['chrom']], filename))
        except AssertionError, e:
          print(e, file=sys.stderr)
          continue
        cnvs.append(cnv)

    return cnvs

class StructVarParser(object):
  def __init__(self, sv_filename):
    self._sv_filename = sv_filename

  def parse(self):
    sv = defaultdict(list)
    used_pos = set()

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

        if (chrom, pos) not in used_pos:
          sv[chrom].append(StructVar(chrom=chrom, pos=pos, svclass=svclass))
          used_pos.add((chrom, pos))

    for chrom in sv.keys():
      sv[chrom].sort(key = lambda S: S.pos)

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
    self.inside_centromere = 0

  def _add(self, positions, centromeres):
    chroms = set(CHROM_LENS.keys())

    for chrom in chroms:
      centromere_start, centromere_end = centromeres[chrom]
      points = {
        'chrom_start': 1,
        'chrom_end': CHROM_LENS[chrom],
        'centromere_start': centromere_start,
        'centromere_end': centromere_end,
      }
      presence = {k: None for k in points.keys()}

      if chrom in positions.keys():
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

  def _move_invalid(self, positions, centromeres):
    for chrom in positions.keys():
      centromere_start, centromere_end = centromeres[chrom]

      for pos in positions[chrom]:
        if pos.pos < 1:
          log('%s is before chromosome start. Moving to 1.' % pos)
          pos.pos = 1
        elif pos.pos > CHROM_LENS[chrom]:
          log('%s is after chromosome end. Moving to %s.' % CHROM_LENS[chrom])
          pos.pos = CHROM_LENS[chrom]
        elif centromere_start < pos.pos < centromere_end:
          closest_dist = float('inf')
          cloeset_point = None
          for point in (centromere_start, centromere_end):
            dist = abs(pos.pos - point)
            if dist < closest_dist:
              closest_dist = dist
              closest_point = point
          assert closest_point is not None
          log('%s is inside centromere (%s, %s). Moving to %s.' % (pos, centromere_start, centromere_end, closest_point))
          self.inside_centromere += 1
          pos.pos = closest_point

  def add_breakpoints(self, positions, centromeres, associate_tracker):
    self._associate_tracker = associate_tracker
    self._move_invalid(positions, centromeres)
    self._add(positions, centromeres)
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

class ConsensusMaker(object):
  def __init__(self, cncalls, window, support_methods, associate_tracker):
    self._cncalls = cncalls
    self._window = window
    self._support_methods = support_methods
    self.cna_pos = defaultdict(dict)
    self._associate_tracker = associate_tracker

  def _sort_intervals(self, intervals):
    return sorted(intervals, key = lambda I: I.start)

  def _make_intervals(self, cncalls, upstream_window, downstream_window):
    all_chroms = set([C for L in cncalls.values() for C in L.keys()])
    intervals = {}

    for chrom in all_chroms:
      intervals[chrom] = []

      for method in cncalls.keys():
        positions = []
        for cnv in cncalls[method][chrom]:
          assert cnv['start'] <= cnv['end']
          positions += [
            Position(chrom=chrom, pos=cnv[coordtype], postype=coordtype, method=method)
            for coordtype in ('start', 'end')
          ]
        for idx in range(len(positions) - 1):
          assert positions[idx].pos <= positions[idx + 1].pos
        intervals[chrom] += self._make_chrom_intervals(positions, chrom, upstream_window, downstream_window)
        self.cna_pos[method][chrom] = positions

      intervals[chrom] = self._sort_intervals(intervals[chrom])

    return intervals

  def _score_breakpoints(self, intervals):
    # Score breakpoints based on the size of the gap they surround. Smaller
    # gaps indicate more certainty in breakpoint placement and thus are
    # preferable.
    method_intervals = defaultdict(list)
    bp_scores = {}

    for chrom in intervals.keys():
      for interval in intervals[chrom]:
        method_intervals[interval.method].append(interval)

    for method in method_intervals.keys():
      sorted_intervals = sorted(method_intervals[method], key = lambda I: I.end - I.start)
      for interval_rank, interval in enumerate(sorted_intervals):
        interval_score = 1.0 - (float(interval_rank) / len(sorted_intervals))
        for bp in interval.breakpoints:
          bp_scores[bp] = interval_score

    return bp_scores

  def _make_chrom_intervals(self, positions, chrom, upstream_window, downstream_window):
    if len(positions) == 0:
      return []

    assert len(positions) > 0 and len(positions) % 2 == 0
    assert positions[0].postype == 'start'
    assert positions[-1].postype == 'end'

    chrom_len = CHROM_LENS[chrom]

    terminal_intervals = set()
    for terminal_pos in (positions[0], positions[-1]):
      terminal_intervals.add(Interval(
        start = max(1, terminal_pos.pos - upstream_window),
        end = min(chrom_len, terminal_pos.pos + downstream_window),
        # Use frozenset rather than set so we can later form sets of intervals.
        # (Sets require their elements to be hashable, but sets themselves
        # aren't hashable.)
        breakpoints = frozenset([terminal_pos]),
        associates = None,
        method = terminal_pos.method,
      ))

    intervals = set()
    for idx in range(1, len(positions) - 1, 2):
      upstream_end, downstream_start = positions[idx], positions[idx + 1]
      assert upstream_end.postype == 'end' and downstream_start.postype == 'start'
      assert upstream_end.pos <= downstream_start.pos
      assert upstream_end.method == downstream_start.method and upstream_end.chrom == downstream_start.chrom == chrom
      intervals.add(Interval(
        start = max(1, upstream_end.pos - upstream_window),
        end = min(chrom_len, downstream_start.pos + downstream_window),
        breakpoints = frozenset([upstream_end, downstream_start]),
        associates = None,
        method = upstream_end.method # Which is same as downstream_start.method
      ))

    all_intervals = terminal_intervals | intervals
    assert len(all_intervals) > 0
    # One interval for every pair in middle, plus separate intervals for first and last point.
    assert len(all_intervals) == (len(positions) - 2)/2 + 2
    all_intervals = sorted(all_intervals, key = lambda I: (I.start, I.end, I.method))

    for interval in all_intervals:
      assert interval.start <= interval.end
    for idx in range(len(all_intervals) - 1):
      first, second = all_intervals[idx], all_intervals[idx + 1]
      assert first.start <= second.start
      assert first.end <= second.end

    return all_intervals

  def _find_first_intersecting_intervals(self, intervals, support_methods):
    # As we halt after we find the first intersection, we don't need to handle
    # resetting these values in the loop after we find an intersection. This is
    # why I wrote both _intersect and _find_first_intersection.
    prev_overlapping = []
    prev_num_overlapping = 0
    threshold_reached = False
    prev_match_better = False
    idx = 0
    overlapping = []

    while idx < len(intervals):
      interval = intervals[idx]
      overlapping = [
        I for I in overlapping
        if I.end > interval.start
        # If segment is smaller than window size, a method can have intervals
        # that overlap with its own intervals, even if none of the segments in
        # the input overlap with each other. If this occurs, discard upstream
        # segments to accommodate the one we're dealing with at the moment.
        and I.method != interval.method
      ]

      overlapping.append(interval)
      num_overlapping = len(overlapping)
      assert num_overlapping > 0

      # Ensure no duplicates.
      assert len(set(overlapping)) == num_overlapping
      # Ensure no method overlaps itself.
      #print(*(['whoa'] + overlapping), sep='\n')
      #print(*(['whoa'] + intervals), sep='\n')
      overlapping_methods = set([I.method for I in overlapping])
      assert len(overlapping_methods) == num_overlapping

      if threshold_reached and num_overlapping < prev_num_overlapping:
        #print('decline', overlapping, overlapping_methods, set([I.method for I in prev_overlapping]))
        prev_match_better = True

      for M in support_methods:
        if M.issubset(overlapping_methods):
          #print('reached', overlapping_methods, M)
          threshold_reached = True
          break
      else:
        if threshold_reached:
          #print('nothing matched', overlapping_methods, set([I.method for I in prev_overlapping]))
          prev_match_better = True

      if prev_match_better:
        assert threshold_reached
        # One or more methods have dropped out of the intersection, so it's time to report the intersection.
        # Suppose you encounter intervals from methods in this order: method_A,
        # method_B, method_C, method_A, method_D. If we report when
        # num_overlapping drops *or stays the same*, we miss the interval from
        # method_D, which gives us smaller intersection and more confidence.
        # Thus, require num_overlapping to drop before reporting.
        for I in prev_overlapping:
          intervals.remove(I)
          assert I not in intervals # Ensure not multiple copies of I in intervals
        #print('taking', prev_overlapping)
        return prev_overlapping

      # Duplicate list.
      #print('prev_overlapping', prev_num_overlapping, prev_overlapping)
      #print('overlapping', num_overlapping, overlapping)
      prev_overlapping = list(overlapping)
      prev_num_overlapping = num_overlapping
      idx += 1

    return None

  def _find_intersecting_intervals(self, intervals, support_methods):
    # Duplicate, as we will be modifying the list.
    intervals = list(intervals)

    assert len(intervals) > 0
    # Ensure no duplicates.
    assert len(set(intervals)) == len(intervals)
    # Ensure no duplicate breakpoints.
    breakpoints = [bp for I in intervals for bp in I.breakpoints]
    assert len(set(breakpoints)) == len(breakpoints)

    for idx in range(len(intervals) - 1):
      assert intervals[idx].start <= intervals[idx + 1].start

    while True:
      intersection = self._find_first_intersecting_intervals(intervals, support_methods)
      if intersection is None:
        return
      else:
        yield intersection

  def _compute_intersection(self, intersecting):
    assert len(intersecting) > 0
    sorted_by_start = sorted(intersecting, key = lambda I: I.start)
    sorted_by_end   = sorted(intersecting, key = lambda I: I.end)
    intersect_start = sorted_by_start[-1].start
    intersect_end = sorted_by_end[0].end

    # BUG: this line formerly used I.start and I.end, rather than
    # intersect_start and intersect_end, as the bounds. The effect of
    # this mistake was that the breakpoints lost some precision, as we
    # effectively select each interval's representative from the *union* of
    # intervals associated with the input breakpoints, not the *intersection*.
    #
    # E.g.: suppose you have input BPs at 1 bp, 60 kb, 60.479 kb, 67.5 kb, and
    # 92.527 kb.  The associated intervals are [1, 50 kb], [10 kb, 60 kb],
    # [10.479 kb, 60.479 kb], [17.5 kb, 67.5 kb], [42.527 kb, 92.527 kb]. Their
    # intersection is [42.527 kb, 50 kb]. None of the breakpoints lie in this
    # interval, so we should take the 5' end (42.527 kb) as the representative.
    # Instead, we took "1 bp" as the representative before the fix.
    #
    # Note that this can also affect intervals that *do* have some of their
    # associated breakpoints lying in the intersection, since the pool of
    # breakpoints from which we select the interval's representative was
    # associated with the union, not the intersection.
    breakpoints = [bp for I in intersecting for bp in I.breakpoints if intersect_start <= bp.pos <= intersect_end]

    associates = [bp for I in intersecting for bp in I.breakpoints]
    # Clone Positions in associates so we can see exactly what positions
    # supported each consensus breakpoint, without them later being
    # repositioned in accordance with SVs.
    associates = [Position(bp.chrom, bp.pos, bp.postype, bp.method) for bp in associates]

    return Interval(
      start = intersect_start,
      end = intersect_end,
      breakpoints = frozenset(breakpoints),
      associates = frozenset(associates),
      method = 'intersection',
    )

  def _make_consensus(self, intervals, bp_scores, support_methods):
    consensus = defaultdict(list)

    assert len(support_methods) > 0

    for chrom in intervals.keys():
      for intersecting_intervals in self._find_intersecting_intervals(intervals[chrom], support_methods):
        intersection = self._compute_intersection(intersecting_intervals)
        associated_pos = [bp for I in intersecting_intervals for bp in I.breakpoints]

        if len(intersection.breakpoints) > 0:
          # Take breakpoint with highest score, which corresponds to smallest gap in associated interval.
          consensus_bp = sorted(intersection.breakpoints, key = lambda bp: (bp_scores[bp], bp.postype == 'start'))[-1]
        else:
          # We have no breakpoints in the interval, so we just take the 5' end of the interval.
          consensus_bp = Position(
            chrom = chrom,
            pos = intersection.start,
            postype = 'undirected',
            method = 'added_from_intersection'
          )
        consensus_bp.interval = intersection
        assert consensus_bp not in consensus[chrom]
        consensus[chrom].append(consensus_bp)

        for apos in associated_pos:
          self._associate_tracker.add(consensus_bp, apos)

      consensus[chrom].sort(key = lambda P: P.pos)
      for idx in range(len(consensus[chrom]) - 1):
        assert consensus[chrom][idx].pos <= consensus[chrom][idx + 1].pos

    return consensus

  def _make_single_method_consensus(self, cncalls, method):
    consensus = defaultdict(list)

    for chrom in cncalls[method].keys():
      positions = []
      for cnv in cncalls[method][chrom]:
        assert cnv['start'] <= cnv['end']
        consensus[chrom] += [
          Position(chrom=chrom, pos=cnv[coordtype], postype=coordtype, method=method)
          for coordtype in ('start', 'end')
        ]
      consensus[chrom].sort(key = lambda P: P.pos)
      for idx in range(len(consensus[chrom]) - 1):
        assert consensus[chrom][idx].pos <= consensus[chrom][idx + 1].pos
    return consensus

  def make_consensus(self):
    if len(self._support_methods) == 0:
      return defaultdict(list)
    elif len(self._support_methods) == 1 and len(list(self._support_methods)[0]) == 1:
      method = list(list(self._support_methods)[0])[0]
      return self._make_single_method_consensus(self._cncalls, method)

    intervals = self._make_intervals(self._cncalls, int(0.5*self._window), int(0.5*self._window))
    bp_scores = self._score_breakpoints(intervals)
    return self._make_consensus(intervals, bp_scores, self._support_methods)

class StructVarIntegrator(object):
  def __init__(self, sv_filename, associate_tracker):
    self._sv = StructVarParser(sv_filename).parse()
    self._associate_tracker = associate_tracker

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
      if chrom not in exemplars:
        exemplars[chrom] = []
      num_initial_exemplars = len(exemplars[chrom])
      chrom_exemplars = set(exemplars[chrom])
      matches, unmatched_svs = self._match_svs_to_exemplar(self._sv[chrom], chrom_exemplars, window)
      num_moved = 0
      num_added = 0

      for matched_sv, matched_exemplar in matches:
        matched_exemplar.pos = matched_sv.pos
        matched_exemplar.method = 'sv_%s' % matched_exemplar.method
        num_moved += 1

      for unmatched_sv in unmatched_svs:
        sv_bp = Position(chrom = chrom, pos = unmatched_sv.pos, postype = 'undirected', method = 'sv')
        if sv_bp in chrom_exemplars:
          # We may have multiple SV breakpoints at the same coordinates (e.g.,
          # a DUP and t2tINV). If this occurs, proceed no further so that our
          # count of the SVs added remains accurate.
          continue
        chrom_exemplars.add(sv_bp)
        num_added += 1

      exemplars[chrom] = list(chrom_exemplars)
      #assert len(exemplars[chrom]) == num_initial_exemplars + num_added, ('Exemplars = %s, num_initial_exemplars = %s, num_added = %s' % (len(exemplars[chrom]), num_initial_exemplars, num_added))

    sort_pos(exemplars)

class AssociateTracker(object):
  def __init__(self):
    self._ass = defaultdict(set)

  def _keyify(self, pos):
    return '_'.join([str(v) for v in (pos.method, pos.postype, pos.chrom, pos.pos)])

  def add(self, pos1, pos2):
    if pos1 == pos2:
      # Silently ignore attempts to associate a position with itself.
      return
    self._ass[pos1].add(pos2)
    self._ass[pos2].add(pos1)

  def remove(self, pos):
    for other in self._ass[pos]:
      self._ass[other].remove(pos)
    del self._ass[pos]

  def get_associates(self):
    stringified = {}
    for pos in self._ass.keys():
      stringified[self._keyify(pos)] = [self._keyify(P) for P in self._ass[pos]]
    return stringified

class OutputWriter(object):
  def __init__(self):
    self._posmap = defaultdict(lambda: defaultdict(list))

  def _add_positions_to_posmap(self, method, positions):
    _make_entry = lambda pos: {
          'postype': pos.postype,
          'pos': pos.pos,
          'method': pos.method,
    }

    for chrom in positions.keys():
      for pos in positions[chrom]:
        entry = _make_entry(pos)
        if method == 'consensus' and pos.interval is not None:
          # Must check if method is consensus, as Positions are modified in
          # place -- most consensus positions also belong to their individual
          # methods.
          entry['interval'] = [pos.interval.start, pos.interval.end],
          associates = [P for P in pos.interval.associates if P != pos]
          associates = [_make_entry(P) for P in associates]
          entry['associates'] = sorted(associates, key = lambda P: (P['pos'], P['method']))
        self._posmap[method][chrom].append(entry)

  def write_details(self, consensus, cna_pos, consensus_methods, avail_methods, stats, params, associate_tracker, outfn):
    self._add_positions_to_posmap('consensus', consensus)
    for method, pos in cna_pos.items():
      self._add_positions_to_posmap(method, pos)

    with open(outfn, 'w') as outf:
      json.dump({
        'consensus_methods': sorted(consensus_methods),
        'avail_methods': sorted(avail_methods),
        'params': params,
        'bp': self._posmap,
        'stats': stats,
        'associates': associate_tracker.get_associates(),
      }, outf)

  def write_consensus(self, exemplars, outfn):
    with open(outfn, 'w') as consensus_bpf:
      print('chrom', 'pos', sep='\t', file=consensus_bpf)
      for chrom in sorted(exemplars.keys(), key = chrom_key):
        for exemplar in exemplars[chrom]:
          print(chrom, exemplar.pos, sep='\t', file=consensus_bpf)

class BreakpointFilter(object):
  def remove_proximal(self, breakpoints, threshold, associate_tracker):
    filtered = defaultdict(list)
    removed = set()

    for chrom in breakpoints.keys():
      points = breakpoints[chrom]

      idx = len(points) - 1
      while idx > 0:
        while idx > 0 and points[idx].pos - points[idx - 1].pos <= threshold:
          # Remove element at idx.
          removed.add(points[idx])
          log('Removing %s because of preceding %s' % (points[idx], points[idx - 1]))
          associate_tracker.remove(points[idx])
          points = points[:idx] + points[(idx + 1):]
          idx -= 1
        idx -= 1

      filtered[chrom] = points

    return (filtered, removed)

  def _retain_only_chroms(self, breakpoints, should_keep):
    retained = defaultdict(list)

    for chrom in breakpoints.keys():
      if should_keep(chrom):
        retained[chrom] = breakpoints[chrom]

    return retained

  def retain_only_chrom(self, breakpoints, to_keep):
    return self._retain_only_chroms(breakpoints, lambda C: C == to_keep)

  def retain_only_autosomes(self, breakpoints):
    return self._retain_only_chroms(breakpoints, lambda C: C not in ('X', 'Y'))


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

def check_sanity(breakpoints, proximity_threshold):
  for chrom in breakpoints.keys():
    last_bp = None
    for idx in range(len(breakpoints[chrom]) - 1):
      assert breakpoints[chrom][idx].pos < breakpoints[chrom][idx + 1].pos
      assert breakpoints[chrom][idx + 1].pos - breakpoints[chrom][idx].pos > proximity_threshold

def count_bp(bp):
  return sum([len(V) for V in bp.values()])

def log(*msgs):
  if log.verbose:
    print(*msgs, file=sys.stderr)
log.verbose = False

def generate_supported_methods(support_masks, consensus_methods):
  if support_masks is None:
    return set()

  support_masks = list(set(support_masks.split(',')))
  S = sorted(consensus_methods)
  support_methods = set()

  for mask in support_masks:
    support = set()
    assert len(mask) == len(consensus_methods)
    for idx, bit in enumerate(mask):
      if bit == '1':
        support.add(S[idx])
    support_methods.add(frozenset(support))

  return support_methods

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
  parser.add_argument('--support-masks', dest='support_masks',
    help='Binary masks indicating which methods must support a breakpoint to place a consensus breakpoint')
  parser.add_argument('--dataset-name', dest='dataset_name', required=True,
    help='Dataset name')
  parser.add_argument('--sv-filename', dest='sv_fn',
    help='Consensus structural variants filename (VCF format)')
  parser.add_argument('--centromere-filename', dest='centromere_fn', required=True,
      help='File containing centromeres (e.g., http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz)')
  parser.add_argument('--include-only-chrom', dest='include_only_chrom',
    help='Include only a certain chromosome. As we use separate strategies for autosomes, X on females, and X & Y on males, we can perform separate runs for each using this option.')
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

  centromere_and_telomere_threshold = 1e6
  # Low proximity threshold just removes breakpoints that are far too close to be real.
  proximity_threshold = 10
  stats = {}

  for method, cnvs in cn_calls.items():
    cn_calls[method] = CnvOrganizer(cnvs).organize()

  associate_tracker = AssociateTracker()
  support_methods = generate_supported_methods(args.support_masks, consensus_methods)
  # In JSON output, if we're not actually using any methods, we shouldn't list
  # any as available.
  if len(support_methods) == 0:
    consensus_methods = set()

  cm = ConsensusMaker(cn_calls, args.window_size, support_methods, associate_tracker)
  consensus = cm.make_consensus()

  if args.sv_fn is not None:
    svi = StructVarIntegrator(args.sv_fn, associate_tracker)
    svi.integrate(consensus, args.window_size)

  centromeres = CentromereParser().load(args.centromere_fn)
  ctb = CentromereAndTelomereBreaker(centromere_and_telomere_threshold)
  ctb.add_breakpoints(consensus, centromeres, associate_tracker)

  if args.include_only_chrom is None:
    consensus = BreakpointFilter().retain_only_autosomes(consensus)
  else:
    consensus = BreakpointFilter().retain_only_chrom(consensus, args.include_only_chrom)

  stats['inside_centromere'] = ctb.inside_centromere
  stats['before_removing_proximal'] = count_bp(consensus)
  consensus, removed = BreakpointFilter().remove_proximal(consensus, proximity_threshold, associate_tracker)
  #stats['proximal_removed'] = list(removed)
  stats['after_removing_proximal'] = count_bp(consensus)

  check_sanity(consensus, proximity_threshold)

  params = {
    'num_needed_methods': args.num_needed_methods,
    'support_masks': args.support_masks,
    'support_methods': [list(E) for E in support_methods],
    'centromere_and_telomere_threshold': centromere_and_telomere_threshold,
    'proximity_threshold': proximity_threshold,
    'window_size': args.window_size,
  }

  ow = OutputWriter()
  ow.write_consensus(consensus, args.consensus_bp_fn)
  ow.write_details(consensus, cm.cna_pos, consensus_methods, avail_methods, stats, params, associate_tracker, args.bp_details_fn)

if __name__ == '__main__':
  main()
