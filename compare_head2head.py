from __future__ import print_function, division
import sys
import os
import json
from make_consensus_breakpoints import CentromereParser, CHROM_LENS
from collections import defaultdict, namedtuple
import glob

Position = namedtuple('Position', ('pos', 'method'))

def exclude_near(breakpoints, around, threshold):
  retained = defaultdict(list)

  for chrom in breakpoints.keys():
    for bp in breakpoints[chrom]:
      closest_dist = float('inf')
      closest_other = None
      for other in around[chrom]:
        dist = abs(bp.pos - other.pos)
        if dist < closest_dist:
          closest_dist = dist
          closest_other = other
      if closest_dist > threshold:
        retained[chrom].append(bp)

  return retained

def extract(breakpoints, prefix, truth):
  retained = defaultdict(list)
  for chrom in breakpoints.keys():
    for bp in breakpoints[chrom]:
      if bp.method.startswith(prefix) is truth:
        retained[chrom].append(bp)
  return retained

def count_bp(bp):
  return sum([len(V) for V in bp.values()])

def parse_centromeres_and_telomeres(centromeres):
  cents_and_telos = defaultdict(list)

  # Exclude chrY, since males lack it.
  chroms = set(CHROM_LENS.keys()) - set(['Y'])
  for chrom in chroms:
    points = {
      'chrom_start': 1,
      'chrom_end': CHROM_LENS[chrom],
      'centromere_start': centromeres[chrom][0],
      'centromere_end': centromeres[chrom][1]
    }
    for K, P in points.items():
      cents_and_telos[chrom].append(Position(pos=P, method=K))
    cents_and_telos[chrom].sort(key = lambda p: p.pos)

  return cents_and_telos

def determine_methods(method_combo, methods):
  assert len(method_combo) == len(methods)
  methods = sorted(methods)
  cmp_methods, consensus_methods = set(), set()

  for indicator, method in zip(method_combo, methods):
    if indicator == '1':
      consensus_methods.add(method)
    elif indicator == '0':
      cmp_methods.add(method)
    else:
      raise Exception('Unknown indicator in %s: %s' % (method_combo, indicator))

  assert len(cmp_methods) == 2
  assert len(consensus_methods) == 4
  return (frozenset(cmp_methods), frozenset(consensus_methods))

def load_breakpoints(fn):
  bp = {}
  with open(fn) as F:
    J = json.load(F)
    for method in J['bp'].keys():
      bp[method] = {}
      for chrom in J['bp'][method]:
        # Note that some methods (e.g., peifer) report duplicate positions (e.g.,
        # start & end coordinate at same position). The set will eliminate these.
        bp[method][chrom] = list(set([Position(pos=B['pos'], method=B['method']) for B in J['bp'][method][chrom]]))
    sort_pos(bp[method])
  return bp

class Matcher(object):
  def _match(self, A, B, max_dist):
    matched = set()
    avail_A = set(A)
    avail_B = set(B)

    while len(avail_A) > 0 and len(avail_B) > 0:
      closest_dist = float('inf')
      closest_A = None
      closest_B = None

      for a in avail_A:
        for b in avail_B:
          dist = abs(a.pos - b.pos)
          if dist < closest_dist:
            closest_A = a
            closest_B = b
            closest_dist = dist

      assert closest_A is not None and closest_B is not None
      if closest_dist > max_dist:
        break

      avail_A.remove(closest_A)
      avail_B.remove(closest_B)
      matched.add((closest_A, closest_B))
    return matched

  def subtract(self, A, B, dist_threshold):
    '''Perform subtraction of A - B based on which positions are mutually
    closer than `dist_threshold`.'''
    difference = {}
    for chrom in A.keys():
      if chrom not in B.keys():
        # Duplicate list.
        difference[chrom] = list(A[chrom])
        continue

      matches = self._match(A[chrom], B[chrom], dist_threshold)
      matched_from_A = [M[0] for M in matches]
      difference[chrom] = list(set(A[chrom]) - set(matched_from_A))
      assert len(difference[chrom]) == len(A[chrom]) - len(matched_from_A)

    sort_pos(difference)
    return difference

  def intersect(self, A, B, dist_threshold):
    '''Perform intersection of A \intersect B based on which positions are
    mutually closer than `dist_threshold`. Return only the points from A.'''
    intersection = {}
    for chrom in A.keys():
      if chrom not in B.keys():
        continue

      matches = self._match(A[chrom], B[chrom], dist_threshold)
      intersection[chrom] = [M[0] for M in matches]

    sort_pos(intersection)
    return intersection

def sort_pos(positions):
  for chrom in positions.keys():
    positions[chrom].sort(key = lambda P: (P.pos, P.method))

def prepare_bp(bp, cents_and_telos, cent_telo_threshold):
  stats = {'all_bp': count_bp(bp)}
  bp = exclude_near(bp, cents_and_telos, cent_telo_threshold)
  stats['away_from_cents_and_telos'] = count_bp(bp)
  bp = extract(bp, 'sv', False)
  stats['away_from_sv_and_cents_and_telos'] = count_bp(bp)
  return (bp, stats)

def compare_to_consensus(consensus_bp, cmp_bp, all_methods):
  # Ensure no centromeres, telomeres, or SVs remain.
  for chrom in cmp_bp.keys():
    for bp in cmp_bp[chrom]:
      assert bp.method in all_methods

  consensus_threshold = 1e5
  # TP: non-SV BPs that appear in both comparison and consensus set
  # FP: non-SV BPs that appear in comparison set but not consensus set
  # FN: non-SV BPs that appear in consensus set but not comparison set
  tp = Matcher().intersect(consensus_bp, cmp_bp, consensus_threshold)
  fp = Matcher().subtract(cmp_bp, consensus_bp, consensus_threshold)
  fn = Matcher().subtract(consensus_bp, cmp_bp, consensus_threshold)
  assert count_bp(tp) + count_bp(fp) == count_bp(cmp_bp)

  stats = {}
  stats['tp'] = count_bp(tp)
  stats['fp'] = count_bp(fp)
  stats['fn'] = count_bp(fn)
  return stats

def compare_methods(bp, guid, cmp_methods, consensus_methods, all_methods, cents_and_telos, indiv_bps):
  assert set(bp.keys()) == consensus_methods | set(['consensus'])
  assert consensus_methods | cmp_methods == all_methods
  best_precision_score = float('-inf')
  best_recall_score = float('-inf')
  most_bp_count = -1
  best_precision_method = None
  best_recall_method = None
  most_bp_method = None

  cent_telo_threshold = 1e6
  stats = {}
  consensus_bp, stats['consensus'] = prepare_bp(bp['consensus'], cents_and_telos, cent_telo_threshold)

  for cmp_method in cmp_methods:
    cmp_bp, stats[cmp_method] = prepare_bp(indiv_bps[cmp_method], cents_and_telos, cent_telo_threshold)

    cmp_stats = compare_to_consensus(consensus_bp, cmp_bp, all_methods)
    stats[cmp_method].update(cmp_stats)

    try:
      precision = cmp_stats['tp'] / float(cmp_stats['tp'] + cmp_stats['fp'])
      recall = cmp_stats['tp'] / float(cmp_stats['tp'] + cmp_stats['fn'])
      stats[cmp_method]['precision'] = precision
      stats[cmp_method]['recall'] = recall

      if precision > best_precision_score:
        best_precision_score = precision
        best_precision_method = cmp_method
      if recall > best_recall_score:
        best_recall_score = recall
        best_recall_method = cmp_method
      if stats[cmp_method]['away_from_sv_and_cents_and_telos'] > most_bp_count:
        most_bp_count = stats[cmp_method]['away_from_sv_and_cents_and_telos']
        most_bp_method = cmp_method
    except ZeroDivisionError:
      stats[cmp_method]['precision'] = None
      stats[cmp_method]['recall'] = None
      return (stats, None, None, None)

  assert 0 <= best_precision_score <= 1
  assert 0 <= best_recall_score <= 1
  assert most_bp_count >= 0
  return (stats, best_precision_method, best_recall_method, most_bp_method)

def load_runs(runs, all_methods):
  results = defaultdict(dict)
  indiv_bps = defaultdict(dict)

  for run in runs:
    resultfns = glob.glob('%s/*-*.json' % run)
    for R in resultfns:
      guid = os.path.basename(R).split('.')[0]

      dirname = os.path.dirname(R)
      assert dirname.startswith('methods.')
      method_vector = dirname.split('.')[1]
      cmp_methods, consensus_methods = determine_methods(method_vector, all_methods)

      results[R] = {
        'guid': guid,
        'cmp_methods': cmp_methods,
        'consensus_methods': consensus_methods,
      }

      bp = load_breakpoints(R)
      for method in consensus_methods:
        if method in indiv_bps[guid].keys():
          continue
        indiv_bps[guid][method] = bp[method]

  return (results, indiv_bps)

def determine_winners(results, all_methods, cents_and_telos, indiv_bps):
  precision_winners = defaultdict(lambda: defaultdict(int))
  recall_winners = defaultdict(lambda: defaultdict(int))
  # TOTAL DOMINATION!!!
  total_dom_winners = defaultdict(lambda: defaultdict(int))
  most_bp = defaultdict(lambda: defaultdict(int))
  stats = {}

  for resultfn, result in results.items():
    cmp_methods = result['cmp_methods']
    consensus_methods = result['consensus_methods']
    guid = result['guid']

    stats[resultfn], precision_winner, recall_winner, most_bp_method = compare_methods(
      load_breakpoints(resultfn),
      guid,
      cmp_methods,
      consensus_methods,
      all_methods,
      cents_and_telos,
      indiv_bps[guid],
    )
    stats[resultfn]['consensus_methods'] = sorted(consensus_methods)
    stats[resultfn]['cmp_methods'] = sorted(cmp_methods)

    if precision_winner is None or recall_winner is None or most_bp_method is None:
      continue
    precision_winners[cmp_methods][precision_winner] += 1
    recall_winners[cmp_methods][recall_winner] += 1
    most_bp[cmp_methods][most_bp_method] += 1
    if precision_winner == recall_winner:
      total_dom_winners[cmp_methods][precision_winner] += 1

  return (stats, precision_winners, recall_winners, total_dom_winners, most_bp)

def generate_output(stats, precision_winners, recall_winners, total_dom_winners, most_bp):
  # Can't dump sets via JSON.
  for result in (precision_winners, recall_winners, total_dom_winners, most_bp):
    keys = list(result.keys())
    for K in keys:
      converted = ','.join(sorted(K))
      result[converted] = result[K]
      del result[K]

  return json.dumps({
    'precision_winners': precision_winners,
    'recall_winners': recall_winners,
    'total_domination_winners': total_dom_winners,
    'most_bp': most_bp,
    'stats': stats,
  })

def main():
  all_methods = sys.argv[1]
  centromeres = CentromereParser().load(sys.argv[2])
  runs = sys.argv[3:]

  cents_and_telos = parse_centromeres_and_telomeres(centromeres)
  all_methods = set(all_methods.split(','))
  datasets = defaultdict(list)

  results, indiv_bps = load_runs(runs, all_methods)
  stats, precision_winners, recall_winners, total_dom_winners, most_bp = determine_winners(results, all_methods, cents_and_telos, indiv_bps)
  print(generate_output(stats, precision_winners, recall_winners, total_dom_winners, most_bp))

main()
