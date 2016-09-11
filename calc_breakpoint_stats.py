from __future__ import print_function, division
import sys
import os
import json
from make_consensus_breakpoints import CentromereParser, CHROM_LENS
from collections import defaultdict

def exclude_near(breakpoints, around, threshold):
  retained = defaultdict(list)

  for chrom in breakpoints.keys():
    for bp in breakpoints[chrom]:
      closest_dist = float('inf')
      closest_other = None
      for other in around[chrom]:
        dist = abs(bp['pos'] - other['pos'])
        if dist < closest_dist:
          closest_dist = dist
          closest_other = other
      if closest_dist > threshold:
        retained[chrom].append(bp)

  return retained

def extract_all_sv(breakpoints):
  return extract(breakpoints, 'sv', True)

def extract_non_sv(breakpoints):
  return extract(breakpoints, 'sv', False)

def extract_replaced_by_sv(breakpoints):
  return extract(breakpoints, 'sv_', True)

def extract(breakpoints, prefix, truth):
  retaine = defaultdict(list)
  for chrom in breakpoints.keys():
    for bp in breakpoints[chrom]:
      if bp['method'].startswith(prefix) is truth:
        retaine[chrom].append(bp)
  return retaine

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
    for P in points.values():
      cents_and_telos[chrom].append({'pos': P})
    cents_and_telos[chrom].sort()

  return cents_and_telos

def count_bp(bp):
  return sum([len(V) for V in bp.values()])

def process(fn, cents_and_telos):
  with open(fn) as F:
    dataset = os.path.basename(fn)
    J = json.load(F)
    bp = J['bp']['consensus']

    away_from_cents_and_telos = exclude_near(bp, cents_and_telos, 1e6)
    svs = extract_all_sv(away_from_cents_and_telos)
    away_from_sv_and_cents_and_telos = exclude_near(away_from_cents_and_telos, svs, 1e5)
    non_svs = extract_non_sv(away_from_cents_and_telos)
    assert count_bp(svs) + count_bp(non_svs) == count_bp(away_from_cents_and_telos)
    replaced_by_sv = extract_replaced_by_sv(away_from_cents_and_telos)
    assert count_bp(replaced_by_sv) <= count_bp(svs)

    print(
      dataset,

      count_bp(away_from_sv_and_cents_and_telos),
      count_bp(away_from_cents_and_telos),
      count_bp(away_from_sv_and_cents_and_telos) / count_bp(away_from_cents_and_telos),

      count_bp(non_svs),
      count_bp(non_svs) + count_bp(svs),
      count_bp(non_svs) / (count_bp(non_svs) + count_bp(svs)),

      count_bp(replaced_by_sv),
      count_bp(bp),

      J['stats']['before_removing_proximal'] - J['stats']['after_removing_proximal'],
      (J['stats']['before_removing_proximal'] - J['stats']['after_removing_proximal']) / J['stats']['before_removing_proximal'],
      count_bp(extract(away_from_cents_and_telos, 'added_from_intersection', True)) / count_bp(away_from_cents_and_telos),

      sep='\t'
    )

def main():
  centromeres = CentromereParser().load(sys.argv[1])
  cents_and_telos = parse_centromeres_and_telomeres(centromeres)

  print(
    'dataset',

    'away_from_sv_and_cents_and_telos',
    'away_from_cents_and_telos',
    'proportion_away_from_sv_and_cents_and_telos',

    'non_sv',
    'blah',
    'blah2',

    'replaced_by_sv',
    'total_bp',

    'proximal_removed',
    'proportion_proximal_removed',
    'proportion_added_from_intersection',

    sep='\t'
  )
  for fn in sys.argv[2:]:
    process(fn, cents_and_telos)

main()
