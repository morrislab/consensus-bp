from __future__ import print_function
from collections import defaultdict
import json
import sys

def main():
  results = json.load(sys.stdin)
  anonymized = {
    'broad': 'C',
    'dkfz': 'F',
    'jabba': 'E',
    'mustonen095': 'D',
    'peifer': 'A',
    'vanloo_wedge_segs': 'B',
  }

  for K in ('most_bp', 'precision_winners', 'recall_winners', 'total_domination_winners'):
    wins = defaultdict(int)
    for comparison in results[K].values():
      for meth, methwins in comparison.items():
        wins[meth] += methwins

    ranked = sorted(wins.items(), key = lambda (M, W): W, reverse = True)

    # Anonymize and rank.
    ranked = [(anonymized[M], W, idx + 1) for idx, (M, W) in enumerate(ranked)]
    # Sort by name.
    ranked = sorted(ranked, key = lambda (M, W, R): -W)

    print(K)
    print(*['%s\t%s\t%s' % T for T in ranked], sep='\n')
    print('---')

main()
