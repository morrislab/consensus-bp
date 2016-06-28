import json
import sys
from collections import defaultdict

import plotly
import plotly.tools as pltools
import plotly.plotly as py
import plotly.graph_objs as go

def plot(scores, method_counts, title):
  methods = sorted(scores.keys())
  bars = []
  N = len(methods)
  titles = []

  for method in methods:
    xvals = sorted(scores[method].keys())
    yvals = [scores[method][S] for S in xvals]

    bars.append(go.Bar(
      x = xvals,
      y = yvals,
    ))
    titles.append('%s (%s tumors, %e breakpoints)' % (method, method_counts[method], sum(yvals)))

  print(titles)
  fig = pltools.make_subplots(
    rows=1,
    cols=N,
    subplot_titles=titles,
    #shared_yaxes=True,
  )
  for idx, bar in enumerate(bars):
    fig.append_trace(bar, 1, idx + 1)

  fig['layout'].update(width = N*400, height=400, title=title)
  for idx in range(N):
    fig['layout']['yaxis%s' % (idx + 1)].update(type='log')

  plotly.offline.plot(fig, filename='%s.html' % title)

def combine_scores():
  combined_scores = defaultdict(lambda: defaultdict(int))
  method_counts = defaultdict(int)

  for line in sys.stdin:
    parsed = json.loads(line.strip())
    for method, scores in parsed['scores'].items():
      method_counts[method] += 1
      for score, count in scores.items():
        combined_scores[method][int(score)] += count

  return (combined_scores, method_counts)

def main():
  scores, method_counts = combine_scores()
  title = sys.argv[1]
  plot(scores, method_counts, title)

if __name__ == '__main__':
  main()
