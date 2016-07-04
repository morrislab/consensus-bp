import json
import sys
from collections import defaultdict

import plotly
import plotly.tools as pltools
import plotly.plotly as py
import plotly.graph_objs as go

def plot(scores, tumor_counts, title):
  bars = []
  N = len(scores.keys())
  titles = []

  for method in sorted(scores.keys()):
    method_bars = []
    summed = 0

    for group in sorted(scores[method].keys()):
      xvals = sorted(scores[method][group].keys())
      yvals = [scores[method][group][S] for S in xvals]
      summed += sum(yvals)

      method_bars.append(go.Bar(
        x = xvals,
        y = yvals,
        name = group,
      ))

    bars.append(method_bars)
    titles.append('%s (%s tumors, %e breakpoints)' % (method, tumor_counts[method], summed))

  print(titles)
  fig = pltools.make_subplots(
    rows=1,
    cols=N,
    subplot_titles=titles,
    #shared_yaxes=True,
  )
  for idx, method_bars in enumerate(bars):
    for group_bar in method_bars:
      fig.append_trace(group_bar, 1, idx + 1)

  fig['layout'].update(width = N*500, height=400, title=title, barmode='stack')
  for idx in range(N):
    #fig['layout']['yaxis%s' % (idx + 1)].update(type='log')
    pass

  plotly.offline.plot(fig, filename='%s.html' % title)

def combine_scores():
  # combined_scores[method][group][score]
  combined_scores = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
  tumor_counts = defaultdict(int)

  for line in sys.stdin:
    parsed = json.loads(line.strip())
    for method in parsed['scores'].keys():
      tumor_counts[method] += 1
      for group in parsed['scores'][method].keys():
        for score, count in parsed['scores'][method][group].items():
          combined_scores[method][group][int(score)] += count

  return (combined_scores, tumor_counts)

def main():
  scores, tumor_counts = combine_scores()
  title = sys.argv[1]
  plot(scores, tumor_counts, title)

if __name__ == '__main__':
  main()
