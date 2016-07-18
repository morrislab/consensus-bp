import json
import sys
from collections import defaultdict

import plotly
import plotly.tools as pltools
import plotly.plotly as py
import plotly.graph_objs as go

def plot_bp_mutual(scores, tumor_counts, title):
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
    titles.append('%s (%s tumors, %.2e breakpoints)' % (method, tumor_counts[method], summed))

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

  plotly.offline.plot(fig, filename='bp_mutual_%s.html' % title)

def plot_sv(bp_sv_scores, sv_bp_scores, title):
  assert sorted(bp_sv_scores.keys()) == sorted(sv_bp_scores.keys())
  N = len(bp_sv_scores.keys())

  bp_sv_bars = []
  sv_bp_bars = []
  titles = []
  num_bp = {}
  num_sv = {}

  for method in sorted(bp_sv_scores.keys()):
    svclasses = sorted(bp_sv_scores[method].keys())
    counts = [bp_sv_scores[method][svclass] for svclass in svclasses]
    bp_sv_bars.append(go.Bar(x = svclasses, y = counts, name = method))
    num_bp[method] = sum(counts)

  for method in sorted(sv_bp_scores.keys()):
    svclasses = sorted(sv_bp_scores[method].keys())
    true_counts = [sv_bp_scores[method][svclass]['True'] for svclass in svclasses]
    false_counts = [sv_bp_scores[method][svclass]['False'] for svclass in svclasses]
    sv_bp_bars.append([
      go.Bar(x = svclasses, y = true_counts, name = 'true'),
      go.Bar(x = svclasses, y = false_counts, name = 'false'),
    ])
    num_sv[method] = sum(true_counts) + sum(false_counts)

  for method in sorted(sv_bp_scores.keys()):
    titles.append('%s (%.2e BPs, %.2e SVs)' % (
      method,
      num_bp[method],
      num_sv[method]
    ))

  fig = pltools.make_subplots(
    rows=2,
    cols=N,
    subplot_titles=titles,
  )

  for idx, method_bars in enumerate(bp_sv_bars):
    fig.append_trace(method_bars, 1, idx + 1)
  for idx, method_bars in enumerate(sv_bp_bars):
    for boolean_bars in method_bars:
      fig.append_trace(boolean_bars, 2, idx + 1)

  fig['layout'].update(width = N*500, height=800, title=title, barmode='stack')
  plotly.offline.plot(fig, filename='sv_%s.html' % title)

def determine_well_supported_bp_prop(bpfn, ws_threshold):
  bp_vs_wsbp = defaultdict(dict)

  with open(bpfn) as bpf:
    for line in bpf:
      parsed = json.loads(line.strip())
      dataset = parsed['dataset']

      for method in parsed['bp_mutual_scores'].keys():
        method_bps = 0
        method_ws_bps = 0
        for group in parsed['bp_mutual_scores'][method].keys():
          method_bps += sum(parsed['bp_mutual_scores'][method][group].values())
          method_ws_bps += sum([
            V for (K, V) in parsed['bp_mutual_scores'][method][group].items() 
            if int(K) >= ws_threshold
          ])
        bp_vs_wsbp[method][dataset] = (method_ws_bps, method_bps)

  return bp_vs_wsbp

def plot_bp_vs_well_supported_bp(bp_vs_wsbp, title, log_axes=False):
  N = len(bp_vs_wsbp.keys())
  titles = []
  traces = []

  for method, points in bp_vs_wsbp.items():
    datasets = sorted(points.keys())
    wsbp, bp = zip(*[points[D] for D in datasets])
    traces.append(go.Scatter(
      x = bp,
      y = wsbp,
      mode = 'markers',
      text = datasets
    ))
    titles.append('%s (%s tumours)' % (method, len(datasets)))

  fig = pltools.make_subplots(
    rows=1,
    cols=N,
    subplot_titles=titles,
  )

  for idx in range(N):
    fig.append_trace(traces[idx], 1, idx + 1)
    fig['layout']['xaxis%s' % (idx + 1)].update(title='Total breakpoints', range=(0, 2000))
    fig['layout']['yaxis%s' % (idx + 1)].update(title='Well-supported breakpoints', range=(0, 600))
    if log_axes:
      for axis in ('xaxis', 'yaxis'):
        fig['layout']['%s%s' % (axis, idx + 1)].update(type='log')

  fig['layout'].update(width = N*500, height=400, title=title)
  plotly.offline.plot(fig, filename='wsbp_scatterplots_%s_%s.html' % (log_axes and 'log' or 'linear', title))

def combine_scores(bpfn):
  # bp_mutual_scores[method][group][score]
  bp_mutual_scores = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
  # bp_sv_scores[method][svclass]
  bp_sv_scores = defaultdict(lambda: defaultdict(int))
  # sv_bp_scores[method][svclass][is_proximal_bp]
  sv_bp_scores = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

  tumor_counts = defaultdict(int)

  with open(bpfn) as bpf:
    for line in bpf:
      parsed = json.loads(line.strip())

      for method in parsed['bp_mutual_scores'].keys():
        tumor_counts[method] += 1
        for group in parsed['bp_mutual_scores'][method].keys():
          for score, count in parsed['bp_mutual_scores'][method][group].items():
            bp_mutual_scores[method][group][int(score)] += count

      for method in parsed['bp_sv_scores'].keys():
        for svclass in parsed['bp_sv_scores'][method].keys():
          count = parsed['bp_sv_scores'][method][svclass]
          bp_sv_scores[method][svclass] += count

      for method in parsed['sv_bp_scores'].keys():
        for svclass in parsed['sv_bp_scores'][method].keys():
          for is_proximal_bp in parsed['sv_bp_scores'][method][svclass].keys():
            count = parsed['sv_bp_scores'][method][svclass][is_proximal_bp]
            sv_bp_scores[method][svclass][is_proximal_bp] += count

  # Print stats.
  for method in bp_sv_scores.keys():
    total = sum(bp_sv_scores[method].values())
    print('bp_sv_scores', method, bp_sv_scores[method]['null'] / float(total))
  for method in sv_bp_scores.keys():
    has_proximal = 0
    total = 0
    for svclass in sv_bp_scores[method].keys():
      has_proximal += sv_bp_scores[method][svclass]['True']
      total += sum(sv_bp_scores[method][svclass].values())
    print('sv_bp_scores', method, has_proximal / float(total))

  return (bp_mutual_scores, bp_sv_scores, sv_bp_scores, tumor_counts)

def main():
  bpfn = sys.argv[1]
  title = sys.argv[2]
  well_supported_threshold = 4

  bp_mutual_scores, bp_sv_scores, sv_bp_scores, tumor_counts = combine_scores(bpfn)
  plot_bp_mutual(bp_mutual_scores, tumor_counts, title)
  plot_sv(bp_sv_scores, sv_bp_scores, title)

  bp_vs_wsbp = determine_well_supported_bp_prop(bpfn, well_supported_threshold)
  for log_axes in (True, False):
    plot_bp_vs_well_supported_bp(bp_vs_wsbp, title, log_axes)

if __name__ == '__main__':
  main()
