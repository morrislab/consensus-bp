from __future__ import print_function, division
import sys
import json
import numpy as np
import csv
import os

import plotly
import plotly.graph_objs as go

def parse_stats(statsfn):
  stats = {
    'tp': [],
    'fp': [],
    'fn': [],
    'datasets': [],
  }

  with open(statsfn) as F:
    reader = csv.DictReader(F, delimiter='\t')
    for row in reader:
      stats['tp'].append(float(row['num_bp_replaced_by_sv']))
      stats['fp'].append(float(row['num_non_sv']))
      stats['fn'].append(float(row['num_lone_sv']))
      stats['datasets'].append(row['dataset'])

  for K in ('tp', 'fp', 'fn'):
    stats[K] = np.array(stats[K], dtype=np.float)

  stats['precision'] = stats['tp'] / (stats['tp'] + stats['fp'])
  stats['recall'] = stats['tp'] / (stats['tp'] + stats['fn'])
  oldlen = len(stats['precision'])
  assert oldlen == len(stats['recall'])

  notnan_idxs = np.logical_not(np.logical_or(np.isnan(stats['precision']), np.isnan(stats['recall'])))
  stats['precision'] = stats['precision'][notnan_idxs]
  stats['recall'] = stats['recall'][notnan_idxs]
  assert len(stats['precision']) == len(stats['recall'])
  print(statsfn, 'has', oldlen - len(stats['precision']), 'nan')

  #stats['precision'] = np.nan_to_num(stats['precision'])
  #stats['recall'] = np.nan_to_num(stats['recall'])

  return stats

def scatter(traces, title, xtitle, ytitle, outfn, logx = False):
  layout = go.Layout(
    title = title,
    hovermode = 'closest',
    xaxis = {
      'title': xtitle,
      'type': logx and 'log' or 'linear',
    },
    yaxis = {
      'title': ytitle,
    },
  )
  fig = go.Figure(data=traces, layout=layout)
  plotly.offline.plot(fig, filename=outfn)

def cdf(arr):
  return (np.sort(arr), np.linspace(0, 1, len(arr), endpoint=False))

def main():
  xvals, yvals, xerrors, yerrors = [], [], [], []
  runs = []

  precision_traces = []
  recall_traces = []
  numbp_traces = []

  for statsfn in sys.argv[1:]:
    run = os.path.basename(statsfn).split('.')[1]
    runs.append(run)
    stats = parse_stats(statsfn)
    xvals.append(np.mean(stats['recall']))
    yvals.append(np.mean(stats['precision']))
    xerrors.append(np.std(stats['recall']))
    yerrors.append(np.std(stats['precision']))

    for collection, vals in ((recall_traces, stats['recall']), (precision_traces, stats['precision']), (numbp_traces, stats['tp'] + stats['fp'])):
      X, Y = cdf(vals)
      collection.append(go.Scatter(
        mode='lines+markers',
        x = X,
        y = Y,
        text = stats['datasets'],
        name = run,
      ))

  show_error_bars = False
  trace = go.Scatter(
    mode = 'markers',
    x = xvals,
    y = yvals,
    text = runs,
    error_x = {
      'type': 'data',
      'array': xerrors,
      'visible': show_error_bars,
    },
    error_y = {
      'type': 'data',
      'array': yerrors,
      'visible': show_error_bars,
    },
  )

  scatter(
    [trace],
    'Performance of different combinations of (broad, dkfz, mustonen095, peifer, vanloo_wedge)',
    'Recall',
    'Precision',
    'method_combos_perf.html',
  )
  scatter(
    precision_traces,
    'Precision ECDF',
    'Precision',
    'ECDF(x)',
    'precision_ecdf.html',
  )
  scatter(
    recall_traces,
    'Recall ECDF',
    'Recall',
    'ECDF(x)',
    'recall_ecdf.html',
  )
  scatter(
    numbp_traces,
    '# BPs ECDF',
    '# BPs',
    'ECDF(x)',
    'numbp_ecdf.html',
    logx = True,
  )

main()
