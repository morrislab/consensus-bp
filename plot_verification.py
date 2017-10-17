from __future__ import print_function, division
import sys
import numpy as np
import csv
import os
import glob

import plotly
import plotly.graph_objs as go

def parse_stats(statsfn):
  stats = {
    'tp': [],
    'fp': [],
    'fn': [],
    'num_input_sv': [],
    'num_bp_away_from_sv_and_cents_and_telos': [],
    'num_bp_away_from_cents_and_telos': [],
    'numbp': [],
    'num_sv_away_from_cents_and_telos': [],
    'datasets': [],
  }
  labels = {}

  with open(statsfn) as F:
    reader = csv.DictReader(F, delimiter='\t')
    for row in reader:
      stats['tp'].append(float(row['num_bp_replaced_by_sv']))
      stats['fp'].append(float(row['num_non_sv']))
      stats['fn'].append(float(row['num_lone_sv']))
      stats['num_sv_away_from_cents_and_telos'].append(float(row['num_bp_replaced_by_sv']) + float(row['num_lone_sv']))
      stats['datasets'].append(row['dataset'])

  for K in ('tp', 'fp', 'fn', 'num_input_sv', 'num_bp_away_from_sv_and_cents_and_telos', 'num_bp_away_from_cents_and_telos', 'num_sv_away_from_cents_and_telos'):
    stats[K] = np.array(stats[K], dtype=np.float)

  stats['numbp'] = stats['tp'] + stats['fp']
  labels['numbp'] = stats['datasets']

  stats['precision'] = stats['tp'] / (stats['tp'] + stats['fp'])
  stats['recall'] = stats['tp'] / (stats['tp'] + stats['fn'])
  labels['precision'] = labels['recall'] = stats['datasets']
  oldlen = len(stats['precision'])
  assert oldlen == len(stats['recall'])

  notnan_idxs = np.logical_not(np.logical_or(np.isnan(stats['precision']), np.isnan(stats['recall'])))
  stats['precision'] = stats['precision'][notnan_idxs]
  stats['recall'] = stats['recall'][notnan_idxs]
  assert len(stats['precision']) == len(stats['recall'])
  #print(statsfn, 'has', oldlen - len(stats['precision']), 'nan')

  stats['nonsv_ratio'] = (stats['fp'] + 1) / (stats['num_sv_away_from_cents_and_telos'] + 1)
  assert np.count_nonzero(np.isnan(stats['nonsv_ratio'])) == 0
  assert np.count_nonzero(stats['nonsv_ratio'] == 0) == 0
  stats['nonsv_ratio'] = np.log2(stats['nonsv_ratio'])
  labels['nonsv_ratio'] = ['%s (nonsv = %s, sv = %s)' % (stats['datasets'][idx], stats['fp'][idx], stats['num_sv_away_from_cents_and_telos'][idx]) for idx in range(len(stats['datasets']))]

  return (stats, labels)

def scatter(traces, title, xtitle, ytitle, outfn, logx = False, xmin = None, xmax = None,):
  xaxis = {
    'title': xtitle,
    'type': logx and 'log' or 'linear',
    'range': [xmin, xmax],
  }

  layout = go.Layout(
    title = title,
    hovermode = 'closest',
    xaxis = xaxis,
    yaxis = {
      'title': ytitle,
    },
  )
  fig = go.Figure(data=traces, layout=layout)
  plotly.offline.plot(fig, filename=outfn)

def cdf(arr, labels=None):
  sorted_idxs = np.argsort(arr)
  ret = [
    arr[sorted_idxs],
    np.linspace(0, 1, len(arr), endpoint=False),
  ]
  if labels is not None:
    ret.append([labels[idx] for idx in sorted_idxs])
  return tuple(ret)

def plot_method_combos(statsfns):
  xvals, yvals, xerrors, yerrors = [], [], [], []
  runs = []

  for statsfn in statsfns:
    run = os.path.basename(statsfn).split('.')[1]
    runs.append(run)
    stats, _ = parse_stats(statsfn)

    xvals.append(np.mean(stats['recall']))
    yvals.append(np.mean(stats['precision']))
    xerrors.append(np.std(stats['recall']))
    yerrors.append(np.std(stats['precision']))

  show_error_bars = False
  method_combo_trace = go.Scatter(
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
    [method_combo_trace],
    'Performance of different combinations of (broad, dkfz, jabba, mustonen095, peifer, vanloo_wedge)',
    'Recall',
    'Precision',
    'method_combos_perf.html',
  )

def plot_ecdfs(run_label, statsfns):
  traces = {
    'precision': [],
    'recall': [],
    'numbp': [],
    'nonsv_ratio': [],
  }

  for statsfn in statsfns:
    run = os.path.basename(statsfn).split('.')[1]
    print(statsfn, run)
    stats, labels = parse_stats(statsfn)

    for plot in traces.keys():
      X, Y, L = cdf(stats[plot], labels[plot])

      if run_label == 'consensus_methods':
        line = {'width': 4}
      else:
        line = {'dash': 'dot', 'width': 4}

      traces[plot].append(go.Scatter(
        mode='lines',
        x = X,
        y = Y,
        text = L,
        name = '%s (%s values)' % (run, len(L)),
        line = line,
        # Comma corresponds to "Y_dkfz,jabba", which uses both (along with SVs).
        #visible = (('any' in run or ',' in run) and True or 'legendonly'),
      ))

  #for T in traces['numbp']:
  #  if T['name'].startswith('any'):
  #    T['visible'] = False

  scatter(
    traces['precision'],
    'Precision ECDF',
    'Precision',
    'ECDF(x)',
    'precision_ecdf.%s.html' % run_label,
  )
  scatter(
    traces['recall'],
    'Recall ECDF',
    'Recall',
    'ECDF(x)',
    'recall_ecdf.%s.html' % run_label,
  )
  scatter(
    traces['numbp'],
    '# BPs ECDF',
    '# BPs',
    'ECDF(x)',
    'numbp_ecdf.%s.html' % run_label,
    logx = True,
    #xmin = 1.9,
    #xmax = 4,
  )
  scatter(
    traces['nonsv_ratio'],
    'Non-SV ratio ECDF',
    'log2((# non-SV BPs + 1) / (# SVs + 1))',
    'ECDF(x)',
    'nonsv_ratio_ecdf.%s.html' % run_label,
    xmin = -10,
    xmax = 10,
  )

def main():
  run_label = sys.argv[1]
  assert run_label in ('consensus_methods', 'indiv_methods')
  plot_method_combos(sys.argv[2:])
  plot_ecdfs(run_label, sys.argv[2:])

main()
