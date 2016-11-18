from __future__ import print_function, division
import sys
import numpy as np
import csv

import plotly
import plotly.graph_objs as go

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

def load_data(datafn):
  vals = {}
  with open(datafn) as dataf:
    reader = csv.DictReader(dataf, delimiter='\t')
    for row in reader:
      vals[row['dataset']] = row
      del row['dataset']
  return vals

def main():
  traces = []

  for run in sys.argv[1:]:
    run_name, statsfn = run.split('=', 1)
    stats = load_data(statsfn)
    datasets = stats.keys()
    num_total_bp = np.array([int(stats[D]['num_total_bp']) for D in datasets])
    print(num_total_bp)

    xvals, yvals = cdf(num_total_bp)
    traces.append(go.Scatter(
      mode='lines+markers',
      x = xvals,
      y = yvals,
      text = datasets,
      name = '%s (%s samples)' % (run_name, len(num_total_bp)),
    ))

  scatter(traces, 'Number of total BPs', 'Number of BPs', 'ECDF(x)', 'num_bp.html', logx = True)

main()
