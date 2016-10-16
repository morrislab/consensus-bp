import plotly
import plotly.graph_objs as go
import sys
import json
import csv

def scatter(xvals, yvals, title, xtitle, ytitle, labels, outfn, draw_diagonal=True):
  if draw_diagonal:
    shapes = [{
        'type': 'line',
        'xref': 'x',
        'yref': 'y',
        'x0': 0,
        'y0': 0,
        'x1': min(max(xvals), max(yvals)),
        'y1': min(max(xvals), max(yvals)),
        'line': {
          'width': 3,
        },
      }]
  else:
    shapes = []

  layout = go.Layout(
    title = title,
    hovermode = 'closest',
    xaxis = {
      'title': xtitle,
      #'range': [0, 1],
    },
    yaxis = {
      'title': ytitle,
    },
    shapes = shapes
  )
  data = [go.Scatter(x = xvals, y = yvals, mode = 'markers', text = labels)]
  fig = go.Figure(data=data, layout=layout)

  plotly.offline.plot(fig, filename=outfn)

def histogram(vals, title, xtitle, ytitle, outfn):
  layout = go.Layout(
    title = title,
    xaxis = {
      'title': xtitle,
      #'range': [0, 1],
    },
    yaxis = {
      'title': ytitle,
    },
  )
  data = [go.Histogram(x = vals)]
  fig = go.Figure(data=data, layout=layout)

  plotly.offline.plot(fig, filename=outfn)


def load_data(datafn):
  vals = {}
  with open(datafn) as dataf:
    reader = csv.DictReader(dataf, delimiter='\t')
    for row in reader:
      vals[row['dataset']] = row
      del row['dataset']
  return vals

def main():
  xrun, yrun = sys.argv[1], sys.argv[2]
  oldfn, newfn = sys.argv[3], sys.argv[4]
  old_vals, new_vals = load_data(oldfn), load_data(newfn)

  common = sorted(set(old_vals.keys()) & set(new_vals.keys()))
  new_tumors = sorted(new_vals.keys())

  for metric in ('prop_bp_away_from_sv_and_cents_and_telos', 'prop_sv_with_proximal_bp', 'num_total_bp'):
    xvals = [float(old_vals[D][metric]) for D in common]
    yvals = [float(new_vals[D][metric]) for D in common]
    assert len(xvals) == len(yvals)
    scatter(
      xvals,
      yvals,
      'Comparison of %s & %s on %s (%s samples)' % (xrun, yrun, metric, len(xvals)),
      '%s in %s' % (xrun, metric),
      '%s in %s' % (yrun, metric),
      common,
      'old_vs_new.%s.html' % metric,
      draw_diagonal=True
    )
  return

  scatter(
    [float(new_vals[D]['away_from_sv_and_cents_and_telos']) for D in common],
    [float(new_vals[D]['proportion_away_from_sv_and_cents_and_telos']) for D in common],
    'Proportion vs. count comparison for away_from_sv_and_cents_and_telos',
    'away_from_sv_and_cents_and_telos',
    'proportion_away_from_sv_and_cents_and_telos',
    common,
    'proportion_vs_count.html',
    draw_diagonal=False
  )

  histogram(
    [float(new_vals[D]['proportion_away_from_sv_and_cents_and_telos']) for D in common],
    'Proportion of breakpoints >100 kb from SVs, centromeres, and telomeres',
    'Proportion of breakpoints >100 kb from SVs, centromeres, and telomeres',
    'Tumors',
    'proportion_away_from_sv_and_cents_and_telos.html',
  )

  histogram(
    [float(new_vals[D]['proportion_proximal_removed']) for D in new_tumors],
    'Proportion of proximal non-SV breakpoints removed',
    'Proportion of proximal non-SV breakpoints removed',
    'Tumors',
    'proportion_proximal_removed.html',
  )

  histogram(
    [float(new_vals[D]['proportion_added_from_intersection']) for D in new_tumors],
    'Proportion of consensus BPs added from intersection without support',
    'Proportion of consensus BPs added from intersection without support',
    'Tumors',
    'proportion_added_from_intersection.html',
  )

  scatter(
    [float(new_vals[D]['proximal_removed']) for D in new_tumors],
    [float(new_vals[D]['proportion_proximal_removed']) for D in new_tumors],
    'Proximal non-SV breakpoints removed',
    'Proximal non-SV breakpoints removed',
    'Proportion of proximal non-SV breakpoints removed',
    new_tumors,
    'proximal_removed.html',
    draw_diagonal=False
  )

main()
