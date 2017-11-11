import plotly
import plotly.tools as pltools
import plotly.plotly as py
import plotly.graph_objs as go

import csv
import sys
import numpy as np
import scipy.stats

def scatter(title, xtitle, ytitle, X, Y, L, outfn, logx = False, xmin = None, xmax = None,):
  print(title, len(X), np.mean(Y), np.std(Y), np.percentile(Y, [25, 50, 75]))
  if logx:
    XY = np.vstack((np.log(X), Y))
  else:
    XY = np.vstack((X, Y))
  density = scipy.stats.gaussian_kde(XY)(XY)
  assert len(X) == len(Y) == len(density)

  trace = go.Scatter(
    x = X,
    y = Y,
    mode = 'markers',
    text = L,
    marker = dict(
      color = density,
      colorscale = 'Viridis',
    ),
  )

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
  fig = go.Figure(data=[trace], layout=layout)
  plotly.offline.plot(fig, filename=outfn)

def extract_vals(rows, key):
  return [R[key] for R in rows]

def extract_floats(rows, key):
  return np.array([float(V) for V in extract_vals(rows, key)])

def pie(title, values, labels, outfn):
  trace = go.Pie(labels=labels, values=values)
  layout = go.Layout(title = title)
  fig = go.Figure(data=[trace], layout=layout)
  plotly.offline.plot(fig, filename=outfn)

def plot(cnvfn):
  with open(cnvfn) as cnvf:
    reader = csv.DictReader(cnvf, delimiter='\t')
    rows = list(reader)
  sampids = extract_vals(rows, 'dataset')
  vals = {K: np.array(extract_floats(rows, K)) for K in (
    'num_bp_replaced_by_sv',
    'num_lone_sv',
    'num_bp_replaced_by_sv',
    'num_bp_away_from_cents_and_telos'
  )}

  num_sv = vals['num_lone_sv'] + vals['num_bp_replaced_by_sv']
  frac_sv_supported_by_cbps = vals['num_bp_replaced_by_sv'] / (vals['num_lone_sv'] + vals['num_bp_replaced_by_sv'])
  sv_indices = num_sv > 0
  print('Discarded %s vals from sv_indices' % (len(sampids) - np.sum(sv_indices)))

  num_bp = vals['num_bp_away_from_cents_and_telos']
  frac_cbp_supported_by_svs = vals['num_bp_replaced_by_sv'] / vals['num_bp_away_from_cents_and_telos']
  bp_indices = num_bp > 0
  print('Discarded %s vals from bp_indices' % (len(sampids) - np.sum(bp_indices)))

  scatter(
    title = 'Support of SVs by consensus BPs',
    xtitle = 'Number of SVs',
    ytitle = 'Fraction of SVs with associated consensus BP',
    X = num_sv[sv_indices],
    Y = frac_sv_supported_by_cbps[sv_indices],
    L = [S for S, use_S in zip(sampids, sv_indices) if use_S is True],
    outfn = 'scatter.sv.html',
    logx = True,
  )
  scatter(
    title = 'Support of consensus BPs by SVs',
    xtitle = 'Number of BPs',
    ytitle = 'Fraction of BPs with associated consensus SV',
    X = num_bp[bp_indices],
    Y = frac_cbp_supported_by_svs[bp_indices],
    L = [S for S, use_S in zip(sampids, bp_indices) if use_S is True],
    outfn = 'scatter.bp.html',
    logx = True,
  )

  pie(
    title = 'Support of consensus SVs by BPs',
    values = [np.sum(vals['num_bp_replaced_by_sv']), np.sum(vals['num_lone_sv'])],
    labels = ['SVs supported by consensus BP', 'SVs not supported by consensus BP'],
    outfn = 'pie.sv.html',
  )

  pie(
    title = 'Support of consensus BPs by SVs',
    values = [np.sum(vals['num_bp_replaced_by_sv']), np.sum(vals['num_bp_away_from_cents_and_telos']) - np.sum(vals['num_bp_replaced_by_sv'])],
    labels = ['Consensus BPs supported by SV', 'Consensus BPs not supported by SV'],
    outfn = 'pie.bp.html',
  )

def main():
  cnvfn = sys.argv[1]
  plot(cnvfn)

main()
