from __future__ import print_function, division
import sys
import os
import json
from make_consensus_breakpoints import CentromereParser, CHROM_LENS
from collections import defaultdict
import numpy as np

#import plotly
#import plotly.graph_objs as go

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

def include_near(breakpoints, around, threshold):
  # Note this method *will* include all breakpoints in `around`, which, for our
  # present purposes, is the desired behaviour.
  retained = defaultdict(list)

  for chrom in breakpoints.keys():
    for bp in breakpoints[chrom]:
      closest_dist = float('inf')
      closest_other = None
      for other in around[chrom]:
        dist = abs(bp['pos'] - other['pos'])
        if dist < closest_dist:
          # This is inefficient -- I should just break out of the parent loop
          # -- but whatever.
          closest_dist = dist
          closest_other = other
      if closest_dist <= threshold:
        retained[chrom].append(bp)

  return retained

def extract_all_sv(breakpoints):
  return extract(breakpoints, 'sv', True)

def extract_non_sv(breakpoints):
  return extract(breakpoints, 'sv', False)

def extract_bp_replaced_by_sv(breakpoints):
  return extract(breakpoints, 'sv_', True)

def extract_lone_svs(breakpoints):
  retained = defaultdict(list)
  for chrom in breakpoints.keys():
    for bp in breakpoints[chrom]:
      if bp['method'] == 'sv':
        retained[chrom].append(bp)
  return retained

def extract(breakpoints, prefix, truth):
  retained = defaultdict(list)
  for chrom in breakpoints.keys():
    for bp in breakpoints[chrom]:
      if bp['method'].startswith(prefix) is truth:
        retained[chrom].append(bp)
  return retained

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

def calc_relpos(chrom, pos):
  chrompos = pos / float(CHROM_LENS[chrom])
  chrom = unicode(chrom)
  if chrom.isnumeric():
    chrom = int(chrom)
  elif chrom == 'X':
    chrom = 23
  elif chrom == 'Y':
    chrom = 24
  else:
    raise Exception('penis')
  return chrom + chrompos

def find_closest(point, candidates):
  closest_point = None
  closest_dist = float('inf')
  for P in candidates:
    dist = abs(point['pos'] - P['pos'])
    if dist < closest_dist:
      closest_dist = dist
      closest_point = P
  return closest_point

def calc_dists_to_svs(bp, svs):
  dists = []
  for chrom in svs.keys():
    for sv in svs[chrom]:
      if chrom in bp.keys():
        closest = find_closest(sv, bp[chrom])
        dists.append(abs(sv['pos'] - closest['pos']))
      else:
        dists.append(-1)
  return dists

def calc_stats(fn, cents_and_telos):
  with open(fn) as F:
    dataset = os.path.basename(fn).split('.', 1)[0]
    J = json.load(F)
    bp = J['bp']['consensus']

    away_from_cents_and_telos = exclude_near(bp, cents_and_telos, 1e6)
    svs = extract_all_sv(away_from_cents_and_telos)
    away_from_sv = exclude_near(bp, svs, 1e5)
    near_sv_away_from_cents_and_telos = include_near(away_from_cents_and_telos, svs, 1e5)
    away_from_sv_and_cents_and_telos = exclude_near(away_from_cents_and_telos, svs, 1e5)
    assert count_bp(away_from_sv_and_cents_and_telos) + count_bp(near_sv_away_from_cents_and_telos) == count_bp(away_from_cents_and_telos)
    non_svs = extract_non_sv(away_from_cents_and_telos)
    assert count_bp(svs) + count_bp(non_svs) == count_bp(away_from_cents_and_telos)
    bp_replaced_by_sv = extract_bp_replaced_by_sv(away_from_cents_and_telos)
    lone_svs = extract_lone_svs(away_from_cents_and_telos)
    assert count_bp(bp_replaced_by_sv) + count_bp(lone_svs) == count_bp(svs)
    assert count_bp(bp_replaced_by_sv) <= count_bp(svs)

    dists_to_svs = {}
    for method in J['bp'].keys():
      if method == 'consensus':
        continue
      dists_to_svs[method] = calc_dists_to_svs(J['bp'][method], extract_all_sv(bp))

    stats = {
      'dataset': dataset,
      'bp_away_from_sv_and_cents_and_telos': away_from_sv_and_cents_and_telos,
      'dists_to_svs': dists_to_svs,
      'num_bp_away_from_sv_and_cents_and_telos': count_bp(away_from_sv_and_cents_and_telos),
      'num_bp_away_from_cents_and_telos': count_bp(away_from_cents_and_telos),
      'num_bp_away_from_sv': count_bp(away_from_sv),
      'num_total_bp': count_bp(bp),
      'num_bp_replaced_by_sv': count_bp(bp_replaced_by_sv),
      'num_lone_sv': count_bp(lone_svs),
      'num_non_sv': count_bp(non_svs),
      'num_bp_near_sv_away_from_cents_and_telos': count_bp(near_sv_away_from_cents_and_telos),
    }

    if count_bp(away_from_cents_and_telos) == 0:
      assert count_bp(away_from_sv_and_cents_and_telos) == 0
      stats['prop_bp_away_from_sv_and_cents_and_telos'] = 0
    else:
      stats['prop_bp_away_from_sv_and_cents_and_telos'] = count_bp(away_from_sv_and_cents_and_telos) / count_bp(away_from_cents_and_telos)

    if count_bp(svs) == 0:
      assert count_bp(bp_replaced_by_sv) == 0
      stats['prop_sv_with_proximal_bp'] = 0
    else:
      stats['prop_sv_with_proximal_bp'] = count_bp(bp_replaced_by_sv) / count_bp(svs)

    return stats

def cdf(arr):
  return (np.sort(arr), np.linspace(0, 1, len(arr), endpoint=False))

def scatter(traces, title, xtitle, ytitle, outfn, extra_shapes=None):
  if extra_shapes is None:
    extra_shapes = []

  layout = go.Layout(
    title = title,
    hovermode = 'closest',
    xaxis = {
      'title': xtitle,
    },
    yaxis = {
      'title': ytitle,
    },
    shapes = extra_shapes,
  )
  fig = go.Figure(data=traces, layout=layout)
  plotly.offline.plot(fig, filename=outfn)

def histogram(xvals, title, xtitle, ytitle, outfn):
  layout = go.Layout(
    title = title,
    xaxis = {
      'title': xtitle,
      #'range': [0, 1],
    },
    yaxis = {
      'title': ytitle,
    }
  )
  traces = [go.Histogram(x = xvals)]
  fig = go.Figure(data=traces, layout=layout)
  plotly.offline.plot(fig, filename=outfn)

def plot_distance_to_svs(dists_to_svs):
  scatter_traces = []
  for method in dists_to_svs.keys():
    cdf_xvals, cdf_yvals = cdf(dists_to_svs[method])
    scatter_traces.append(go.Scatter(
      x = cdf_xvals,
      y = cdf_yvals,
      mode = 'lines+markers',
      name = method,
    ))
  scatter(
    scatter_traces,
    'Distance from SVs to BPs',
    'Distance from SVs to nearest BP',
    'ECDF(x)',
    'sv_to_bp_dist.html'
  )

def plot_bp_away_props(prop_bp_away):
  datasets = sorted(prop_bp_away.keys())
  counts = [prop_bp_away[D][0] for D in datasets]
  props = [prop_bp_away[D][1] for D in datasets]
  cdf_xvals, cdf_yvals = cdf(props)

  scatter_traces = [go.Scatter(x = cdf_xvals, y = cdf_yvals, mode = 'lines+markers', name='CDF', text=datasets)]
  scatter(
    scatter_traces,
    'Proportion of consensus breakpoints unsupported by SVs', 
    'Proportion of BPs distal from cents & telos that are unsupported by SVs',
    'Tumors with this proportion or lower',
    'bp_unsupported_by_sv.ecdf.html'
  )

  histogram(
    props,
    'Proportion of consensus breakpoints unsupported by SVs', 
    'Proportion of BPs distal from cents & telos that are unsupported by SVs',
    'Number of tumors',
    'bp_unsupported_by_sv.hist.html'
  )

  scatter(
    [go.Scatter(x = counts, y = props, mode='markers', text=datasets)],
    'Breakpoints unsupported by SVs', 
    'Number of BPs distal from cents & telos that are unsupported by SVs',
    'Proportion of BPs distal from cents & telos that are unsupported by SVs',
    'bp_unsupported_by_sv.scatter.html'
  )

def plot_sv_away_from_bp_props(prop_sv_away):
  datasets = sorted(prop_sv_away.keys())
  props = [prop_sv_away[D] for D in datasets]
  histogram(
    props,
    'Proportion of consensus SVs without proximal BP (%s tumors)' % len(datasets),
    'Proportion of consensus SVs distal from cents & telos that are unsupported by BPs',
    'Number of tumors',
    'sv_unsupported_by_bp.hist.html'
  )

def plot_bp_away_positions(away_points, cents_and_telos):
  datasets, positions = zip(*away_points)
  xvals, yvals = cdf(positions)
  scatter_traces = [go.Scatter(x = xvals, y = yvals, mode = 'lines+markers', name='CDF', text=datasets)]

  shapes = []
  for chrom in cents_and_telos.keys():
    for pos in cents_and_telos[chrom]:
      relpos = calc_relpos(chrom, pos['pos'])
      shapes.append({
        'type': 'line',
        'x0': relpos,
        'y0': 0,
        'x1': relpos,
        'y1': 1,
        'line': {
          'width': 2,
          'dash': 'dot',
          'color': 'rgba(128, 128, 128, 0.5)',
        }
      })

  scatter(
    scatter_traces,
    'Positions of breakpoint distant from SVs',
    'Position of breakpoint distant from SVs',
    'Proportion of positions before location',
    'away_from_cents_and_telos.html',
    extra_shapes = shapes
  )



def main():
  centromeres = CentromereParser().load(sys.argv[1])
  cents_and_telos = parse_centromeres_and_telomeres(centromeres)
  bp_away_from_sv_and_cents_and_telos = []

  prop_bp_away = {}
  prop_sv_away = {}

  stat_types = (
    'dataset',
    'num_bp_away_from_sv_and_cents_and_telos',
    'prop_bp_away_from_sv_and_cents_and_telos',
    'num_bp_away_from_cents_and_telos',
    'num_bp_away_from_sv',
    'num_total_bp',
    'prop_sv_with_proximal_bp',
    'num_bp_replaced_by_sv',
    'num_lone_sv',
    'num_non_sv',
  )
  print(*stat_types, sep='\t')

  for fn in sys.argv[2:]:
    stats = calc_stats(fn, cents_and_telos)
    print(*[stats[T] for T in stat_types], sep='\t')
    continue

    dataset = stats['dataset']
    prop_bp_away[dataset] = (stats['num_bp_away_from_sv_and_cents_and_telos'], stats['prop_bp_away_from_sv_and_cents_and_telos'])
    prop_sv_away[dataset] = 1 - stats['prop_sv_with_proximal_bp']

    if prop_bp_away[dataset] >= 0.95:
      for chrom in stats['bp_away_from_sv_and_cents_and_telos']:
        for ap in stats['bp_away_from_sv_and_cents_and_telos'][chrom]:
          relpos = calc_relpos(chrom, ap['pos'])
          bp_away_from_sv_and_cents_and_telos.append((dataset, relpos))
  return

  plot_distance_to_svs(stats['dists_to_svs'])
  plot_bp_away_positions(bp_away_from_sv_and_cents_and_telos, cents_and_telos)
  plot_sv_away_from_bp_props(prop_sv_away)
  plot_bp_away_props(prop_bp_away)

main()
