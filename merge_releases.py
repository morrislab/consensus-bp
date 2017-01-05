from __future__ import print_function
import sys
import os
import glob
import json
from collections import defaultdict, OrderedDict
from run_comparison import parse_sex
import numpy as np

# Counting methods for a given run:
#   (D=6; head -n1 list.txt | cut -f$D; cat list.txt | awk '$2=="male"{print $0}' | cut -f$D | python2 -c 'import sys; print("\n".join([str(len(L.split(","))) for L in sys.stdin]))' | sort | uniq -c | sort -nr | awk -v OFS='\t' '{print $1,$2}')

def merge(txt_inputs, inputs_used, run_names, outdir):
  for sampidx, guid in enumerate(sorted(txt_inputs.keys())):
    header = None
    outpath = os.path.join(outdir, '%s.txt' % guid)
    assert not os.path.exists(outpath)

    with open(outpath, 'w') as outf:
      for ridx, R in enumerate(run_names):
        if R not in txt_inputs[guid]:
          continue
        bpfn = txt_inputs[guid][R]
        inputs_used[sampidx, ridx] += 1

        with open(bpfn) as bpf:
          fh = next(bpf).strip()
          if header is None:
            header = fh
            print(header, file=outf)
          else:
            assert header == fh
          for L in bpf:
            print(L.strip(), file=outf)

def load_sample_list(indirs, ext):
  inputs = defaultdict(dict)

  for I in indirs:
    for D in glob.glob('%s/*-*.%s' % (I, ext)):
      guid = os.path.basename(D).split('.')[0]
      assert I not in inputs[guid]
      inputs[guid][I] = D

  return inputs

def list_methods(json_inputs, sex, run_names):
  header = '\t'.join(('guid', 'sex') + tuple(run_names))
  print(header)

  for guid in sorted(json_inputs.keys()):
    methods = []
    for R in run_names:
      if R in json_inputs[guid]:
        with open(json_inputs[guid][R]) as JF:
          J = json.load(JF)
        methods.append(J['consensus_methods'])
      else:
        methods.append([])
    methods = [','.join(sorted(M)) for M in methods]
    methods = [str(len(M) > 0 and M or None) for M in methods]
    line = '\t'.join((guid, sex[guid]) + tuple(methods))
    print(line)

def verify(inputs_used, run_names):
  assert np.all(0 <= inputs_used) and np.all(inputs_used <= 1)
  times_sampled_used = np.sum(inputs_used, axis=1)
  assert np.all(2 <= times_sampled_used) and np.all(times_sampled_used <= 3)

  assert set(run_names) == set(('bp.autosomes', 'bp.X_female', 'bp.X_male', 'bp.Y'))
  indices = {R: I for (I, R) in enumerate(run_names)}
  assert np.all(inputs_used[:,indices['bp.autosomes']] == 1)
  assert np.all(inputs_used[:,indices['bp.X_male']] == inputs_used[:,indices['bp.Y']])
  assert np.all(inputs_used[:,indices['bp.X_female']] + inputs_used[:,indices['bp.X_male']] == 1)

def main():
  sexfn = sys.argv[1]
  outdir = sys.argv[2]
  indirs = sys.argv[3:]

  sex = parse_sex(sexfn)

  txt_inputs = load_sample_list(indirs, 'txt')
  json_inputs = load_sample_list(indirs, 'json')
  assert set(txt_inputs.keys()) == set(json_inputs.keys())

  inputs_used = np.zeros((len(txt_inputs), len(indirs)))

  os.makedirs(outdir)
  merge(txt_inputs, inputs_used, indirs, outdir)
  verify(inputs_used, indirs)

  list_methods(json_inputs, sex, indirs)

if __name__ == '__main__':
  main()
