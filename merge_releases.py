from __future__ import print_function
import sys
import os
import glob
import json
from collections import defaultdict, OrderedDict
from run_comparison import parse_sex

# Counting methods for a given run:
#   (D=6; head -n1 list.txt | cut -f$D; cat list.txt | awk '$2=="male"{print $0}' | cut -f$D | python2 -c 'import sys; print("\n".join([str(len(L.split(","))) for L in sys.stdin]))' | sort | uniq -c | sort -nr | awk -v OFS='\t' '{print $1,$2}')

def merge(txt_inputs, times_sample_used, run_names, outdir):
  for guid in sorted(txt_inputs.keys()):
    header = None
    outpath = os.path.join(outdir, '%s.txt' % guid)
    assert not os.path.exists(outpath)

    with open(outpath, 'w') as outf:
      for R in run_names:
        if R not in txt_inputs[guid]:
          continue
        bpfn = txt_inputs[guid][R]
        times_sample_used[guid] += 1

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

  for guid in json_inputs.keys():
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

def main():
  sexfn = sys.argv[1]
  outdir = sys.argv[2]
  indirs = sys.argv[3:]

  sex = parse_sex(sexfn)

  times_sample_used = defaultdict(int)
  txt_inputs = load_sample_list(indirs, 'txt')
  json_inputs = load_sample_list(indirs, 'json')

  os.makedirs(outdir)
  merge(txt_inputs, times_sample_used, indirs, outdir)
  assert set(times_sample_used.values()) == set([2, 3])
  list_methods(json_inputs, sex, indirs)

if __name__ == '__main__':
  main()
