#!/bin/bash
set -euo pipefail

export LC_ALL=C
module load gnu-parallel/20150822
OUTDIR=~/work/exultant-pistachio/data/consensus_bp
SVDIR=~/work/exultant-pistachio/data/sv
BLACKLIST=~/work/exultant-pistachio/data/misc/blacklist.20160815.txt
CENTROMERES=~/work/exultant-pistachio/data/misc/cytoBand.txt.gz

cd ~/work/exultant-pistachio/data/cnvs

rm -f $OUTDIR/*.{json,txt,stderr}
python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/run_comparison.py \
  --blacklist $BLACKLIST \
  --window-size 100000 \
  --centromere-filename $CENTROMERES \
  --optional-methods vanloo_wedge_segs,mustonen095,peifer,dkfz,broad \
  $SVDIR \
  $OUTDIR \
  broad dkfz mustonen095 peifer vanloo_wedge_segs \
  | parallel -j16
