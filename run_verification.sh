#!/bin/bash
set -euo pipefail

#export LC_ALL=C
#module load gnu-parallel/20150822

SVDIR=~/work/exultant-pistachio/data/sv
BLACKLIST=~/work/exultant-pistachio/data/misc/blacklist.20160906.txt
CENTROMERES=~/work/exultant-pistachio/data/misc/cytoBand.txt.gz
OUTDIR=~/work/exultant-pistachio/data/consensus_bp.verify.post_sv
PLOTDIR=$OUTDIR/plots
methods="broad dkfz mustonen095 peifer vanloo_wedge"

function create {
  mkdir -p $OUTDIR && rm -rf $OUTDIR/methods.*

  cd ~/work/exultant-pistachio/data/cnvs.verify

  python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/run_verification.py \
    --blacklist $BLACKLIST \
    --window-size 100000 \
    --centromere-filename $CENTROMERES \
    --methods $(echo $methods | sed 's/ /,/g') \
    $SVDIR \
    $OUTDIR \
    $methods \
    | parallel -j40
}

function evaluate {
  mkdir -p $PLOTDIR && rm -f $PLOTDIR/stats.* $PLOTDIR/*.html
  cd $OUTDIR

  for run in methods.*; do
    runtype=$(echo $run | cut -d . -f2-)
    echo "python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/calc_breakpoint_stats.py" \
      "$CENTROMERES" \
      "$run"/'*-*.json' \
      "> $PLOTDIR/stats.$runtype.txt"
  done | parallel -j40
}

function main {
  create
  evaluate
}

main
