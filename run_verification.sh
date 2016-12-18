#!/bin/bash
set -euo pipefail

#export LC_ALL=C
#module load gnu-parallel/20150822

SVDIR=~/work/exultant-pistachio/data/sv
BLACKLIST=~/work/exultant-pistachio/data/misc/blacklist.20160906.txt
CENTROMERES=~/work/exultant-pistachio/data/misc/cytoBand.txt.gz
SEX=~/work/exultant-pistachio/data/misc/inferred_sex.all_samples.20161209.txt
OUTDIR=~/work/exultant-pistachio/data/consensus_bp.verify.post_sv_sex
PLOTDIR=$OUTDIR/plots
methods="broad_x dkfz jabba mustonen095 vanloo_wedge_segs"

function create {
  mkdir -p $OUTDIR && rm -rf $OUTDIR/methods.*

  cd ~/work/exultant-pistachio/data/cnvs.pre_consensus_bp

  python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/run_verification.py \
    --blacklist $BLACKLIST \
    --window-size 100000 \
    --centromere-filename $CENTROMERES \
    --sex-filename $SEX \
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
