#!/bin/bash
set -euo pipefail

#export LC_ALL=C
#module load gnu-parallel/20150822

SVDIR=~/work/exultant-pistachio/data/sv
BLACKLIST=~/work/exultant-pistachio/data/misc/blacklist.20160906.txt
CENTROMERES=~/work/exultant-pistachio/data/misc/cytoBand.txt.gz
SEX=~/work/exultant-pistachio/data/misc/inferred_sex.all_samples.20161209.txt
OUTDIR=~/work/exultant-pistachio/data/consensus_bp.verify.post_sv
PLOTDIR=$OUTDIR/plots
methods="broad dkfz jabba mustonen095 peifer vanloo_wedge_segs"
PARALLEL=40

function create {
  mkdir -p $OUTDIR && rm -rf $OUTDIR/methods.* $OUTDIR/bp.*

  cd ~/work/exultant-pistachio/data/cnvs.pre_consensus_bp

  python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/run_verification.py \
    --blacklist $BLACKLIST \
    --window-size 100000 \
    --centromere-filename $CENTROMERES \
    --sex-filename $SEX \
    $SVDIR \
    $OUTDIR \
    $methods \
    | parallel -j$PARALLEL --halt 1

  cd $OUTDIR
  for foo in bp.*; do
    mv $foo $(echo $foo | sed 's/^bp/methods/')
  done
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
  done | parallel -j$PARALLEL --halt 1
}

function plot {
  cd $PLOTDIR
  python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/plot_verification.py "consensus_methods" stats.any*.txt
  python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/plot_verification.py "indiv_methods" $(ls stats*txt | grep -v "^stats.any")
  python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/plot_supported.py stats.any3_any2_conservative.txt
}

function main {
  create
  evaluate
  plot
}

main
