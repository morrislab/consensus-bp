#!/bin/bash
set -euo pipefail

export LC_ALL=C
module load gnu-parallel/20150822
OUTDIR=~/work/exultant-pistachio/data/tmp/breaks/
SVDIR=~/work/exultant-pistachio/data/sv
cd ~/work/exultant-pistachio/data/cnvs

cd $OUTDIR
for foo in *.json; do
  python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/plot_breakpoint_scores.py $foo $(echo $foo | cut -d . -f1-2)
done
exit

for W in 1000 5000 10000 20000 50000 100000 300000 1000000; do
  rm -f $OUTDIR/{un,}directed.$W.{json,stderr,html}

  python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/run_comparison.py \
    --window-size $W \
    --optional-methods vanloo_wedge,mustonen095,peifer,dkfz,broad \
    $SVDIR \
    broad dkfz mustonen095 peifer vanloo_wedge \
    | parallel -j16 \
    > $OUTDIR/directed.$W.json \
    2>$OUTDIR/directed.$W.stderr

  continue

  python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/run_comparison.py \
    --window-size $W \
    --optional-methods vanloo_wedge,mustonen095,peifer,dkfz,broad \
    --undirected \
    $SVDIR \
    broad dkfz mustonen095 peifer vanloo_wedge \
    | parallel -j16 \
    > $OUTDIR/undirected.$W.json \
    2>$OUTDIR/undirected.$W.stderr
done

