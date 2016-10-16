#!/bin/bash
set -euo pipefail

export LC_ALL=C
module load gnu-parallel/20150822

SVDIR=~/work/exultant-pistachio/data/sv
BLACKLIST=~/work/exultant-pistachio/data/misc/blacklist.20160906.txt
CENTROMERES=~/work/exultant-pistachio/data/misc/cytoBand.txt.gz
NUM_NEEDED=$1
methods="broad dkfz mustonen095 peifer vanloo_wedge_segs jabba"

OUTDIR=~/work/exultant-pistachio/data/consensus_bp.thresh${NUM_NEEDED}

cd ~/work/exultant-pistachio/data/cnvs.pre_consensus_bp

mkdir -p $OUTDIR && rm -f $OUTDIR/*.{json,txt,stderr}
python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/run_comparison.py \
  --blacklist $BLACKLIST \
  --window-size 100000 \
  --centromere-filename $CENTROMERES \
  --optional-methods $(echo $methods | sed 's/ /,/g') \
  --num-needed-methods $NUM_NEEDED \
  $SVDIR \
  $OUTDIR \
  $methods \
  | parallel -j16

releaseid=$(date '+%Y%m%d').thresh${NUM_NEEDED}
cd $OUTDIR
for foo in *.json; do echo $(echo $foo | cut -d . -f1)$'\t'$(cat $foo | jq -r '.methods|sort|join(",")'); done > ~/work/exultant-pistachio/data/archives/consensus_bp_methods.$releaseid.txt
cd ~/work/exultant-pistachio/data
tar czf ~/work/exultant-pistachio/data/archives/consensus_bp.$releaseid.tar.gz consensus_bp/*.txt

cd ~/work/bp-witness
rm -f data/index.json
python2 index_data.py $OUTDIR
tar czf ~/work/exultant-pistachio/data/archives/bp_witness.$releaseid.tar.gz css index* js README.txt data/*.json --transform "s|^|bp-witness/|"
