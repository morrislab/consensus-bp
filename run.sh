#!/bin/bash
set -euo pipefail

#export LC_ALL=C
#module load gnu-parallel/20150822

SVDIR=~/work/exultant-pistachio/data/sv
BLACKLIST=~/work/exultant-pistachio/data/misc/blacklist.20160906.txt
CENTROMERES=~/work/exultant-pistachio/data/misc/cytoBand.txt.gz
methods="broad dkfz jabba mustonen095 peifer vanloo_wedge_segs"

NUM_NEEDED=3
suffix=thresh_${NUM_NEEDED}
OUTDIR=~/work/exultant-pistachio/data/consensus_bp.$suffix
releaseid=$(date '+%Y%m%d').$suffix
sample_list=~/work/exultant-pistachio/data/archives/consensus_bp_methods.$releaseid.txt

function make_sample_list {
  for foo in *.json; do
    echo $(echo $foo | cut -d . -f1)$'\t'$(cat $foo | jq -r '.methods|sort|join(",")')
  done > $sample_list
}

function main {
  cd ~/work/exultant-pistachio/data/cnvs.pre_consensus_bp

  mkdir -p $OUTDIR && rm -f $OUTDIR/*.{json,txt,stderr}
  python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/run_comparison.py \
    --blacklist $BLACKLIST \
    --window-size 100000 \
    --centromere-filename $CENTROMERES \
    --num-needed-methods $NUM_NEEDED \
    $SVDIR \
    $OUTDIR \
    $methods \
    | parallel -j40

  cd $OUTDIR
  rm -rf too_few && mkdir too_few
  make_sample_list
  for foo in $(cat $sample_list | python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/filter.py | cut -d ' ' -f 1); do
    mv $foo* too_few
  done
  make_sample_list

  cd ~/work/exultant-pistachio/data
  tar czf ~/work/exultant-pistachio/data/archives/consensus_bp.$releaseid.tar.gz consensus_bp.$suffix/*.txt
  tar czf ~/work/exultant-pistachio/data/archives/bp_witness.$releaseid.tar.gz consensus_bp.$suffix/*.json

  exit

  cd ~/work/bp-witness
  rm -f data/index.json
  python2 index_data.py $OUTDIR
  tar czf ~/work/exultant-pistachio/data/archives/bp_witness.$releaseid.tar.gz css index* js README.txt data/*.json --transform "s|^|bp-witness/|"
}

main
