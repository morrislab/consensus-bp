#!/bin/bash
set -euo pipefail

export LC_ALL=C
module load gnu-parallel/20150822
OUTDIR=~/work/exultant-pistachio/data/consensus_bp
SVDIR=~/work/exultant-pistachio/data/sv
BLACKLIST=~/work/exultant-pistachio/data/misc/blacklist.20160906.txt
CENTROMERES=~/work/exultant-pistachio/data/misc/cytoBand.txt.gz

cd ~/work/exultant-pistachio/data/cnvs

rm -f $OUTDIR/*.{json,txt,stderr}
python2 ~/work/exultant-pistachio/protocols/compare-breakpoints/run_comparison.py \
  --blacklist $BLACKLIST \
  --window-size 100000 \
  --centromere-filename $CENTROMERES \
  --optional-methods vanloo_wedge_segs,mustonen095,peifer,broad,dkfz \
  $SVDIR \
  $OUTDIR \
  broad dkfz mustonen095 peifer vanloo_wedge_segs \
  | parallel -j16

today=$(date '+%Y%m%d')
cd $OUTDIR
for foo in *.json; do echo $(echo $foo | cut -d . -f1)$'\t'$(cat $foo | jq -r '.methods|sort|join(",")'); done > ~/work/exultant-pistachio/data/archives/consensus_bp_methods.$today.txt
cd ~/work/exultant-pistachio/data
tar czf ~/work/exultant-pistachio/data/archives/consensus_bp.$today.tar.gz consensus_bp/*.txt

cd ~/work/bp-witness
rm -f data/index.json
python2 index_data.py
tar czf ~/work/exultant-pistachio/data/archives/bp_witness.$today.tar.gz css index* js README.txt data/*.json --transform "s|^|bp-witness/|"
