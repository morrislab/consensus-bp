#!/bin/bash
set -euo pipefail

#export LC_ALL=C
#module load gnu-parallel/20150822

BASEDIR=~/work/exultant-pistachio
PROTDIR=$BASEDIR/protocols/compare-breakpoints

SVDIR=$BASEDIR/data/sv
WHITELIST=$BASEDIR/data/misc/whitelist.20160831.txt
GREYLIST=$BASEDIR/data/misc/greylist.20160831.txt
CENTROMERES=$BASEDIR/data/misc/cytoBand.txt.gz
SEX=$BASEDIR/data/misc/inferred_sex.all_samples.20161209.txt
WINDOW_SIZE=100000
INPUT_SEGS=$BASEDIR/data/cnvs.pre_consensus_bp

PARALLEL=40
OUTDIR=$BASEDIR/data/consensus_bp
MERGEDDIR=$OUTDIR/bp.merged
releaseid=$(date '+%Y%m%d')
sample_list=$BASEDIR/data/archives/consensus_bp_methods.$releaseid.txt

METHODS_AUTOSOMES="broad dkfz jabba mustonen095 peifer vanloo_wedge_segs"
METHODS_X_MALE="broad_x dkfz jabba mustonen095 peifer"
METHODS_X_FEMALE="broad_x dkfz jabba mustonen095 peifer vanloo_wedge_segs"
METHODS_Y="dkfz jabba"

function make_consensus {
  rm -rf $OUTDIR && mkdir -p $OUTDIR

  cd $INPUT_SEGS
  (
  python2 $PROTDIR/run_comparison.py \
    --whitelist $WHITELIST \
    --greylist $GREYLIST \
    --window-size $WINDOW_SIZE \
    --centromere-filename $CENTROMERES \
    --sex-filename $SEX \
    --mode autosomes \
    $SVDIR \
    $OUTDIR \
    $METHODS_AUTOSOMES
  python2 $PROTDIR/run_comparison.py \
    --whitelist $WHITELIST \
    --greylist $GREYLIST \
    --window-size $WINDOW_SIZE \
    --centromere-filename $CENTROMERES \
    --sex-filename $SEX \
    --mode X_female \
    $SVDIR \
    $OUTDIR \
    $METHODS_X_FEMALE
  python2 $PROTDIR/run_comparison.py \
    --whitelist $WHITELIST \
    --greylist $GREYLIST \
    --window-size $WINDOW_SIZE \
    --centromere-filename $CENTROMERES \
    --sex-filename $SEX \
    --mode X_male \
    $SVDIR \
    $OUTDIR \
    $METHODS_X_MALE
  python2 $PROTDIR/run_comparison.py \
    --whitelist $WHITELIST \
    --greylist $GREYLIST \
    --window-size $WINDOW_SIZE \
    --centromere-filename $CENTROMERES \
    --sex-filename $SEX \
    --mode Y \
    $SVDIR \
    $OUTDIR \
    $METHODS_Y
  ) | parallel -j$PARALLEL
}

function merge_releases {
  cd $OUTDIR
  rm -rf $MERGEDDIR && python2 $PROTDIR/merge_releases.py $SEX bp.merged bp.{autosomes,X_female,X_male,Y} > $sample_list
}

function package_release {
  cd $MERGEDDIR/data
  tar czf $BASEDIR/data/archives/consensus_bp.$releaseid.tar.gz *-*.txt --transform "s|^|consensus_bp.$releaseid.tar.gz/|"
  tar czf $BASEDIR/data/archives/consensus_bp.extended.$releaseid.tar.gz *-*.json --transform "s|^|consensus_bp.extended.$releaseid.tar.gz/|"
}

function main {
  make_consensus
  merge_releases
  package_release
  exit

  rm -rf too_few && mkdir too_few
  for foo in $(cat $sample_list | python2 $PROTDIR/filter.py | cut -d ' ' -f 1); do
    mv $foo* too_few
  done

  cd ~/work/bp-witness
  rm -f data/index.json
  python2 index_data.py $OUTDIR
  tar czf $BASEDIR/data/archives/bp_witness.$releaseid.tar.gz css index* js README.txt data/*.json --transform "s|^|bp-witness/|"
}

main
