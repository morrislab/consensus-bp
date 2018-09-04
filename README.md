ICGC PCAWG-11 consensus breakpoint (segmentation) pipeline
==========================================================

General procedure
-----------------
This pipeline uses segmentations produced by six callers (ABSOLUTE, Aceseq,
Battenberg, cloneHD, JaBbA and Sclust) to create consensus segmentations
consisting of breakpoints that denote where in the genome copy number status
changes for each patient.

Independent processes are run for each sample, permitting the pipeline to be
parallelized to decrease execution time. IO speed affects runtime because the
pipeline writes several thousand files to disk. On a 20-core Xeon E5-2630 with
128 GB of RAM and a fast SSD, the pipeline finishes in approximately five
minutes.

Dependencies
------------
* Python 2.7 (2.7.12 used)
* NumPy (1.11.0 used)
* Plotly (optional, only required for some plots) (2.0.14 used)
* SciPy (optional, only required for some plots) (0.17.0 used)

Running the pipeline
--------------------
1. Create a base directory `$BASEDIR` from which to run the analysis. For this, I will use `BASEDIR=~/consensus_bp`.

2. Open `run.sh` in a text editor and change `BASEDIR` to reference `~/consensus_bp`.

3. Download `consensus_bp_data.20171111.tar.gz` from Synapse at https://www.synapse.org/#!Synapse:syn16580326.

4. Change directories to `$BASEDIR`. Extract `consensus_bp_data.20171111.tar.gz` to create the subdirectory `data/`.

```
cd ~/consensus_bp
tar xzf ~/Downloads/consensus_bp_data.20171111.tar.gz
```

5. Run `bash run.sh`.

Outputs
-------
Upon pipeline completion, results will be available in
`$BASEDIR/results/archives`. Three files are written:

1. `consensus_bp.basic.<YYYYMMDD>.tar.gz`: this archive contains `.txt` files
   for each sample listing consensus breakpoints in a two-column format, with
   columns corresponding to chromosome and locus, respectively. The first line
   is a header.

2. `consensus_bp.extended.<YYYYMMDD>.tar.gz`: this archive contains four
   directories corresponding to autosome, female X, male X, and male Y
   breakpoints. Each directory in turn contains JSON files for each patient,
   with all patients having files in the autosome directory, and the other
   three directories segregated by sex. These JSON files contain more in-depth
   information concerning how each breakpoint in the `basic` archive was
   derived from the inputs. There are also `.stderr` files for each patient
   containing log messages.

3. `consensus_bp_methods.<YYYYMMDD>.tar.gz`: this archive contains a
   tab-delimited file listing what input methods were used for each stage of
   the consensus (autosome, male X, male Y, female X). Different methods and
   agreement thresholds were used for each set, as some methods could not
   produce outputs for some of these stages.
