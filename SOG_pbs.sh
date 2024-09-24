#!/bin/bash
#PBS -q normal
#PBS -l mem=80GB
#PBS -l walltime=24:00:00
#PBS -l jobfs=100GB
#PBS -l ncpus=96
#PBS -P fa63
#PBS -l wd
#PBS -l storage=scratch/fa63+gdata/fa63
set -xue
module load python3
cd /g/data/fa63/zl1602/SOG
mkdir -p $PBS_JOBFS/SOG
cp -r genome $PBS_JOBFS/SOG
cp primer_tmp.tar $PBS_JOBFS/SOG
cp -r /g/data/fa63/zl1602/primersearch_tmp/ $PBS_JOBFS/SOG
cp  SOG_primersearch.py $PBS_JOBFS/SOG
cd $PBS_JOBFS/SOG
tar -xvf primer_tmp.tar
python3 SOG_primersearch.py >/g/data/fa63/zl1602/primersearch_tmp/SOG_primersearch.log 2>&1




cat *.report.csv|grep 'AU3_CHR'|sort -k11,11 -k4,4nr|awk '{if($4 >20)print $0}'|awk '{if($5 <=3) print $0}'