#!/bin/bash
#$ -P moru-batty.prjc
#$ -wd /well/moru-batty/users/gka292/PLATCOV/Ivermectin
#$ -N res
#$ -pe shmem 4
#$ -o /well/moru-batty/users/gka292/PLATCOV/Ivermectin/o_and_e_files
#$ -e /well/moru-batty/users/gka292/PLATCOV/Ivermectin/o_and_e_files
#$ -q short.qf
#$ -t 1-9
#$ -tc 32

echo started=`date`
module purge
module load R/4.1.2-foss-2021b

echo "job=$JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`
Rscript /well/moru-batty/users/gka292/PLATCOV/Ivermectin/run_models.R ${SGE_TASK_ID} --no-save --no-restore
echo "finished="`date`
exit 0
