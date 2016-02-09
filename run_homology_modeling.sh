#!/bin/bash
#$ -S /bin/bash
#$ -N homology_modeling_setup
#$ -cwd
#$ -e logs
#$ -o logs
#$ -l h_vmem=20G
#$ -v FAST_CM
#$ -v ROSETTA_BIN
#$ -v ROSETTA_DB
#$ -v ROBETTA
##$ -p -1
#$ -tc 50 


## run with qsub -t 1:N run_homology_modeling.sh myseq.fasta
## where N is the values in the list that you want to model

export PYTHONPATH=$PYTHONPATH:/home/bertolan/geo
export PATH="/home/bertolan/envs/my_anaconda_clone/bin:$PATH"
export PATH=/home/bertolan/software/hmmer-3.1b1/bin:$PATH
export PATH=$PATH:/home/bertolan/software/hmmer-3.1b1/easel/miniapps

working_dir=$(pwd)
grep ">" $1 |cut -d ">" -f2 > list

infile=`awk 'NR==n' n=$SGE_TASK_ID $working_dir/list`
echo ~/envs/my_anaconda_clone/bin/python /home/bertolan/portable_scripts/new/homology_model_cm.py $1 $infile
~/envs/my_anaconda_clone/bin/python /home/bertolan/portable_scripts/new/homology_model_cm.py $1 $infile
