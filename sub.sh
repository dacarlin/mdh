#!/bin/sh
#
#SBATCH --job-name=mdh
#SBATCH --output=log.txt 

/share/work/alex/rosetta/source/bin/rosetta_scripts.linuxgccrelease @flags -user_tag $SLURM_ARRAY_TASK_ID

