#!/usr/bin/env bash
 
#SBATCH -c 4                  # number of CPUs (here 4)
#SBATCH --mem=32000          # memory pool for all cores (here 32GB)
#SBATCH --partition=short      # Partition to submit to
#SBATCH -o slurm.%j.out        # STDOUT (the standard output stream)
#SBATCH -e slurm.%j.err        # STDERR (the output stream for errors)
#SBATCH --mail-type=fail       # notifications for job to abort
#SBATCH --mail-type=time_limit # notifications for job to abort
#SBATCH --mail-use=fagernaes@shh.mpg.de # these notifications will be sent to your shh.mpg.de email-address.
#SBATCH --array=0-2%3         # The number of jobs in the array (first job indicated as 0, so here we are submitting 23 jobs), and how many to submit at once (after %, so in here we are running 4 jobs in parallel)
#SBATCH -J "EAGER"             # Job name
 
SAMPLES=( $(find /projects1/microbiome_calculus/pacific_calculus/03-Preprocessing/eager_HG19/*/ -name '2020-05-23-14-13-EAGER.xml' -type f) ) # The outermost set of brackets defines this as a bash "list". This steps finds all the xml files in /PATH/TO/EAGER_RESULTS/ and save them as a bash "list"
SAMPLENAME=${SAMPLES[$SLURM_ARRAY_TASK_ID]} #The variable is being set every time that the SLURM_ARRAY_TASK_ID changes, so it is going through the list

unset DISPLAY # This will empty the DISPLAY variable, which can make DnaDamage throw an error and stop execution.

eagercli "${SAMPLENAME}" #This submits an eager job with the current xml from the list
