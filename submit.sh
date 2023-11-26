#!/bin/bash

#SBATCH --job-name=dynamice
#SBATCH --account=penny
#SBATCH --cpus-per-task=1
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

############################################################
# BASH SUBMIT
#
# Submit array of cluster tasks.
#
# Arguments:
#  JOB_TYPE: Type of job we're submitting
#  LOG_FILE: Text log file that stores IDs of completed tasks
#
# Various sbatch options can be set in options.R (see 'cluster settings' 
# section). For example, cluster partition, job memory, and slurm queue.
#
############################################################

# Load R
module purge
ml R/4.3.0-foss-2021a

# Extract inputs
job_type=$1
log_file=$2

# Extract array job ID
job_id=$(expr ${SLURM_ARRAY_TASK_ID})

# Run R script which calls simulation function
Rscript submit.R $job_type $job_id

# Write job ID to log file
echo "$job_id" >> $log_file

