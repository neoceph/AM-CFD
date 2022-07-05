#!/bin/sh
#PBS -V
#PBS -N ShowChk
#PBS -q batch
#PBS -l nodes=node020:ppn=1
#PBS -l walltime=120:00:00
#PBS -I
#PBS -X

echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR 
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo The following processors are allocated to this job:
echo `cat $PBS_NODEFILE`
NP=`wc -l < $PBS_NODEFILE`

# These environment variables are available in every batch job
#
# $PBS_ENVIRONMENT set to PBS_BATCH to indicate that the job is a batch job; otherwise,
#                  set to PBS_INTERACTIVE to indicate that the job is a PBS interactive job
# $PBS_JOBID       the job identifier assigned to the job by the batch system
# $PBS_JOBNAME     the job name supplied by the user
# $PBS_NODEFILE    the name of the file that contains the list of nodes assigned to the job
# $PBS_QUEUE       the name of the queue from which the job is executed
# $PBS_O_HOME      value of the HOME variable in the environment in which qsub was executed
# $PBS_O_LANG      value of the LANG variable in the environment in which qsub was executed
# $PBS_O_LOGNAME   value of the LOGNAME variable in the environment in which qsub was executed
# $PBS_O_PATH      value of the PATH variable in the environment in which qsub was executed
# $PBS_O_MAIL      value of the MAIL variable in the environment in which qsub was executed
# $PBS_O_SHELL     value of the SHELL variable in the environment in which qsub was executed
# $PBS_O_TZ        value of the TZ variable in the environment in which qsub was executed
# $PBS_O_HOST      the name of the host upon which the qsub command is running
# $PBS_O_QUEUE     the name of the original queue to which the job was submitted
# $PBS_O_WORKDIR   the absolute path of the current working directory of the qsub command
#
# End of example PBS script
