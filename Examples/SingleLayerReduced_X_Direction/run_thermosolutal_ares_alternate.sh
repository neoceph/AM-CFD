#!/bin/sh

# Set some input variables to determine input files to run
inputNames="inputfile_MultiLayer.k"
logString="Log_MultiLayer"
export OMP_NUM_THREADS=2

# Provide appropiate paths
execDirec="../../build/src/"
execFile="AM-CFD_run"
inputDirec="./resources/"

# Define the name of the logfile
logAppend=".log"

# Set up the full file names/paths
fileI=$inputNames
printf "Running %s file\n" $fileI
logFileName="$outputDirectory$logString$logAppend"
inputI=$inputDirec$fileI

# run the simulation
cmd="$execDirec$execFile $inputI > $logFileName &"
eval $cmd
