#!/bin/sh
echo "Starting the simulation..."
################################################################################
#### make sure the nodes and processor numbers match in the following script####
################################################################################
inputNames="inputfile.k"
logString="Log_MultiLayer"
export OMP_NUM_THREADS=18

# ulimit -v 8,388,608,000
ulimit -v 8388608000

# Provide appropiate paths
execDirec="./"
execFile="../../../build/AM-CFD"
inputDirec="./"

# Define the name of the logfile
logAppend=".log"

# Set up the full file names/paths
fileI=$inputNames
printf "Running %s file\n" $fileI
logFileName="$logString$logAppend"
inputI=$inputDirec$fileI

# delete last analysis files
rm -R Calibration/
rm -R *.pvd *.log *.txt *.output
pkill -9 execFile

# run the simulation
cmd="$execDirec$execFile $inputI >> $logFileName 2>&1"
eval $cmd
