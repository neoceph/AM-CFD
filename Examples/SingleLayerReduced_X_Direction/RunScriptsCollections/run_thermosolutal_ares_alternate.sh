#!/bin/sh

# Set some input variables to determine input files to run 
inputNames="inputfile_MultiLayer.k"
logString="Log_MultiLayer"
export OMP_NUM_THREADS=4

# Provide appropiate paths
execDirec="/home/kkj7810/ThermalCFD/Thermosolutal_Play/ThermalCFD_Master/build/"
execFile="ThermosolutalCFD"
inputDirec="/home/aaa6979/test_ThermalCFD_Examples/testRun/InputFiles/"

# Define the name of the logfile
logAppend=".log"

# Set up the full file names/paths 
fileI=$inputNames
printf "Running %s file\n" $fileI
logFileName="$logString$logAppend"
inputI=$inputDirec$fileI
 
# run the simulation
cmd="$execDirec$execFile $inputI > $logFileName &"
eval $cmd
