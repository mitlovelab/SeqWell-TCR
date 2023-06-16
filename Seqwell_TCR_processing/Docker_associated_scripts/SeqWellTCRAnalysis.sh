#!/bin/bash
#### $1 = input text file containing location of aligned bam files and optional matching barcode file
#### $2 = # of desired threads
settings="/TCRAnalysis/bin/TCRsettings.txt"
bamFiles=($(awk -F $'\t' '{print "\""$1"\""}' $1 | sed 's/[[:space:]]*//g'))
species=($(awk -F $'\t' '{print "\""$2"\""}' $1 | sed 's/[[:space:]]*//g'))
BCFiles=($(awk -F $'\t' '{print "\""$3"\""}' $1 | sed 's/[[:space:]]*//g'))
timestamp="$(date +"%T")"
for ((i=0;i<${#bamFiles[@]};i++)); do
    /TCRAnalysis/bin/run_SeqWellTCRAnalysisKickoffv1.sh /usr/local/MATLAB/MATLAB_Compiler_Runtime/v84 ${bamFiles[$i]} ${species[$i]} ${BCFiles[$i]} $1 $settings $2 > log.out
    infoFile=($(grep 'infoFile =' log.out | sed 's/infoFile =//' | sed 's/[[:space:]]*//g'))
    echo $infoFile
    for ((k=0;k<$2;k++)); do
        /TCRAnalysis/bin/run_SeqWellTCRMapv1.sh /usr/local/MATLAB/MATLAB_Compiler_Runtime/v84 $infoFile $k &
    done
    wait
    /TCRAnalysis/bin/run_SeqWellTCRAnalysisAggregatev1.sh /usr/local/MATLAB/MATLAB_Compiler_Runtime/v84 $infoFile
done
