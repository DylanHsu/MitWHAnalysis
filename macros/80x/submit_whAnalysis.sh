#!/bin/bash
function splitPath {
  local IFS=/
  local pieces
  # Disable glob expansion
  set -f
  pieces=( $@ ) 
  set +f
  #printf '%d\n' "${#pieces[@]}"
  #printf '%s\n' "${pieces[@]}"
  echo ${pieces[${#pieces[@]}-1]}
}
jobsPerSkim=100
receptacle=/data/t3home000/${USER}/whAnalysisPlots
weights_dir=/home/dhsu/CMSSW_8_0_26_patch1/src/weights


lastJob=`expr ${jobsPerSkim} - 1`
mkdir -p ${receptacle}/${1}
to_run="/home/dhsu/bin/condor-run ./batch_whAnalysis.sh"
weightsFileBase=$3
libs="${weights_dir}/${weightsFileBase} ${CMSSW_BASE}/src/MitWHAnalysis/macros/80x/batch_whAnalysis.sh ${CMSSW_BASE}/src/MitWHAnalysis/macros/80x/whAnalysis.C ${CMSSW_BASE}/src/MitWHAnalysis/macros/80x/whAnalysis.h"
while read line
do
    cd ${receptacle}/${1}
    options=( $line )
    batchFile=${options[0]}
    batchCategory=${options[1]}
    batchFileBase=`splitPath ${batchFile}`
    echo "Setting up ${jobsPerSkim} jobs with input files: ${libs} ${batchFile}"
    echo "Will pick up output files wh_batchTree_*_${batchFileBase}"
    for job in `seq 0 ${lastJob}`
    do
      outputFile="wh_batchTree_${job}_${batchFileBase}"
      args="--auxiliary-input ${libs} ${batchFile} --output ${outputFile}"
      args="${args} --task-name $1"
      args="${args} --job-args \"${batchFileBase} ${batchCategory} ${weightsFileBase} ${job} ${jobsPerSkim}\" "
      #echo "$to_run $args"
      eval "$to_run $args"
    done
done < "$2"
cd ${CMSSW_BASE}/src
