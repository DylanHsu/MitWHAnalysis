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

receptacle=/data/t3home000/${USER}/whAnalysisPlots
mkdir -p ${receptacle}/${1}
weights_dir=/home/dhsu/CMSSW_8_0_26_patch1/src/weights
to_run="/home/dhsu/bin/condor-run ./batch_whAnalysis.sh"
weightsFileBase=$3
libs="${weights_dir}/${weightsFileBase} ${CMSSW_BASE}/src/MitWHAnalysis/macros/80x/batch_whAnalysis.sh ${CMSSW_BASE}/src/MitWHAnalysis/macros/80x/whAnalysis.C"
while read line
do
    cd ${receptacle}/${1}
    options=( $line )
    batchFile=${options[0]}
    batchCategory=${options[1]}
    batchFileBase=`splitPath ${batchFile}`
    outputFile="histo_wh_nice_${batchFileBase}"
    echo "Setting up job with input files: ${libs} ${batchFile}"
    echo "Will pick up output file ${outputFile}"

    args="--auxiliary-input ${libs} ${batchFile} --output histo_wh_nice_${batchFileBase}"
    args="${args} --task-name $1"
    args="${args} --job-args \"${batchFileBase} ${batchCategory} ${weightsFileBase}\" "
    #echo "$to_run $args"
    eval "$to_run $args"
done < "$2"
cd ${CMSSW_BASE}/src
