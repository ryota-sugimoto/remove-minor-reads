#!/usr/bin/env bash

while getopts 't:' OPTION
do
  case $OPTION in
  t) THRESHOLD="-t ${OPTARG}" ;;
  esac
done
shift $((OPTIND - 1))

[ $# == 3 ] || { echo "${0} <bam/sam> <ref_fasta> <out_dir>"; exit 1; }
[ -f ${1} ] || { echo "${1} not found."; exit 1; }
[ -f ${2} ] || { echo "${2} not found."; exit 1; }
[ -d ${3} ] || { echo "${3} not found."; exit 1; }

BAM=$(cd $(dirname ${1}); pwd)/$(basename ${1})
REF_FASTA=$(cd $(dirname ${2}); pwd)/$(basename ${2})

WORKDIR=$(cd ${3}; pwd)
BAM_FILENAME=$(echo $(basename ${BAM}) | sed 's/\.[sb]am$//')


#sort sam
SORTED_SAM=${WORKDIR}/${BAM_FILENAME}.sorted.sam
SORTED_BAM=${WORKDIR}/${BAM_FILENAME}.sorted.bam
if [[ ${BAM} == *.sam ]]
then
  COMMAND_1=(samtools view -Shb ${BAM})
elif [[ ${BAM} == *.bam ]]
then
  COMMAND_1=(samtools view -hb ${BAM})
fi
${COMMAND_1[@]} \
  | samtools sort -o ${BAM} /dev/null \
  | tee ${SORTED_BAM} \
  | samtools view -h /dev/stdin > ${SORTED_SAM} || exit 1

#create pileup
PILEUP=${WORKDIR}/${BAM_FILENAME}.pileup
COMMAND=(samtools mpileup -f ${REF_FASTA}
                             ${SORTED_BAM})
${COMMAND[@]} > ${PILEUP} 2>/dev/null|| exit 1

#remove minor reads
CLEANED_BAM=${WORKDIR}/${BAM_FILENAME}.cleaned.bam
COMMAND_1=($(dirname ${0})/filter_hetero_reads.py ${THRESHOLD}
                                                  ${SORTED_SAM}
                                                  ${PILEUP})
${COMMAND_1[@]} \
  | samtools view -Shb /dev/stdin \
  | samtools sort -o /dev/stdin /dev/null > ${CLEANED_BAM} || exit 1

samtools index ${CLEANED_BAM}
rm ${PILEUP}
rm ${SORTED_BAM} ${SORTED_SAM}
