#!/usr/bin/env bash

while getopts 't:' OPTION
do
  case $OPTION in
  t) THRESHOLD="-t ${OPTARG}" ;;
  esac
done
shift $((OPTIND - 1))

[ $# == 3 ] || { echo "${0} <bam/sam> <ref_fasta>"; exit 1; }
[ -f ${1} ] || { echo "${1} not found."; exit 1; }
[ -f ${2} ] || { echo "${2} not found."; exit 1; }
[ -d ${3} ] || { echo "${3} not found."; exit 1; }

BAM=$(cd $(dirname ${1}); pwd)/$(basename ${1})
REF_FASTA=$(cd $(dirname ${2}); pwd)/$(basename ${2})

PICARD_SORTSAM="/usr/local/share/java/SortSam.jar"

WORKDIR=$(cd ${3}; pwd)
BAM_FILENAME=$(echo $(basename ${BAM}) | sed 's/\.[sb]am$//')

echo ${WORKDIR} ${BAM_FILENAME}

#sort sam
SORTED_SAM=${WORKDIR}/${BAM_FILENAME}.sorted.sam
SORTED_BAM=${WORKDIR}/${BAM_FILENAME}.sorted.bam

COMMAND=(java -jar ${PICARD_SORTSAM}
         SORT_ORDER=coordinate
         CREATE_INDEX=true
         INPUT=${BAM}
         OUTPUT=${SORTED_BAM}
         VALIDATION_STRINGENCY=LENIENT)
${COMMAND[@]} || exit 1
samtools view -h ${SORTED_BAM} > ${SORTED_SAM}
echo

#create pileup
PILEUP=${WORKDIR}/${BAM_FILENAME}.pileup
COMMAND=(samtools mpileup -f ${REF_FASTA}
                             ${SORTED_BAM})
${COMMAND[@]} > ${PILEUP} 2>/dev/null|| exit 1

#remove minor reads
CLEANED_SAM=${WORKDIR}/${BAM_FILENAME}.cleaned.sam
CLEANED_BAM=${WORKDIR}/${BAM_FILENAME}.cleaned.bam
COMMAND_1=($(dirname ${0})/remove_minor_reads.py ${THRESHOLD}
                                                 ${SORTED_SAM}
                                                 ${PILEUP})
COMMAND_2=(java -jar ${PICARD_SORTSAM}
           SORT_ORDER=coordinate
           CREATE_INDEX=true
           INPUT=/dev/stdin
           OUTPUT=${CLEANED_BAM}
           VALIDATION_STRINGENCY=LENIENT)
${COMMAND_1[@]} | ${COMMAND_2[@]} || exit 1

rm ${PILEUP}
rm ${SORTED_BAM} ${SORTED_BAM%m}i
rm ${SORTED_SAM}
