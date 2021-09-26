#!/bin/bash

if [ -z "$4" ]; then
	echo
	echo "Program: Initialize reference files for XCAVATOR"
	echo
	echo "Usage: bash $0 <assembly> <target-name> <window-size> <ref.fa> [ref.bw]"
	echo
	echo "  e.g: bash $0 hg19 bin-500kb 500000 /path/to/ucsc.hg19.fasta"
	echo
	exit 1
fi

ASSEMBLY=$1
TARGET_NAME=$2
WINDOW_SIZE=$3
REF_FA=$4
REF_BW=$5

#----------------------------------------------------------#

APP_DIR=$(dirname $(realpath $0))

if [ -z "${REF_BW}" ]; then
	if [ "${ASSEMBLY}" == "hg19" ]; then
		REF_BW=${APP_DIR}/data/ucsc.hg19.bw
	else
		echo 1>&2 "Error: <ref.wig> file should be specified when <assembly> is not 'hg19'!"
		exit 1
	fi
fi

for FILE in "${APP_DIR}/ReferenceWindowInitialize.pl" "${REF_FA}" "${REF_BW}"; do
	if [ ! -e "${FILE}" ]; then
		echo 1>&2 "Error: File '${FILE}' does not exists!"
		exit 1
	fi
done

#----------------------------------------------------------#

if [ -e "${APP_DIR}/data/targets/${ASSEMBLY}/${TARGET_NAME}/SourceTarget.txt" ]; then
	echo 1>&2 "Error: XCAVATOR reference '${ASSEMBLY}/${TARGET_NAME}' has existed!"
	exit 1
fi

#----------------------------------------------------------#

mkdir -pv ${APP_DIR}/data/targets/${ASSEMBLY}/${TARGET_NAME}/log/

echo 1>&2 "Generate '${APP_DIR}/data/targets/${ASSEMBLY}/${TARGET_NAME}/SourceTarget.txt'"
echo "${REF_BW} ${REF_FA}" > ${APP_DIR}/data/targets/${ASSEMBLY}/${TARGET_NAME}/SourceTarget.txt

#----------------------------------------------------------#

LOG_PREFIX="${APP_DIR}/data/targets/${ASSEMBLY}/${TARGET_NAME}/log/$(date '+%Y%m%d.%H%M%S').$$.ref-win-init"

echo 1>&2 "Run 'ReferenceWindowInitialize.pl' (${ASSEMBLY}, ${TARGET_NAME}, ${WINDOW_SIZE})"
perl ${APP_DIR}/ReferenceWindowInitialize.pl \
	${APP_DIR}/data/targets/${ASSEMBLY}/${TARGET_NAME}/SourceTarget.txt \
	${TARGET_NAME} ${WINDOW_SIZE} ${ASSEMBLY} \
	>${LOG_PREFIX}.log 2>${LOG_PREFIX}.err

if [ "$?" != "0" ]; then
	exit $?
fi

#----------------------------------------------------------#
echo 1>&2 "done."
