#!/bin/bash

################################################
#
# .bam file filtering 
#
################################################

Bam_File=$1
MAPQ=$2
Output_Folder=$3
Program_Folder=$4
Sample_Name=$5
Target_Name=$6
Assembly=$7




if [ ! -d "$Output_Folder/.tmp" ]; then
	mkdir $Output_Folder/.tmp
fi

mychr=$(cat $Program_Folder/data/targets/$Assembly/$Target_Name/*_chromosome.txt)



for i in $mychr
do
		grep -w $i $Program_Folder/data/targets/$Assembly/$Target_Name/TargetWindow.bed > $Output_Folder/.tmp/.temp.bed
		R --slave --args $Sample_Name,$Output_Folder,$Program_Folder,$i,$Target_Name,$Bam_File < $Program_Folder/lib/R/MakeDOC.R
done

