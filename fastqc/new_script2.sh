#!/bin/bash

runs=( '230221_JABBA_0141_AH3MCMDSX5/UnalignedI10I10SR'  '220729_JABBA_0099_AHF22GDMXY/UnalignedI8I8' '220804_JABBA_0101_BHF377DRX2/UnalignedI8I8' '230424_AVITIEB_0001_KOL-0472/Unaligned' '230505_JABBA_0160_AHTTFKDSX5/Unaligned10X' '230515_JABBA_0163_AHW2F7DRX2/UnalignedPE150I10I10'  ) 

for e in ${runs[@]}
do
	path=/work/PLTF/FBM/GTF/pltfgtf/pltf_gtf/ILLUMINA/DEMUX/$e/GTFQC/FastQC
	echo $path
	name=`echo $e | cut -d/ -f 1`
	echo $name
	python3 ./extract_fqc_metric2.py $path 'Sequence Duplication Levels' ./fastqc_data2/${name}_seq_dupl_lev.txt
        python3 ./extract_fqc_metric2.py $path	Adapter ./fastqc_data2/${name}_adapters.txt
	python3 ./extract_fqc_metric2.py $path 'Basic Statistics' ./fastqc_data2/${name}_basic_stats.txt
	python3 ./extract_fqc_metric2.py $path 'Per base sequence content' ./fastqc_data2/${name}_per_base_sequence_content.txt
	python3 ./extract_fqc_metric2.py $path 'Per base sequence quality' ./fastqc_data2/${name}_perBase_sequencing_quality.txt
	python3 ./extract_fqc_metric2.py $path  GC ./fastqc_data2/${name}_per_seq_gc.txt

done


#python3 extract_fqc_metric.py /work/PLTF/FBM/GTF/pltfgtf/pltf_gtf/ILLUMINA/DEMUX/230221_JABBA_0141_AH3MCMDSX5/UnalignedI10I10SR/GTFQC/FastQC 'Sequence Duplication Levels' ./230221_JABBA_0141_AH3MCMDSX5_seq_dupl_lev.txt
#python3 extract_fqc_metric.py /work/PLTF/FBM/GTF/pltfgtf/pltf_gtf/ILLUMINA/DEMUX/230221_JABBA_0141_AH3MCMDSX5/UnalignedI10I10SR/GTFQC/FastQC Adapter ./230221_JABBA_0141_AH3MCMDSX5_adapters.txt
#python3 extract_fqc_metric.py /work/PLTF/FBM/GTF/pltfgtf/pltf_gtf/ILLUMINA/DEMUX/230221_JABBA_0141_AH3MCMDSX5/UnalignedI10I10SR/GTFQC/FastQC 'Basic Statistics' ./230221_JABBA_0141_AH3MCMDSX5_basic_stats.txt
#python3 extract_fqc_metric.py /work/PLTF/FBM/GTF/pltfgtf/pltf_gtf/ILLUMINA/DEMUX/230221_JABBA_0141_AH3MCMDSX5/UnalignedI10I10SR/GTFQC/FastQC 'Per base sequence content' ./230221_JABBA_0141_AH3MCMDSX5_per_base_sequence_content.txt
#python3 extract_fqc_metric.py /work/PLTF/FBM/GTF/pltfgtf/pltf_gtf/ILLUMINA/DEMUX/230221_JABBA_0141_AH3MCMDSX5/UnalignedI10I10SR/GTFQC/FastQC 'Per base sequence quality' ./230221_JABBA_0141_AH3MCMDSX5_perBase_sequencing_quality.txt
#python3 extract_fqc_metric.py /work/PLTF/FBM/GTF/pltfgtf/pltf_gtf/ILLUMINA/DEMUX/230221_JABBA_0141_AH3MCMDSX5/UnalignedI10I10SR/GTFQC/FastQC GC ./230221_JABBA_0141_AH3MCMDSX5_per_seq_GC_cont.txt
