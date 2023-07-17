#!/bin/bash

export PATH=/fshare2/Softwares/cellranger/cellranger-3.0.1/cellranger-cs/3.0.1/bin/:$PATH
path_to_reference_genome="/fshare2/sharedData/Single_Cell/10X"
fastq_path="/Sshare/home/lanLessonPublic/scRNA-seq/1.0.rawdata/pbmc_1k_v3_fastqs"

cellranger count \
	--id=pbmc_1k_v3 \
       	--transcriptome=${path_to_reference_genome}/refdata-cellranger-GRCh38-3.0.0 \
       	--fastqs=${fastq_path} \
	--sample=pbmc_1k_v3 \
	 --localmem=50 \
	 --localcores=4
