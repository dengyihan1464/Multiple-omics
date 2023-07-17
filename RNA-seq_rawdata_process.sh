#!/bin/bash

# *************************
# FastQC
# *************************
export PATH=/fshare2/JiFansen/Software/Fastqc/FastQC/:$PATH
mkdir 2.fastQC
inputData_dir="/Sshare/home/lanLessonPublic/RNA-seq/1.rawdata/"
mkdir tmp

for sample_name in HBR_Rep1 HBR_Rep2 HBR_Rep3 
do
	fastqc -d ./tmp ${inputData_dir}${sample_name}_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq -O 2.fastQC
	fastqc -d ./tmp ${inputData_dir}${sample_name}_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq -O 2.fastQC
done

for sample_name in UHR_Rep1 UHR_Rep2 UHR_Rep3
do
	fastqc -d ./tmp ${inputData_dir}${sample_name}_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq -O 2.fastQC
	fastqc -d ./tmp ${inputData_dir}${sample_name}_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq -O 2.fastQC
done
# -------------------------------------------------------------------
# **************************
# Reads Trimming
# **************************
mkdir 3.trimReads
export trimmomatic=/Sshare/home/lanLessonPublic/Softwares/Trimmomatic-0.33/trimmomatic-0.33.jar
export adaptor_file=/Sshare/home/lanLessonPublic/Softwares/Trimmomatic-0.33/test.fa
cd 3.trimReads
mkdir -p paired_fastq
mkdir -p unpaired_fastq
HEAD=13

for sample_name in HBR_Rep1 HBR_Rep2 HBR_Rep3
do
	java -jar $trimmomatic \
		PE -phred33 \
		${inputData_dir}${sample_name}_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq \
		${inputData_dir}${sample_name}_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq \
		paired_fastq/${sample_name}_R1.fastq unpaired_fastq/${sample_name}_R1.fastq \
		paired_fastq/${sample_name}_R2.fastq unpaired_fastq/${sample_name}_R2.fastq \
		ILLUMINACLIP:${adaptor_file}:2:30:10 \
		HEADCROP:${HEAD} \
		SLIDINGWINDOW:4:15 \
		MINLEN:65
done

for sample_name in UHR_Rep1 UHR_Rep2 UHR_Rep3
do
	        java -jar $trimmomatic \
			PE -phred33 \
			${inputData_dir}${sample_name}_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq \
			${inputData_dir}${sample_name}_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq \
			paired_fastq/${sample_name}_R1.fastq unpaired_fastq/${sample_name}_R1.fastq \
			paired_fastq/${sample_name}_R2.fastq unpaired_fastq/${sample_name}_R2.fastq \
			ILLUMINACLIP:${adaptor_file}:2:30:10 \
			HEADCROP:${HEAD} \
			SLIDINGWINDOW:4:15 MINLEN:65
done
cd ../
# -----------------------------------------------
 ****************
 Alignment
 ****************
export PATH=/fshare2/JiFansen/Software/STAR/STAR-2.5.3a/bin/Linux_x86_64/:$PATH
mkdir 4.alignment
GENOME="/fshare2/JiFansen/1.1.Reference/2.4.ReferenceGenome/Ensemble_GRCh37"
GENOME_FA="/fshare2/JiFansen/1.1.Reference/2.4.ReferenceGenome/Ensemble_GRCh37/Ensemble_GRCh37_genome.fa"
GENOME_GTF="/fshare2/JiFansen/1.1.Reference/2.4.ReferenceGenome/Ensemble_GRCh37/genes.gtf"
STAR_INDEX="$GENOME/star/STAR_2.4.1c/"
N_CPUS=4

for sample_name in HBR_Rep1 HBR_Rep2 HBR_Rep3 UHR_Rep1 UHR_Rep2 UHR_Rep3
do
	STAR \
		--genomeDir $STAR_INDEX \
		--sjdbGTFfile $GENOME_GTF \
		--runThreadN $N_CPUS \
		--outSAMstrandField intronMotif \
		--outFilterIntronMotifs RemoveNoncanonical \
		--outFileNamePrefix 4.alignment/$sample_name \
		--readFilesIn \
		3.trimReads/paired_fastq/${sample_name}_R1.fastq \
		3.trimReads/paired_fastq/${sample_name}_R2.fastq \
		--outSAMtype BAM SortedByCoordinate \
		--outReadsUnmapped Fastx \
		--outSAMmode Full
done
# -------------------------------------------------
# **********************
# Counts Quantification
# **********************
export PATH=/fshare2/JiFansen/Software/subread/subread-1.6.3-Linux-x86_64/bin/:$PATH
mkdir 5.featureCounts
suffix="Aligned.sortedByCoord.out.bam"

for sample_name in HBR_Rep1 HBR_Rep2 HBR_Rep3 UHR_Rep1 UHR_Rep2 UHR_Rep3
do
	outname="$sample_name.count.txt"
	bam="4.alignment/$sample_name$suffix"
	featureCounts \
		-T $N_CPUS \
		-t exon \
		-g gene_id \
		-a $GENOME_GTF \
		-o 5.featureCounts/$outname \
		$bam
	cut -f 1,7 5.featureCounts/$outname |grep -v '^#' >5.featureCounts/${sample_name}_featureCounts.txt
done

cd 5.featureCounts
paste *featureCounts.txt > Counts.txt
