###############################################################
#
#
# Simpson et al. (2020) Methods 6 Quality Control of Sequencing Results
#
#       Input: Raw paired-end read files located in $SEQDATA
#
#       Output: FastQc files of reads before and after Trimmomatic along with Trimmed reads.
#               All output files are stored in $SEQDATA/
#
#       
#
################################################################


FASTQC=/home/shared/tools/FastQC
TRIMMOMATIC=/home/shared/tools/Trimmomatic-0.38
TrimGalore=/home/shared/tools/TrimGalore-0.060/trim_galore
PATH=$PATH:$FASTQC
PATH=$PATH:$TRIMMOMATIC
SEQDATA=/home/data/raw/ReefData/WaterSamples


cd $SEQDATA




# Running FastQc on the raw paired-end reads for Maycocks and output files

fastqc $SEQDATA/Maycocks_Reads/NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R1.fastq.gz -o $SEQDATA/Maycocks_Reads/

fastqc $SEQDATA/Maycocks_Reads/NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R2.fastq.gz -o $SEQDATA/Maycocks_Reads/

#Maycocks FastQc output files:
# 1) NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R1_fastqc.zip
# 2) NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R1_fastqc.html
# 3) NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R2_fastqc.zip
# 4) NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R2_fastqc.html



# Running FastQc on the raw paired end reads for Bellairs and output files. 

fastqc $SEQDATA/Bellairs_Reads/NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R1.fastq.gz -o $SEQDATA/Bellairs_Reads/

fastqc $SEQDATA/Bellairs_Reads/NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R2.fastq.gz -o $SEQDATA/Bellairs_Reads/

#Bellairs FastQc output files:
# 1) NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R1_fastqc.zip
# 2) NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R1_fastqc.html
# 3) NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R2_fastqc.zip
# 4) NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R2_fastqc.html



# Running Trimmomatic on the raw paired-end reads for Maycocks

java -jar trimmomatic-0.38.jar PE $SEQDATA/Maycocks_Reads/NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R1.fastq.gz $SEQDATA/Maycocks_Reads/NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R2.fastq.gz $SEQDATA/Maycocks_Reads/trimmed1.CS.p.fastq.gz $SEQDATA/Maycocks_Reads/trimmed1.CS.u.fastq.gz $SEQDATA/Maycocks_Reads/trimmed2.CS.p.fastq.gz $SEQDATA/Maycocks_Reads/trimmed2.CS.u.fastq.gz ILLUMINACLIP: NovaSeq.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Maycocks Trimmomatic Output files
# 1) trimmed1.CS.p.fastq.gz
# 2) trimmed1.CS.u.fastq.gz
# 3) trimmed2.CS.p.fastq.gz
# 4) trimmed2.CS.u.fastq.gz


# Running Trimmomatic on the raw paired-end reads for Bellairs

java -jar trimmomatic-0.38.jar PE $SEQDATA/Bellairs_Reads/NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R1.fastq.gz $SEQDATA/Bellairs_Reads/NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R2.fastq.gz trimmed1.BS.p.fastq.gz $SEQDATA/Bellairs_Reads/trimmed1.BS.u.fastq.gz  $SEQDATA/Bellairs_Reads/trimmed2.BS.u.fastq.gz ILLUMINACLIP: NovaSeq.fa:2:30:10 LEADING:3 TRAILING:3SLIDINGWINDOW:4:15 MINLEN:36

# Bellairs Trimmomatic Output files
# 1) trimmed1.BS.p.fastq.gz
# 2) trimmed1.BS.u.fastq.gz
# 3) trimmed2.BS.p.fastq.gz
# 4) trimmed2.BS.u.fastq.gz




# Running FastQc on trimmed and filtered paired-end reads for Maycocks and output files.

fastqc $SEQDATA/Maycocks_Reads/trimmed1.CS.p.fastq.gz -o $SEQDATA/Maycocks_Reads/

fastqc $SEQDATA/Maycocks_Reads/trimmed2.CS.p.fastq.gz -o $SEQDATA/Maycocks_Reads/

# Maycocks trimmed & filter FastQc output files: 
# 1) trimmed1.CS.p_fastqc.zip
# 2) trimmed1.CS.p_fastqc.html
# 3) trimmed2.CS.p_fastqc.zip
# 4) trimmed2.CS.p_fastqc.html


# Running FastQc on trimmed and filtered paired-end reads for Bellairs and output files.

fastqc $SEQDATA/Bellairs_Reads/trimmed1.BS.p.fastq.gz -o $SEQDATA/Bellairs_Reads/

fastqc $SEQDATA/Bellairs_Reads/trimmed2.BS.p.fastq.gz -o $SEQDATA/Bellairs_Reads/

# Bellairs trimmed & filter FastQc output files: 
# 1) trimmed1.BS.p_fastqc.zip
# 2) trimmed1.BS.p_fastqc.html
# 3) trimmed2.BS.p_fastqc.zip
# 4) trimmed2.BS.p_fastqc.html


# Running TrimGalore on raw paired-end reads for Maycocks and output files

TrimGalore-0.6.0/trim_galore --path_to_cutadapt ~/.local/bin/cutadapt --paired --phred33 --length 60 -q 5 --stringency 1 -e 0.1 -o $SEQDATA/Maycocks_trimgalore/ $SEQDATA/Maycocks_Reads/NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R1.fastq.gz $SEQDATA/Maycocks_Reads/NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R2.fastq.gz

# Maycocks TrimGalore output files: 
# 1) NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R1.fastq.gz_trimming_report.txt
# 2) NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R1_val_1.fq.gz
# 3) NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R2.fastq.gz_trimming_report.txt
# 4) NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R2_val_2.fq.gz




# Running TrimGalore on raw paired-end reads for Bellairs and output files 

TrimGalore-0.6.0/trim_galore --path_to_cutadapt ~/.local/bin/cutadapt --paired --phred33 --length 60 -q 5 --stringency 1 -e 0.1  -o $SEQDATA/Bellairs_trimgalore/ $SEQDATA/Bellairs_Reads/NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R1.fastq.gz $SEQDATA/Bellairs_Reads/NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R2.fastq.gz

# Bellairs TrimGalore output files:
# 1) NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R1.fastq.gz_trimming_report.txt
# 2) NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R1_val_1.fq.gz
# 3) NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R2.fastq.gz_trimming_report.txt
# 4) NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R2_val_2.fq.gz



# Running FastQc on TrimGalore paired-end reads for Maycocks and output files.

fastqc $SEQDATA/Maycocks_trimgalore/NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R1_val_1.fq.gz -o $SEQDATA/Maycocks_trimgalore/

fastqc $SEQDATA/Maycocks_trimgalore/NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R2_val_2.fq.gz -o $SEQDATA/Maycocks_trimgalore/

# Maycocks Fastqc of TrimGalore output files:
#1) NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R1_val_1_fastqc.html
#2) NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R1_val_1_fastqc.zip
#3) NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R2_val_2_fastqc.html
#4) NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R2_val_2_fastqc.zip





# Running FastQc on TrimGalore paired-end reads for Bellairs and output files.

fastqc $SEQDATA/Bellairs_trimgalore/NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R1_val_1.fq.gz -o $SEQDATA/Bellairs_trimgalore/

fastqc $SEQDATA/Bellairs_trimgalore/NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R2_val_2.fq.gz -o $SEQDATA/Bellairs_trimgalore/

# Bellairs Fastqc of TrimGalore output files:
# 1) NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R1_val_1_fastqc.html
# 2) NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R1_val_1_fastqc.zip
# 3) NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R2_val_2_fastqc.html
# 4) NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R2_val_2_fastqc.zip







