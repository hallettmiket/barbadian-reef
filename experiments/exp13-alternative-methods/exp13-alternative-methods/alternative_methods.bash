###############################################################
#
#
# Simpson et al. (2020) Methods 6 Quality Control of Sequencing Results
#
#       Input: Trimmed & Filtered paired-end read files located in $SEQDATA
#
#       Output: Metaspades Files, Prokka Files, Diamond Files and BBMAP Files
#               All output files are stored in $RES/{metaspades, prokka, diamond, bbmap}
#
#       
#
################################################################


METASPADES=/home/shared/tools/MetaSpades/SPAdes-3.11.1-Linux/bin/
PROKKA=/home/shared/tools/prokka/bin/
DIAMOND=/home/shared/tools/Diamond/
BBMAP=/home/shared/tools/BBMAP/bbmap/
BOWTIE=/home/shared/tools/bowtie2-2.3.4.1/
SAMTOOLS=/home/shared/tools/samtools-1.9/
KRAKEN=/home/data/analysis-tools/kraken/
BRACKEN=/home/data/analysis-tools/Bracken/
PATH=$PATH:$METASPADES
PATH=$PATH:$PROKKA
PATH=$PATH:$DIAMOND
PATH=$PATH:$BBMAP
PATH=$PATH:$BOWTIE
PATH=$PATH:$SAMTOOlS
PATH=$PATH:$KRAKEN
PATH=$PATH:$BRACKEN
SEQDATA=/home/data/raw/ReefData/WaterSamples
RES=/home/data/refined/reef/
DDB=/home/shared/tools/Diamond/NCBI_DB/TaxNr.dmnd
MagpurifyDB=$RES/magpurify/db/


cd $METASPADES




# Running Metaspades on the trimmed & filter paired-end reads for Maycocks and output files

./metaspades.py -1 $SEQDATA/Maycocks_Reads/trimmed1.CS.p.fastq.qz -2 $SEQDATA/Maycocks_Reads/trimmed2.CS.p.fastq.qz -o $RES/metaspades/maycocks

#Maycocks Metasaapades output files:
# 1) assembly_graph.fastg
# 2) dataset.info
# 3) misc
# 4) params.txt
# 5) input_dataset.yaml
# 6) before_rr.fasta
# 7) assembly_graph_with_scaffolds.gfa
# 8) first_pe_contigs.fasta
# 9) contigs.fasta
# 10) contigs.paths 
# 11) corrected
# 12) K21
# 13) K33
# 14) K55
# 15) spades.log
# 16) scaffolds.fasta
# 17) tmp



# Running Metaspades on the trimmed & filter paired-end reads for Bellairs and output files 

./metaspades.py -1 $SEQDATA/Bellairs_Reads/trimmed1.BS.p.fastq.qz -2 $SEQDATA/Bellairs_Reads/trimmed2.BS.p.fastq.qz -o $RES/metaspades/bellairs 

#Bellairs Metaspades output files:
# 1) assembly_graph.fastg
# 2) dataset.info
# 3) misc
# 4) params.txt
# 5) input_dataset.yaml
# 6) before_rr.fasta
# 7) assembly_graph_with_scaffolds.gfa
# 8) first_pe_contigs.fasta
# 9) contigs.fasta
# 10) contigs.paths 
# 11) corrected
# 12) K21
# 13) K33
# 14) K55
# 15) spades.log
# 16) scaffolds.fasta
# 17) tmp


cd $PROKKA

# Running Prokka on the Metaspades assembled contigs for Maycocks and output files

prokka --outdir $RES/prokka/maycocks --centre may_sample $RES/metaspades/maycocks/contigs.fasta 


# Maycocks Prokka Output files
# 1) PROKKA_12212018.faa
# 2) PROKKA_12212018.fna
# 3) PROKKA_12212018.gff
# 4) PROKKA_12212018.tbl
# 5) PROKKA_12212018.txt
# 6) PROKKA_12212018.ffn
# 7) PROKKA_12212018.fsa
# 8) PROKKA_12212018.log
# 9) PROKKA_12212018.tsv


# Running Prokka on the Metaspades assembled contigs for Bellairs and output files

prokka --outdir $RES/prokka/bellairs --centre bel_sample $RES/metaspades/bellairs/contigs.fasta


# Bellairs Prokka Output files
# 1) PROKKA_12222018.faa
# 2) PROKKA_12222018.fna
# 3) PROKKA_12222018.gff
# 4) PROKKA_12222018.tbl
# 5) PROKKA_12222018.txt
# 6) PROKKA_12222018.ffn
# 7) PROKKA_12222018.fsa
# 8) PROKKA_12222018.log
# 9) PROKKA_12222018.tsv


cd $DIAMOND


# Running Diamond on prokka protein fasta file of the translated CDS sequences for Maycocks and output files.
 
diamond blastp -d $DDB -q $RES/prokka/maycocks/PROKKA_12212018.faa -o $RES/diamond/maycocks/Maycocks_faa_Blastp -f 5


# Maycocks Diamond output file: 
1) Maycocks_faa_Blastp





# Running Diamond on prokka protein fasta file of the translated CDS sequences for Bellairs and output file.

diamond blastp -d $DDB -q $RES/prokka/bellairs/PROKKA_12222018.faa -o $RES/diamond/bellairs/Bellairs_faa_Blastp -f 5

# Bellairs Diamond output files: 
1) Bellairs_faa_Blastp




################################## JUNE 12th ###################################

# Running Metaspades on the trimgalore paired-end reads for Maycocks and output files

./metaspades.py -1 $SEQDATA/Maycocks_trimgalore/NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R1_val_1.fq.gz -2 $SEQDATA/Maycocks_trimgalore/NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R2_val_2.fq.gz -o $RES/metaspades/trimgl_maycocks

#Maycocks Trimgalore Metaspades output files:
# 1) assembly_graph.fastg
# 2) dataset.info
# 3) misc
# 4) params.txt
# 5) input_dataset.yaml
# 6) before_rr.fasta
# 7) assembly_graph_with_scaffolds.gfa
# 8) first_pe_contigs.fasta
# 9) contigs.fasta
# 10) contigs.paths 
# 11) corrected
# 12) K21
# 13) K33
# 14) K55
# 15) spades.log
# 16) scaffolds.fasta
# 17) tmp



# Running Metaspades on the trimgalore paired-end reads for Bellairs and output files 

./metaspades.py -1 $SEQDATA/Bellairs_trimgalore/NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R1_val_1.fg.gz -2 $SEQDATA/Bellairs_trimgalore/NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R2_val_2.fq.gz -o $RES/metaspades/trimgl_bellairs


#Bellairs Trimgalore Metaspades output files :
# 1) assembly_graph.fastg
# 2) dataset.info
# 3) misc
# 4) params.txt
# 5) input_dataset.yaml
# 6) before_rr.fasta
# 7) assembly_graph_with_scaffolds.gfa
# 8) first_pe_contigs.fasta
# 9) contigs.fasta
# 10) contigs.paths 
# 11) corrected
# 12) K21
# 13) K33
# 14) K55
# 15) spades.log
# 16) scaffolds.fasta
# 17) tmp

#Running bowtie2 on the Maycocks assembled contigs from the trimgalore reads.
cd $BOWTIE

  #forming the reference index from the contigs
    bowtie2-build , $RES/metaspades/trimgl_maycocks/contigs.fasta  $RES/bowtie/contigs_index_maycocks
    
    #Outfile:
  # 1) *contigs_index_maycocks*.

    
 #mapping Maycocks raw reads back to their assemblies
    
    bowtie2 -p 12 -q -1 $SEQDATA/Maycocks_trimgalore/NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R1_val_1.fq.gz -2 $SEQDATA/Maycocks_trimgalore/NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R2_val_2.fq.gz -x $RES/bowtie/contigs_index_maycocks -S $RES/bowtie/maycocks.sam
  
  #Outfile:
  # 1) maycocks.sam

#Running bowtie2 on the Bellairs assembled contigs from the trimgalore reads.
cd $BOWTIE

  #forming the reference index from the contigs
    bowtie2-build , $RES/metaspades/trimgl_bellairs/contigs.fasta  $RES/bowtie/contigs_index_bellairs
    
  #Outputfile:
  # 1) *contigs_index_bellairs*.
    
 #mapping Bellairs raw reads back to their assemblies
   
   bowtie2 -p 12 -q -1 $SEQDATA/Bellairs_trimgalore/NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R1_val_1.fq.gz -2 $SEQDATA/Bellairs_trimgalore/NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R2_val_2.fq.gz -x $RES/bowtie/contigs_index_bellairs -S $RES/bowtie/bellairs.sam
   
   #Outfile:
  # 1) Bellairs.sam
  
#Converting  Maycocks Bowtie samfile to Bam using samtools

cd $RES/bowtie

 $SAMTOOLS/samtools sort maycocks.sam -@ 8 -o $RES/bowtie/maycocks.bam

  #Outfile:
  # 1) maycocks.bam
  
#Converting Bellairs Bowtie samfile to Bam using samtools

cd $RES/bowtie

$SAMTOOLS/samtools sort bellairs.sam -@ 8 -o $RES/bowtie/bellairs.bam

  #Outfile:
  # 1) bellairs.bam

#Binning Maycocks assemblies using MetaBat2 (within a conda environment)

cd $RES/bowtie

  #Creating Coverge depth file
  jgi_summarize_bam_contig_depths --outputDepth $RES/metabat/maycocks_depth.txt maycocks.bam
  
  #Outputfile:
  # 1) maycocks_depth.txt
  
  #Running Metabat2
  
  metabat2 -i  $RES/metaspades/trimgl_maycocks/contigs.fasta -a  $RES/metabat/maycocks_depth.txt -o  $RES/metabat/maycocks/bin
  
  #Outputfiles(Bins):
  # (1 bin.1.fa
  # (2 bin.2.fa
  # (3 bin.3.fa
  # (4 bin.4.fa
  # (5 bin.5.fa
  # (6 bin.6.fa
  # (7 bin.7.fa
  # (8 bin.8.fa
  # (9 bin.9.fa
  # (10 bin.10.fa
  # (11 bin.11.fa
  # (12 bin.12.fa
  # (13 bin.13.fa
  # (14 bin.14.fa
  # (15 bin.15.fa
  # (16 bin.16.fa
  # (17 bin.17.fa
  # (18 bin.18.fa
  # (19 bin.19.fa
  # (20 bin.20.fa
  # (21 bin.21.fa
  # (22 bin.22.fa
  # (23 bin.23.fa
  # (24 bin.24.fa
  # (25 bin.25.fa
  # (26 bin.26.fa
  # (27 bin.27.fa
  # (28 bin.28.fa
  # (29 bin.29.fa
  # (30 bin.30.fa
  
#Binning Bellairs assemblies using MetaBat2 (within a conda environment)  

cd $RES/bowtie

#Creating Coverge depth file
  jgi_summarize_bam_contig_depths --outputDepth $RES/metabat/bellairs_depth.txt bellairs.bam
  
  #Outputfile:
  # 1) bellairs_depth.txt
  
#Running Metabat2
  
  metabat2 -i  $RES/metaspades/trimgl_bellairs/contigs.fasta -a  $RES/metabat/bellairs_depth.txt -o  $RES/metabat/bellairs/bin
  
  #Outputfiles(Bins):
  # (1 bin.1.fa
  # (2 bin.2.fa
  # (3 bin.3.fa
  # (4 bin.4.fa
  # (5 bin.5.fa
  # (6 bin.6.fa
  # (7 bin.7.fa
  # (8 bin.8.fa
  # (9 bin.9.fa
  # (10 bin.10.fa
  # (11 bin.11.fa
  # (12 bin.12.fa
  # (13 bin.13.fa
  # (14 bin.14.fa
  # (15 bin.15.fa
  # (16 bin.16.fa
  # (17 bin.17.fa
  # (18 bin.18.fa
  
#Bin refinement of Maycocks Metabat2 bins using Magpurify via an internal bashscript that process each bin (within a conda environment)
cd $RES/metabat/maycocks
kl=1
for file in *.fa ; do
        dir="${file}"
        mkdir -p /home/data/refined/reef/magpurify/maycocks/"$dir"
        MAGpurify phylo-markers $RES/metabat/maycocks/$file $RES/magpurify/maycocks/"$dir"/
        MAGpurify clade-markers $RES/metabat/maycocks/$file $RES/magpurify/maycocks/"$dir"/
        MAGpurify tetra-freq  $RES/metabat/maycocks/$file $RES/magpurify/maycocks/"$dir"/
        MAGpurify gc-content $RES/metabat/maycocks/$file $RES/magpurify/maycocks/"$dir"/
        MAGpurify known-contam $RES/metabat/maycocks/$file $RES/magpurify/maycocks/"$dir"/
        MAGpurify clean-bin $RES/metabat/maycocks/$file $RES/magpurify/maycocks/"$dir"/ 
        echo "Bin number:$kl "
        kl=$((kl+1))
done

  
  
  
  #Outputfiles(Each bin (1-30) was represented by filename, within each bin displays the specific refinement):
  # 1) clade-markers
  # 2) cleaned_bin.fna
  # 3) gc-content 
  # 4) known-contam
  # 5) phylo-markers 
  # 6) tetra-freq
  
  
  
  
  
  
#Bin refinement of Bellairs Metabat2 bins using Magpurify (within a conda environment)
cd $RES/metabat/bellairs
kl=1
for file in *.fa ; do
        dir=${file}
        mkdir -p $RES/magpurify/maycocks/"$dir"
        MAGpurify phylo-markers $RES/metabat/bellairs/$file $RES/magpurify/bellairs/"$dir"/
        MAGpurify clade-markers $RES/metabat/bellairs/$file $RES/magpurify/bellairs/"$dir"/
        MAGpurify tetra-freq $RES/metabat/bellairs/$file $RES/magpurify/bellairs/"$dir"/
        MAGpurify gc-content $RES/metabat/bellairs/$file $RES/magpurify/bellairs/"$dir"/
        MAGpurify known-contam $RES/metabat/bellairs/$file $RES/magpurify/bellairs/"$dir"/
        MAGpurify clean-bin $RES/metabat/bellairs/$file $RES/magpurify/bellairs/"$dir"/  
        kl=$((kl+1))   
        echo "Bin number:$kl"
done
  

  #OutputFiles (Each bin (1-18) was represented by filename, within each bin displays the specific refinement):
  # 1) clade-markers
  # 2) cleaned_bin.fna
  # 3) gc-content 
  # 4) known-contam
  # 5) phylo-markers 
  # 6) tetra-freq
  
#Filtering incomplete bins with dRep and CheckM and the identify groups of organism with similar DNA content based on Average Nucleotide Identity from each bin's cleaned_bins.fna for Maycocks (within a conda environment) (Each bin's cleaned_bins.fna was placed in one folder with the respective bin# in the filename)

dRep dereplicate $RES/dRep/maycocks/drep_bins -g $RES/magpurify/maycocks/may_cleaned_bins/*.fna

#OutputFiles:
<<<<<<< HEAD
1) data
2) data_tables
3) dereplicated_genomes
4) figures
5) log




#Filtering incomplete bins with dRep and CheckM and the identify groups of organism with similar DNA content based on Average Nucleotide Identity from each bin's cleaned_bins.fna for Bellairs (within a conda environment) (Each bin's cleaned_bins.fna was placed in one folder with the respective bin# in the filename)


dRep dereplicate $RES/dRep/bellairs/drep_bins -g $RES/magpurify/bellairs/bell_cleaned_bins/*.fna  

#OutputFiles:

1) data
2) data_tables
3) dereplicated_genomes
4) figures
5) log


#Taxnomically Labelling cleaned bins using Kraken2(as outlined in exp4 script_bracken.bash)(Bellairs)


cd $KRAKEN

binDIR="$RES/magpurify/bellairs/bell_cleaned_bins"
files=`find $binDIR -name "*.fna"`
echo $files

for i in $files
do
out=`basename "${i%%.fna}"`
kraken2 --db $RES/kr_br $i \
--report $RES/kr_br/assemblies_krbr/bellairs/$out.kreport > $RES/kr_br/assemblies_krbr/bellairs/$out.kraken
done

#Outputfiles:
#Kraken files and reports for all cleaned bins



#Taxnomically Labelling cleaned bins using Kraken2(as outlined in exp4 script_bracken.bash)(Maycocks)

cd $KRAKEN
binDIR="$RES/magpurify/maycocks/may_cleaned_bins"
files=`find $binDIR -name "*.fna"`
echo $files

for i in $files
do
out=`basename "${i%%.fna}"`
kraken2 --db $RES/kr_br $i \
--report $RES/kr_br/assemblies_krbr/maycocks/$out.kreport > $RES/kr_br/assemblies_krbr/maycocks/$out.kraken
done


=======
