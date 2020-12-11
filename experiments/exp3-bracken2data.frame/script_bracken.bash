###############################################################
#
#
# Simpson et al. (2020) Methods 7  Read alignment and taxa classification
#
#       Input: post-Trimmomatic read files located in $SEQDATA
#
#       Output: Kraken indexed database followed by Bracken assignments to taxa.
#               Bracken output files stored in $RES/{bellairs.bracken, maycocks.bracken}
#
#       These output files are used by the R routines in exp4/bracken2data.frame.R to produce an R datastructure
#
################################################################

PATH=$PATH:~/bin:~/ncbi-blast-2.9.0+/bin/
KRAKEN=/home/data/analysis-tools/kraken
BRACKEN=/home/data/analysis-tools/Bracken
PATH=$PATH:$KRAKEN
PATH=$PATH:$BRACKEN
SEQDATA=/home/data/raw/ReefData/WaterSamples
RES=/home/data/refined/reef/kr_br
KDB=$RES

cd $RES

# Download list of genomes to be included in the Kraken Mar reference
wget -nv -O dl_list_mardb_mags.txt https://s1.sfb.uit.no/public/mar/Resources/kraken/dl_list_mardb_mags.txt
wget -nv -O dl_list_mardb_no_mags.txt https://s1.sfb.uit.no/public/mar/Resources/kraken/dl_list_mardb_no_mags.txt
wget -nv -O dl_list_marref.txt https://s1.sfb.uit.no/public/mar/Resources/kraken/dl_list_marref.txt

# Check that downloads succeeded
[ -r dl_list_mardb_mags.txt ] || { echo Download of dl_list_mardb_mags.txt has failed; exit 1; }
[ -r dl_list_mardb_no_mags.txt ] || { echo Download of dl_list_mardb_no_mags.txt has failed; exit 1; }
[ -r dl_list_marref.txt ] || { echo Download of dl_list_marref.txt has failed; exit 1; }

# Download all files from the MAR suite of marine databases to subfolder marine_micro
cat dl_list_mardb_mags.txt | xargs -P 40 wget -P marine_micro/library -nv
cat dl_list_mardb_no_mags.txt | xargs -P 40 wget -P marine_micro/library -nv
cat dl_list_marref.txt | xargs -P 40 wget -P marine_micro/library -nv

# Then we download the NCBI's Taxonomy database
 kraken2-build --download-taxonomy -use-ftp  --db $RES
 
# Download Libraries
 kraken2-build --download-library plasmid -use-ftp  --db  $KDB
 kraken2-build --download-library viral -use-ftp  --db  $KDB
 kraken2-build --download-library protozoa -use-ftp  --db  $KDB
 kraken2-build --download-library plant -use-ftp  --db  $KDB
 kraken2-build --download-library UniVec -use-ftp  --db  $KDB
 kraken2-build --download-library bacteria -use-ftp  --db  $KDB 
 kraken2-build --download-library archaea -use-ftp  --db  $KDB 
 kraken2-build --download-library fungi -use-ftp  --db  $KDB 
 kraken2-build --download-library human -use-ftp --db $KDB
 kraken2-build --download-library nt -use-ftp  --db  $KDB
 kraken2-build --download-library nr -use-ftp --protein  --db  $KDB 

# Add each genome to kraken kmer reference using kraken-build 
find mar_db/marine_micro  -name '*.fa' -print0 | \
        xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db $KDB
        
# Build Kraken Database & Inspect
kraken2-build --build --threads=48 --db  $KDB
kraken2-inspect  --db  $KDB | head -20

# Running Trimmed Reads into Kraken
kraken2  --db  $KDB --report $RES/bellairs.kreport --paired $SEQDATA/Bellairs_Reads/trimmed1.BS.p.fastq.gz  $SEQDATA/Bellairs_Reads/trimmed2.BS.p.fastq.gz  >   $RES/bell.kraken 
 
kraken2  --db  $KDB --report $RES/maycocks.kreport  --paired $SEQDATA/Maycocks_Reads/trimmed1.CS.p.fastq.gz    $SEQDATA/Maycocks_Reads/trimmed2.CS.p.fastq.gz > $RES/maycocks.kraken 
 
 
# Build Bracken Database 
bracken-build -d  $KDB -t 48 -k 35 -l 100 
#-x $KRAKEN
 
# Running Kraken Reports in Bracken 
bracken -d  $KDB -i  $RES/maycocks.kreport -o  $RES/maycocks.bracken -r 100 -l 'S' -t 10
bracken -d  $KDB -i  $RES/bellairs.kreport -o  $RES/bellairs.bracken -r 100 -l 'S' -t 10


#################################################################

# Running Trimgalore Reads into Kraken

kraken2 --db $KDB --report $RES/temp_trimgl_krbr/bellairs_trimgl.kreport --paired $SEQDATA/Bellairs_trimgalore/NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R1_val_1.fq.gz $SEQDATA/Bellairs_trimgalore/NS.A00266_0049.002.IDT_i7_49---IDT_i5_49.BS_R2_val_2.fq.gz > $RES/temp_trimgl_krbr/bell_trimgl.kraken

kraken2 --db $KDB --report $RES/temp_trimgl_krbr/maycocks_trimgl.kreport --paired $SEQDATA/Maycocks_trimgalore/NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R1_val_1.fq.gz $SEQDATA/Maycocks_trimgalore/NS.A00266_0049.002.IDT_i7_61---IDT_i5_61.CS_R2_val_2.fq.gz > $RES/temp_trimgl_krbr/may_trimgl.kraken

# Running Trimgalore kraken reports into Bracken 

bracken -d  $KDB -i  $RES/temp_trimgl_krbr/bellairs_trimgl.kreport -o  $RES/temp_trimgl_krbr/bellairs_trimgl.bracken -r 100 -l 'S' -t 10

bracken -d  $KDB -i  $RES/temp_trimgl_krbr/maycocks_trimgl.kreport -o  $RES/temp_trimgl_krbr/maycocks_trimgl.bracken -r 100 -l 'S' -t 10
