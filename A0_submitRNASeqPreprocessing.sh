#' ---
#' title: "submitRNASeqPreprocessing"
#' author: "Nicolas Delhomme & Thomas Riquelme"

#!/bin/bash -l

# stop on error and be verbose
set -ex 

# environment variables
proj=u2015030

mail="thri0013@gapps.umu.se"

in=/mnt/picea/projects/spruce/sjansson/seasonal-needles/20170808/raw

out=/mnt/picea/projects/spruce/sjansson/seasonal-needles/20170808/

genome=/mnt/picea/storage/reference/Picea-abies/v1.1/indices/STAR/v2.5.2b/Pabies01-genome
#/mnt/picea/storage/reference/Populus-tremula/v1.1/indices/STAR/2.5.2b/Potra01

gff3=/mnt/picea/storage/reference/Picea-abies/v1.1/gff3/Eugene.gff3 
#/mnt/picea/storage/reference/Populus-tremula/v1.1/gff3/Potra01-gene-synthetic-transcripts-wo-intron.gff3

gtf=/mnt/picea/storage/reference/Picea-abies/v1.1/gtf/Eugene.gtf
#/mnt/picea/storage/reference/Populus-tremula/v1.1/gtf/Potra01-gene-mRNA-wo-intron.gtf

kallistoFasta=/mnt/picea/storage/reference/Picea-abies/v1.0/fasta/GenePrediction/phased/Pabies1.0-all.phase.gff3.CDS.fa
#/mnt/picea/storage/reference/Populus-tremula/v1.1/fasta/Potra01-mRNA.fa
			  
kallistoIndex=/mnt/picea/storage/reference/Picea-abies/v1.0/indices/kallisto/Pabies1.0-all.phase.gff3.CDS.fa.inx
#/mnt/picea/storage/reference/Populus-tremula/v1.1/indices/kallisto/Potra01-mRNA.fa.inx

# all steps from 1 to 9 (cf runRNASeqPreprocessing.sh)
start=1
end=9

#RAM
mem=256 # minimum 160 Go because indexed genomes'size is of this range, plus we take a bit larger for buffer

#load necessary modules
module load bioinfo-tools samtools fastQvalidator FastQC sortmerna Trimmomatic star htseq kallisto 

# checks if UPSCb variable is defined in user environment
if [ -z $UPSCb ]; then 
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir." 
    exit 1
fi

# loop to get every PE fastq file and apply the script runRNASeqPreprocessing.sh to it
for f in `find $in -name "*_[1,2].fastq.gz"`; do echo "${f//_[1,2].fastq.gz/}" ; done | sort | uniq | while read line;
do
bash $UPSCb/pipeline/runRNASeqPreprocessing.sh -s $start -e $end -m $mem -g $genome -G $gtf -H $gff3 -I 70000 -f $kallistoFasta -K $kallistoIndex $proj $mail ${line}_1.fastq.gz ${line}_2.fastq.gz $out
done