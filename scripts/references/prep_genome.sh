#!/usr/bin/bash
# Download and process reference genome files and aligner indices
# Software Requirements:
#   - samtools
#   - bowtie2-build
#   - STAR
# Example Usage: bash prep_genome.sh hg38 mm10 PhiX 

set -euo pipefail
genomes=$@
###############################################
#################### fasta ####################
###############################################
for genome in $@
do
    echo "Downloading $genome"
    
    # download_UCSC_chromosomes
    if [ "$genome" = "PhiX" ]
    then
        dir_tar=$genome
        dir_unzipped=$genome/IlluminaDownload
        mkdir -p $dir_tar $dir_unzipped
        wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz -O $dir_tar/PhiX_Illumina_RTA.tar.gz
        tar -xzf $dir_tar/PhiX_Illumina_RTA.tar.gz -C $dir_unzipped
        cp $dir_unzipped/PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa $genome/$genome.fa
    else
        dir_chrom=$genome/UCSC-chromosomes
        mkdir -p $dir_chrom
        rsync -a rsync://hgdownload.cse.ucsc.edu/goldenPath/$genome/chromosomes/ $dir_chrom

        # filter out haplotype or variant chromosomes, then concatenate into a full fasta 
        ls $dir_chrom/*.fa.gz | grep -v _hap | grep -v _alt | xargs zcat > $genome/$genome.fa
    fi
    
    # index fastas
    samtools faidx $genome/$genome.fa
    
done

###############################################
################### bowtie2 ###################
###############################################
for genome in $@
do
    mkdir -p $genome/bowtie2
    bowtie2-build -f --threads 16 $genome/$genome.fa $genome/bowtie2/$genome
done

###############################################
##################### star ####################
###############################################
for genome in $@
do
    if [ $genome = "hg38" ] # TODO support other genomes
    then
        gtf_link=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz
        gtf_out=gencode.v41.annotation.filtered.gtf
        Rscript prep_gtf_reference.R $gtf_link $genome/$gtf_out

        mkdir -p $genome/star/
        STAR --runMode genomeGenerate --runThreadN 16 --genomeDir $genome/star/ --genomeFastaFiles $genome/$genome.fa --sjdbGTFfile $genome/$gtf_out --sjdbOverhang 99
    fi
done

