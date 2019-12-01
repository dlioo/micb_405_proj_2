cd ~
if [ ! -d project_2 ]; then
  mkdir project_2
fi
cd ~/project_2
if [ ! -d transcriptome ]; then
  mkdir transcriptome
fi

fastq=/projects/micb405/resources/project_2/2019/SaanichInlet_100m/MetaT_reads/7705.2.81665.GTGAAA.qtrim.3ptrim.artifact.rRNA.clean.fastq.gz

mag="bacteria"

bwa index ~/project_2/ffn/total_bacteria.ffn

  ffn=~/project_2/ffn/total_bacteria.ffn

    cd ~/project_2/transcriptome
    if [ ! -d $mag ]; then
      mkdir $mag
    fi

    bwa mem -t 8 -p $ffn $fastq > ~/project_2/transcriptome/${mag}/${mag}_MAG_ORFs.sam

    /projects/micb405/resources/project_2/2019/rpkm \
    -c $ffn \
    -a ~/project_2/transcriptome/${mag}/${mag}_MAG_ORFs.sam \
    -o ~/project_2/transcriptome/${mag}/${mag}_MAG_ORFs_RPKM.csv \
    > ~/project_2/transcriptome/${mag}/${mag}_MAG_ORFs_RPKM.log

    rm  ~/project_2/transcriptome/${mag}/${mag}_MAG_ORFs.sam
