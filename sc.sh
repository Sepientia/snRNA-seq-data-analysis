sample=$1
mv ${sample}_f1.fq.gz ${sample}_S1_R1_001.fastq.gz
mv ${sample}_r2.fq.gz ${sample}_S1_R2_001.fastq.gz
/D5/NGS/software/cellranger-7.1.0/bin/cellranger count --localcores=4 --nosecondary --disable-ui --id=${sample} --fastqs=./ --sample=${sample} --transcriptome=/D5/NGS/software/refdata-gex-mm10-2020-A
