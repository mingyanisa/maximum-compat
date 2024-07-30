#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# for filename in ./data/*.fasta; do
#     echo $filename
#     minimap2 -a ./data/GCA_000439415.1.fasta $filename > $(basename "$filename" .fasta).sam
# done
# for filename in $DIR/sars-cov-2/sam/*.sam; do
#     echo $filename
#     samtools view -bS $filename > $DIR/sars-cov-2/sam/$(basename "$filename" .sam).bam
#     samtools sort $filename -o $DIR/sars-cov-2/sorted/$(basename "$filename" .sam).sorted.bam
# done
for filename in $DIR/sars-cov-2/voivoc/sorted/*.sorted.bam; do
    echo $filename
    # samtools index $filename $DIR/sars-cov-2/sam/$(basename "$filename" .sorted.bam).bai
    # samtools bam2fq $filename| seqtk seq -A > $DIR/sars-cov-2/sorted-fasta/$(basename "$filename" .sorted.bam).fasta
    # bcftools mpileup -Ou -f $DIR/sars-cov-2/ref/NC_045512.2.fa $filename | bcftools call -mv -Oz -o $DIR/sars-cov-2/vcf/$(basename "$filename" .sorted.bam).vcf.gz
    bcftools view --exclude-types indels $DIR/sars-cov-2/voivoc/vcf/$(basename "$filename" .sorted.bam).vcf.gz > $DIR/sars-cov-2/voivoc/snp-vcf/$(basename "$filename" .sorted.bam).vcf
    bgzip $DIR/sars-cov-2/voivoc/snp-vcf/$(basename "$filename" .sorted.bam).vcf
    tabix $DIR/sars-cov-2/voivoc/snp-vcf/$(basename "$filename" .sorted.bam).vcf.gz
    cat $DIR/sars-cov-2/ref/NC_045512.2.fa | bcftools consensus $DIR/sars-cov-2/voivoc/snp-vcf/$(basename "$filename" .sorted.bam).vcf.gz > $DIR/sars-cov-2/voivoc/snp-consensus/$(basename "$filename" .sorted.bam).fa
done
# cat ./snp_consensus/*.fa > ./combine_salmonella.fa
# compat ./trimmed_salmonella.fa ./trimmed_snp_salmonella_tree.netwick