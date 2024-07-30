#!/bin/bash
# for filename in ./data/*.fasta; do
#     echo $filename
#     minimap2 -a ./data/GCA_000439415.1.fasta $filename > $(basename "$filename" .fasta).sam
# done
# for filename in ./*.bam; do
#     echo $filename
#     # samtools view -bS $filename > $(basename "$filename" .sam).bam
#     samtools sort  $filename -o $(basename "$filename" .bam).sorted.bam
# done
for filename in ./sorted/*.sorted.bam; do
    echo $filename
    # samtools index $filename $(basename "$filename" .sorted.bam).bai
    # samtools bam2fq $filename| seqtk seq -A > ./fasta/$(basename "$filename" .sorted.bam).fasta
    # bcftools mpileup -Ou -f ./data/GCA_000439415.1.fasta $filename | bcftools call -mv -Oz -o ./vcf/$(basename "$filename" .sorted.bam).vcf.gz
    # bcftools view --exclude-types indels ./vcf/$(basename "$filename" .sorted.bam).vcf.gz > ./vcf/filtered_$(basename "$filename" .sorted.bam).vcf
    bgzip ./vcf/$(basename "$filename" .sorted.bam).vcf
    bcftools index ./vcf/$(basename "$filename" .sorted.bam).vcf.gz
    cat ./data/GCA_000439415.1.fasta | bcftools consensus ./vcf/$(basename "$filename" .sorted.bam).vcf.gz > ./consensus/$(basename "$filename" .sorted.bam).fa
done
# cat ./snp_consensus/*.fa > ./combine_salmonella.fa
# compat ./trimmed_salmonella.fa ./trimmed_snp_salmonella_tree.netwick