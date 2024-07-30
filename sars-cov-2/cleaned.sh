DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
tsv_file="$DIR/dataset.tsv"
ref_file=${DIR}/ref/NC_045512.2
# Check if the file exists
if [ ! -f "$tsv_file" ]; then
  echo "Error: File '$tsv_file' not found."
  exit 1
fi
# bowtie2-build "$ref_file.fa" ${DIR}/ref/NC_045512.2
# Loop through each line of the TSV file (excluding header)
while IFS=$'\t' read -r GISAID_assession lineage srArun_acc other_columns; do
    if [[ "$srArun_acc" == "SRArun_acc" || "$srArun_acc" == 'NA' || "$srArun_acc" != 'ERR5405022' ]]; then
        continue
    fi
    echo $srArun_acc
    fastp --in1 ${DIR}/reads/${srArun_acc}_1.fastq.gz --in2 ${DIR}/reads/${srArun_acc}_2.fastq.gz --out1 ${DIR}/trimmed/${srArun_acc}_1.fastq.gz --out2 ${DIR}/trimmed/${srArun_acc}_2.fastq.gz -h ${DIR}/${srArun_acc}_report.html -j ${DIR}/${srArun_acc}_report.json
    (bowtie2 -x ${ref_file} -1 ${DIR}/trimmed/${srArun_acc}_1.fastq.gz -2 ${DIR}/trimmed/${srArun_acc}_2.fastq.gz -S ${DIR}/sam/${srArun_acc}.sam --al-conc-gz ${DIR}/aligned/${srArun_acc}_%.fastq.gz) 2>${DIR}/${srArun_acc}_stat.log 

done < "$tsv_file"
