
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# Define the TSV file path (modify as needed)
tsv_file="$DIR/cleaned_dataset.tsv"

# Check if the file exists
if [ ! -f "$tsv_file" ]; then
  echo "Error: File '$tsv_file' not found."
  exit 1
fi

# Loop through each line of the TSV file (excluding header)
while IFS=$'\t' read -r GISAID_assession lineage srArun_acc biosample_acc ftpRead1 ftpRead2 other_columns; do
    # Skip header line (assuming first line is header)
    echo "$srArun_acc" 
    if [[ "$srArun_acc" == "SRArun_acc" || "$srArun_acc" == 'NA' ]]; then
        continue
    fi
    if [[ -f "$DIR/reads/${srArun_acc}_1.fastq.gz" || -f "$DIR/reads/${srArun_acc}_2.fastq.gz" ]]; then
        echo "File '$srArun_acc' already exists, skipping..."
    else
        # Print the extracted SRArun_acc
       
        # echo "$ftpRead1" "$ftpRead2"
        wget -t 0 -O $DIR/reads/${srArun_acc}_1.fastq.gz ${ftpRead1}
        wget -t 0 -O $DIR/reads/${srArun_acc}_2.fastq.gz ${ftpRead2}
        # Add error handling for wget if needed
    fi

done < "$tsv_file"
# < <(tail -n +2 "$tsv_file")

echo "SRArun_acc extraction complete."