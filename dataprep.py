import pandas as pd
df=pd.read_csv('./salmonella/dataset.tsv',delimiter='\t', index_col=False)
# print(df['SRArun_acc'])
print(df.keys())
df['ftpRead1']=df['SRArun_acc'].apply(lambda srr: f'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{srr[:6]}/0{srr[-2:]}/{srr}/{srr}_1.fastq.gz')
df['ftpRead2']=df['SRArun_acc'].apply(lambda srr: f'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{srr[:6]}/0{srr[-2:]}/{srr}/{srr}_2.fastq.gz')
df.to_csv('./salmonella/cleaned_dataset.tsv',sep='\t', index=False)
