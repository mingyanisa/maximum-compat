from os import listdir
from os.path import isfile, join
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

onlyfiles = [f for f in listdir('./sars-cov-2/voivoc/snp-consensus') if isfile(join('./sars-cov-2/voivoc/snp-consensus', f))]
min_len=None
for file in onlyfiles:
    print(file)
    for rec in SeqIO.parse(f'./sars-cov-2/voivoc/consensus/{file}', 'fasta'):
        name = rec.id
        seq = rec.seq
        seqLen = len(rec)
        print(f"{name}, {seqLen}")
        # if 'CP006053.1' in name:
        min_len=len(rec.seq) if min_len==None or min_len>len(rec.seq) else min_len
print(min_len)

#salmonella: ENA|CP006053|CP006053.1, 4730612
#voivoc: NC_045512.2|ERR5181310, 29850
#non-voivoc: NC_045512.2|SRR13195265, 29897

with open('./sars-cov-2/voivoc/combine_snp_sarscov2.fa', 'r') as f_in, open('./sars-cov-2/voivoc/trimmed_snp_sarscov2.fa', 'w') as f_out:
    for record in SeqIO.parse(f_in, "fasta"):
        name=record.id
        print(f"name:{name}\n")
        # if 'CP006053.1' in name:
        seq=Seq(record.seq[:min_len])
        # else: 
        #     seq=record.seq
        record = SeqRecord(
            seq,
            id=record.id,
            name=record.name,
            description=record.description,
        )
        SeqIO.write(record, f_out, "fasta")  # Write trimmed record

for rec in SeqIO.parse(f'./sars-cov-2/voivoc/trimmed_snp_sarscov2.fa', 'fasta'):
    name = rec.id
    seq = rec.seq
    seqLen = len(rec)
    print(f"{name}, {seqLen}")
