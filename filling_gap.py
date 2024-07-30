from pysam import VariantFile
from os import listdir
from os.path import isfile, join
# read in original VCF file

onlyfiles = [f for f in listdir('./vcf') if isfile(join('./vcf', f)) if 'filtered' not in f and 'csi' not in f]
for file in onlyfiles:
    print(file)
    vcf_in = VariantFile(f"./vcf/{file}")

    # initializing output file
    vcf_out = VariantFile(f'./vcf/filled_{file}', 'w', header=vcf_in.header)

    for rec in vcf_in.fetch():
        # adding omni tag to record
        if len(rec.alleles[0])<len(rec.alleles[1]):
            continue
        if len(rec.alleles[0])>len(rec.alleles[1]):
            new_rec =  rec.alleles[1]+'N'* (len(rec.alleles[0])-len(rec.alleles[1]))
            rec.alleles = tuple([rec.alleles[0],new_rec])
        vcf_out.write(rec)
    vcf_out.close()   