#sys.args = [vcf_in] [CSV_master_location_file] [CSV_out]
import sys, pandas as pd
def main():
    g=open(sys.argv[3],'w')
    taxa=taxList(sys.argv[1])
    data=popper(taxa, sys.argv[2])
    for dat in data:
        g.write(dat+','+data.get(dat)+'\n')
    g.close()

#extract taxon list from VCF file
def taxList(vcfFile):
    taxa=[]
    for line in open(vcfFile):
        z=line.strip().split('\t')
        if line.startswith('#CHROM'):
            for tax in z:
                if z.index(tax)>8: # taxa listed after 8th column
                    taxa.append(tax.strip()) #list taxon name
            break
    return taxa

#associate taxa with population data
def popper(taxList, datFile):
    isls=[None]*len(taxList)
    file = pd.read_csv(datFile) #pandas assumes headers from CSV
    for tax in taxList:
        if tax in list(file['dna_extraction_code']): # i.e. header = dna_extraction_code
            isls[taxList.index(tax)]=file['location'][list(file['dna_extraction_code']).index(tax)] #appropriate index of location for taxon
    dik=dict(zip(taxList,isls))
    return dik

if __name__ == '__main__':
    main()
