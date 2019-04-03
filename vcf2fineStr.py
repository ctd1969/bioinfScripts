#converts VCF from Stacks:populations to fineStructureRAD format
#sys.args = [vcf_in] [CSV_master_location_file] [txt_out]
#NB DOES NOT CONSIDER PHASED DATA - see fineRADstructure manual

import sys, pandas as pd

def main():
    g = open(sys.argv[3],'w')
    taxa = taxList(sys.argv[1])
    isls = popper(taxa, sys.argv[2])
    setIsls = set(isls)
    #for taxon list with no root names it will name each taxon with a population name (and integer)
    for i in setIsls:
        cnt, cnt2 = 0, 1
        for j in isls:
            if i == j:
                isls[cnt] = j+repr(cnt2)
                cnt2 += 1
            cnt += 1
    for isl in isls[:-1]: g.write(isl+'\t')
    g.write(isls[-1]+'\n')
    #convert VCF to fineStr format
    for line in open(sys.argv[1]):
        if line[0] != '#':
            als, z = [], line.strip().split('\t')
            als.append(z[3])
            for _ in z[4].split(','): als.append(_) #all alt alleles
            if len(als) > 8: continue #from 9th column onwards
            for t in z[9:-1]:
                if t[:3] == './.': #id missing data
                    g.write('\t')
                    continue
                elif t[0] == t[2]: #or write allele nucleotides
                    g.write(als[int(t[0])] + '\t')
                    continue
                elif t[0] != t[2]:
                    g.write('{}/{}\t'.format(als[int(t[0])], als[int(t[2])]))
                    continue
            if z[-1][:3] == './.': #deal with line endings
                g.write('\t\n')
                continue
            elif z[-1][0] == z[-1][2]:
                g.write('{}\n'.format(als[int(z[-1][0])]))
                continue
            elif z[-1][0] != z[-1][2]:
                g.write(als[int(z[-1][0])] + '/'+als[int(z[-1][2])] + '\n')
                continue
    g.close()

#extract taxon list from VCF file
def taxList(vcfFile):
    taxa=[]
    for line in open(vcfFile):
        z=line.strip().split('\t')
        if line.startswith('#CHROM'):
            [taxa.append(tax.strip()) for tax in z if z.index(tax) > 8]
            break
    return taxa

#determine population
def popper(taxList, datFile):
    isls, file = [None] * len(taxList), pd.read_csv(datFile)
    for tax in taxList:
        if tax in list(file['dna_extraction_code']): # i.e. header = dna_extraction_code
            isls[taxList.index(tax)]=file['location'][list(file['dna_extraction_code']).index(tax)] #appropriate index of location for taxon
    return isls

if __name__ == '__main__':
    main()
