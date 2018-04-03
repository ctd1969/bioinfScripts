#sys.args = [vcf_in] [csv_out]
import sys, random
print('\nconverts .vcf into CSV w/ 0s, 1s, NAs')

def main():
    swch=0
    swch2=0
    g=open(sys.argv[2],'w')
    #write taxa list
    taxa=taxList(sys.argv[1])
    g.write('locus')
    for tax in taxa:
        g.write(','+tax)
    g.write('\n')
    for line in open(sys.argv[1]):
        z=line.split('\t')
        if line=='\n':
            continue
        #ignore metadata
        if line.startswith('#CHROM'):
            swch=1
            continue
        #only use single SNP loci
        if swch==1 and len(z[4])==1:
            loc=z[0]+'_'+z[1]+'_'+z[3]+'_'+z[4]+','
            g.write(loc)
            snps=z[9:] #ignore non-SNP data
            #cycle thru SNPs
            for snp in snps:
                if snp[-1]=='\n': #watch for line ending
                    swch2=1
                snp=snp[:3] #ignore metadata at each SNP
                #convert to 0s & 1s
                snp=snp.replace('./.','NA')
                snp=snp.replace('0/0','0')
                snp=snp.replace('1/1','1')
                snp=snp.replace('0/1',repr(random.randint(0,1)))
                snp=snp.replace('1/0',repr(random.randint(0,1)))
                if swch2==0:
                    g.write(snp+',')
                elif swch2==1:
                    g.write(snp+'\n')
                    swch2=0
    g.close()

#extract taxon list from VCF file
def taxList(vcfFile):
    taxa=[]
    for line in open(vcfFile):
        z=line.strip().split('\t')
        if line.startswith('#CHROM'):
            for tax in z:
                if z.index(tax)>8:
                    taxa.append(tax.strip()) #list taxon name
            break
    return taxa

if __name__ == '__main__':
    main()

