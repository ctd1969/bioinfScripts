#THIS WRITES FASTAS FOR EVERY RAD LOCUS FROM POPULATIONS .FA FILE
#sys.args = [vcf_in] [stacksFasta_in] [path2newFolder - e.g. "fastas"]
import sys, os

def main():
    taxa=taxList(sys.argv[1])
    try:
        os.mkdir(sys.argv[3])
    except:
        print('\n*** FOLDER "'+sys.argv[3]+'" ALREADY EXISTS!!! ***')
        return
    locs=[] #locus names
    haps=[] #grouped loci
    seqs=[]
    cnts=[]
    cnt=0
    prev='>CLocus_1' #dummy variable to allow first id of when we move to next locus
    for line in open(sys.argv[2]): #go thru fasta file, store data and list of all loci
        if line.startswith('>'):
            locs.append(line)
            haps.append(line.split('_')[0]+'_'+line.split('_')[1]) #list all loci
            cnts.append(cnt)
            cnt+=1
        else:
            seqs.append(line)

    cnt=0
    taxCnt=0
    subseqs=[]
    sublocs=[]
    tmpTaxa=[]
    #process all loci & id which taxa are missing
    for hap in haps:
        if locs[cnt].split(' ')[1][1:-2] not in tmpTaxa and hap==prev:
            tmpTaxa.append(locs[cnt].split(' ')[1][1:-2])
            subseqs.append(seqs[cnt])
            sublocs.append(locs[cnt])
            taxCnt+=1
        if hap!=prev and taxCnt> int(len(taxa)/2) - 1: #min 50% taxa per locus
            print(taxCnt)
            snps=[]
            for subseq in subseqs:
                cnt2=0
                for bp in subseq[:-1]:
                    if bp!=subseqs[0][cnt2]:
                        snps.append(cnt2)
                    cnt2+=1
                    cnt3=len(subseq[:-1])
            snpset=set(snps)
            if len(snpset)>0: #Min. no. SNPS per locus
                g=open(sys.argv[3]+'/'+prev[1:]+'.fas','w')
                for tax in taxa:
                    if tax not in tmpTaxa:
                        g.write('>'+tax+'\n')
                        g.write('-'*cnt3)
                        g.write('\n')
                    else:
                        g.write('>'+tax+'\n'+subseqs[tmpTaxa.index(tax)])
                g.close()
            #reset variables
            tmpTaxa=[]
            subseqs=[]
            sublocs=[]
            tmpTaxa.append(locs[cnt].split(' ')[1][1:-2])
            subseqs.append(seqs[cnt])
            sublocs.append(locs[cnt])
            taxCnt=0
        if hap!=prev and taxCnt < int(len(taxa)/2): #reset variables for loci with too few taxa
            tmpTaxa=[]
            subseqs=[]
            sublocs=[]
            tmpTaxa.append(locs[cnt].split(' ')[1][1:-2])
            subseqs.append(seqs[cnt])
            sublocs.append(locs[cnt])
            taxCnt=0
        cnt+=1
        prev=hap

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

if __name__ == '__main__':
    main()
