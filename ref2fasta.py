import numpy as np
import linecache

ref_gen = 'GCA_004354405.1_ASM435440v1_genomic.fna' #From NCBI
chrom_no = 3 #i.e. Chromosome 3
bps = [20000, 22000] #between bases 20000-22000

def main(chrom_no):
    print(f"Chrom: {chrom_no}:{bps}")
    chrom_no -= 1
    prev_lines, linelens2 = refSumm(ref_gen, chrom_no)
    sequence = seq_out(prev_lines, linelens2, chrom_no, bps)
    print(sequence)

#summarise reference genome
def refSumm(ref_gen, chrom_no):
    linelens, tot, cnt = [], 0, 0
    for line in open(ref_gen):
        if cnt > chrom_no + 1: break
        if line.startswith('>'):
            linelens.append([])
            tot = 0
            cnt += 1
        else:
            linelens[-1].append(tot + len(line.strip()))
            tot += len(line.strip())
        
    prev_lines = 1
    for l in linelens[:chrom_no]: prev_lines += len(l) + 1 #no. of lines in ref preceeding this chr

    return prev_lines, np.array(linelens[chrom_no]) #line end counts for ea line in this chr

#chop out sequence
def seq_out(prev_lines, linelens2, chrom_no, bps):
    start, end = bps[0] - 1, bps[-1] #75bp before & after 1st & last SNP
    seq_lines1 = np.where(linelens2 > start)[0][0] #get line # of seq data containing start posn
    seq_lines2 = np.where(linelens2 > end)[0][0] #get line # of seq data containing end posn
        
    #work out vals for trimming excess from seqs
    trim1 =  start - linelens2[seq_lines1 - 1]
    trim2 = -1 * (linelens2[seq_lines2] - end)
        
    #get the seq data
    start_seq = seq_lines1 + prev_lines + 1 #MAY GET WRONG LINES HERE!!!
    end_seq = seq_lines2 + prev_lines + 1
    seq_dat = ''
    for txt in range(start_seq, end_seq + 1):
        seq_dat += linecache.getline(ref_gen, txt).strip()

    seq_dat = seq_dat[trim1:] #trim to gap value either side of 1st/last SNP
    seq_dat = seq_dat[:trim2]

    return seq_dat

if __name__ == '__main__':
    main(chrom_no)



