#sys.argv = [folder_containing_fastas] [outFile_fas]
#This pgm concatenates FASTA files (that have identical taxa e.g. from Stacks loci)
#NB very old script: maybe odd syntax - but works

import sys, os
fileList=os.listdir(sys.argv[1])
os.chdir(sys.argv[1])

taxa=[]
y=0
z=1
for line in open(fileList[0]):
    if line.startswith('>'):
        taxa.append(line)

f=open(sys.argv[2],'w')
for name in taxa: #cycle thru and write each taxon name
    x=0
    seq=""
    f.write(name)
    while x<len(fileList): #cycle thru all fastas
        y=0
        g=open(fileList[x])
        while y<z:
            g.readline()
            newSeq=g.readline()
            y+=1

        seq=seq+newSeq[0:-1]
        x+=1

    z+=1
    f.write(seq+'\n')
    
print("************ HUNKY DORY! ***************")
f.close()
