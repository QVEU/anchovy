'''
consensus.py
Patrick T. Dolan, Ph.D.
Unit Chief Quantitative Virology and Evolution Unit

Usage: python consensus.py fastq_pass/
'''

## Imports
import sys
import pandas as pd
import numpy as np
import os


## read directory
pileupDir=sys.argv[1]#fastq_pass directory

## read barcodes in from input fastq directory
barcodes = [f for f in os.listdir(pileupDir) if os.path.isfile(pileupDir+"/"+f) is False]

print(barcodes)

for bc in barcodes:#for BC directory
    print(bc)
    outfile=pileupDir+"/"+bc+"_cons.fa"
    print(outfile)
    with open(outfile,"w") as OUTFA:
        pass
    files=[i for i in os.listdir(pileupDir+bc) if i.endswith(".pile") and not i.startswith("._")]
    for f in files:
        pileupFile=pileupDir+bc+"/"+f
        if os.path.getsize(pileupFile)>0:
            pileup=pd.read_table(pileupFile,header=None)
            print(pileup)
            allPiles=[]
            consensus=[]
            mod=False
            for i in pileup[4]:
                #print(i)
                pile=""
                n=""
                for c in i:
                    c=c.upper()
                    if c not in ["^","]","$"]:
                        #print(c)
                        if c in ["+","-"]: #new mod
                            mod=True
                            n=""
                        elif mod!=True and (c in ["A","C","T","G","*"]):
                            pile=pile+c #Add base to pileup at each site.
                        elif mod==True and (c in ["A","C","T","G","N"]):#inserted or deleted base
                            if n != "":
                                #print("n: "+str(n))
                                n=int(n)-1  #count down.
                                #print(str(n))
                                if n==0:    #if countdown ended, mod is off.
                                    mod=False
                                    n=""    #n is reset
                        elif mod==True and c!="*":     #number gets added to n
                            n=str(n)+c
                            #print("n: "+str(n))

                allPiles=allPiles+[pile]
                #print(pile)
                print(pile)
                pileCounts=np.unique([p for p in pile],return_counts=True)
                print(pileCounts)
                print(pileCounts[0][np.argmax(pileCounts[1])])
                consensus=consensus+[pileCounts[0][np.argmax(pileCounts[1])]]
    pileup["clNT"]=allPiles
    pileup["cons"]=consensus
    consensusString=""

    for i in range(max(pileup[1])):
    #    print(i)
        if i in pileup[1].values:
            #print(pileup["cons"][pileup[1]==i])
            consensusString=consensusString+pileup["cons"].loc[pileup[1]==i].values
        else:
            consensusString=consensusString+"-"
    with open(outfile,"a") as OUTFA:
        print(outfile)
        OUTFA.write(">"+str(pileupFile.replace("_sort_pile.pile",""))+"\n"+str(consensusString))

#print(consensusString)

#print(pileup)

#print([i for i in pileup.itertuples()])
