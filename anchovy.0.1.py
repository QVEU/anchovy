#!/usr/bin/env python
'''
anchovy.0.1.py
  ><>  ><>
    ><>
Patrick T. Dolan
Unit Chief, Quantitative Virology and Evolution Unit
7/18/22
USAGE: python anchovy.0.1.py pathto/input.sam query
Identify barcodes from nanopore sequencing to map back onto cells .
v0.1: Identify ligation products in the sequencing.
'''

##### Imports #####
import pandas as pd
import Levenshtein
import time
import re
import numpy as np
import sys
from multiprocessing import Pool
import pysam

##### Functions #####
def Initialize(args):
    '''
    Function: Initialize()
        Reads in and filters data and initializes a bunch of parameters for the mapping.

    Arguments:
        "seq_query" is a *TUPLE* of the template sequence and the query. Need to be sent as a package for multi-threading with Pool.
    Note: replaces position if alternative minimum sites are found. Producing a 3' bias? Need to consider improvements here.
    '''
    #Input SAM
    print("\n\n----------------=============-----------------")
    print("==--==--==--==--    ><>      --==--==--==--==-")
    print("--==--==--==--==  ><>   ><>  ==--==--==--==--=")
    print("==--==--==--==--   anchovy   --==--==--==--==-")
    print("----------------=============-----------------\n\n")


    #Check for Args
    if len(args)<3:
        print("USAGE: python anchovy.0.1.py pathto/input.sam pathto/whitelist.txt query")
        exit()

    #Query
    query=args[3]
    #5' query: CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNTTTCTTATAT
    whitelist=args[2]

    quL=len(query)
    print("Query Length: {}".format(quL))

    #input SAM
    inputData=sys.argv[1]
    print("\nLoading SAM: {}".format(inputData))
    pdSam=loadSAM(inputData,quL)
    pdCBCs=loadBC(whitelist)
    outfile=inputData.replace(".sam",'_anchovy_v1.csv')
    print("\nOutfile: {}".format(outfile))

    return(str(query), pdCBCs, pdSam, outfile)

def loadSAM(inputData,quL):
    '''
    Function: loadSAM()
        Reads in SAM file.
        Since this is the slow step, I broke it out to work on later.

    Arguments:
        "inputData" is a string of the path to the input SAM file.
    Note: quite slow, how can we speed up? pd.read seems to have issues with the headers on the file.
    '''
    with open(inputData,"r") as IF:
        samInput=[lin.split("\t")[0:11] for lin in IF if lin.startswith("@") is False]#Filter out headers with @

    names=["read","FLAG","template","pos","mapq","cigar","Rnext","Pnext","Tlen","seq","Qscore"]#SAM columns
    samFile    = sys.argv[1];
    print("Parsing Cigars...")
    samFP = pysam.Samfile(samFile, "rb");
    cigars=[]
    for read in samFP:
        readLen=0#match
        CreadLen=0#match
        if( not( read.is_unmapped ) ):   #if it's mapped
            cigarLine=read.cigar;
            switch=0
            for (cigarType,cigarLength) in cigarLine:
                try:
                    if(cigarType == 0):
                        readLen+=cigarLength#match
                        CreadLen+=cigarLength#match
                        switch=1
                    elif(cigarType == 1):
                        readLen+=cigarLength#insertions
                        CreadLen+=cigarLength#insertions
                        switch=1
                    elif(cigarType == 2): pass #deletion
                    elif(cigarType == 3):
                        readLen+=cigarLength#skip
                        CreadLen+=cigarLength#skip
                        switch=1
                    elif(cigarType == 4):
                        if switch==0:#only count opening soft clipping
                            readLen+=cigarLength#soft clipping
                        switch=1
                    elif(cigarType == 5): pass#hard clipping
                    elif(cigarType == 6): pass#padding
                    else:
                        print ("Wrong CIGAR number");
                        sys.exit(1);
                except:
                    print("Problem")
            cigarLengths=[readLen,CreadLen,readLen-CreadLen]
            cigars.append(cigarLengths)

    print("Done.")
    minL=quL
    size=len(samInput)
    print("Total Candidate Reads: {}".format(size))
    pdSam=pd.DataFrame(samInput,columns=names)
    pdSam=pdSam[pdSam.template!="*"]
    print(len(pdSam))
    print(len(cigars))
    pdSam[['readLen','clipReadLen','offset']]=cigars
    pdSam['length']=pdSam.seq.apply(lambda s: len(s))
    pdSam=pdSam[(pdSam.length>minL)&(pdSam.template!="*")]
    return(pdSam)

def loadBC(whitelist):
    '''
    Function: loadBC()
        Reads in 10X whitelist.

    Arguments:
        "whitelist" is a string of the path to the input 10x whitelist of cell barcodes.
    '''

    with open(whitelist,"r") as IF:
        CBCs=[lin.split("\t")[0] for lin in IF]
    size=len(CBCs)
    print("Total Cell Barcodes: {}".format(size))

    pdCBC=pd.DataFrame(CBCs,columns=["CBC"])

    return(pdCBC)

def blockDist(seq_query):
    '''
    Function: blockDist()
        unpacks the sequence query and template, makes chunks of query length from template and then computes the distance for each kmer.
        Uses numpy vectorize to speed up the computation (formerly for-loop)

    Arguments:
        `seq_query` is a tuple of read and query string, passed this way for 'pool'
    '''
    #Unpack Values
    seq=seq_query[0]
    query=seq_query[1]

    #Prepare Query
    quL=len(query)
    minD=quL #set maximum distance to length of query
    minPos=0
    matchseq=""
    nN=np.sum([1 for i in query if i=="N"])
    # Build vectorized BlockDistance function
    BlockDistance = lambda Block: Levenshtein.distance(Block,query)
    vectorBlockDist = np.vectorize(BlockDistance,otypes=[int])

    # Make blocks of read sequence
    blocks=np.array([seq[i:(i+quL)] for i in range(len(seq)-quL)])
    bDist = vectorBlockDist(blocks)
    #hitPosList=np.argsort(bDist)[np.sort(bDist)<=(1.20*(nN))]
    #print(hitPosList)
    #compute minimum and position of minimum
    try:
        minD=min(bDist.astype(int))
        minPos=np.argmin(bDist.astype(int))
        matchseq=blocks[minPos] #Grab this matched sequence.
    except:
        print("No Hit:"+str(seq_query))

    return(minD, minPos, matchseq)

def poolBlocks(query,pdSam,nthreads=16):#uses multithreading to compute the matches
    with Pool(nthreads) as p:
        print("\n1. Computing minimum distance hit position for {} reads.".format(len(pdSam)))
        quL=len(query)
        print("Query: {}".format(query))
        print( [c[10][(c[14]-100):(c[14]+5)].upper() for c in pdSam.itertuples()] )#ONLY MAP WITHIN 100 bases upstream of the hit site. 
        pdSam['minD'],   pdSam['minPos'],   pdSam['matchseq'] = zip(*p.map(blockDist, [(c[10][(c[14]-100):(c[14]+5)].upper(),query.upper()) for c in pdSam.itertuples()]))
    return(pdSam)

def cellMatch(input): #TUPLE including (readID, readSeq, matchSeq, matchseq, CBCs, blocks):
    '''
    Function: cellMatch()

    Arguments:
        `input` is tuple containing the following:
            readID=input[0]
            readSeq=input[1]
            matchPos=input[2]
            matchseq=input[3]
            CBCs=input[4]
            blocks=input[5]
            #f[1], f[9], f[16], f[17], pdCBCs, blocks
    '''
    readID=input[0]
    readSeq=input[1]
    matchPos=input[2]
    print("Pos of Index Match:"+str(matchPos))
    matchseq=input[3]
    print("Index Match Seq:"+str(matchseq))
    offset=input[4]
    print("offset:"+str(offset))
    CBCs=input[5]
    blocks=input[6]

    #Define Levenshtein distance function
    BlockDistance = lambda Block: Levenshtein.distance(Block,matchseq)
    #Vectorize the BlockDistance function
    vectorBlockDist = np.vectorize(BlockDistance,otypes=[int])

    #Generate Block Distance vector, 'bDist'
    bDist = vectorBlockDist(blocks)

    #compute minimum and position of minimum
    minD=min(bDist.astype(int))
    minPos=np.argmin(bDist.astype(int))
    matchblock=blocks[minPos]

    return(CBCs.iloc[minPos][0], minD, minPos, matchPos, offset, matchblock, matchseq, readID, readSeq)

#Pooled cell ID function
def cellIDPool(pdSam, pdCBCs, nthreads=16):
    blocks=np.array(["CTACACGACGCTCTTCCGATCT"+i+"NNNNNNNNNNTTTCTTATAT" for i in pdCBCs.CBC])
    pdSam=pdSam[pdSam.minD<33]#filter by distance of index sequences cassette to template
    anchOut=pd.DataFrame()
    print(list(pdSam))
    print(pdSam)
    with Pool(nthreads) as p:
        print("\n2. Identifying Cell Barcodes...")
        anchOut['CBC'], anchOut['minD'], anchOut['BC_ID'], anchOut['readPos'], anchOut['cigarClip'], anchOut['reconQuery'], anchOut['matchseq'], anchOut['read'], anchOut['readSeq'] = zip(*p.map(cellMatch, [(f[1], f[10], f[17], f[18], f[14], pdCBCs, blocks) for f in pdSam.itertuples()]))
    print(anchOut)
    return(anchOut)

if __name__=="__main__":
    t0=(time.time())
    query, pdCBCs, pdSam, outfileName = Initialize(sys.argv)
    pdSam = poolBlocks(query,pdSam)
    print("\nMapped hits in {} reads.".format(np.sum([int(i)>0 for i in pdSam.minPos])))
    anchOut = cellIDPool(pdSam, pdCBCs)
    anchOut.to_csv(outfileName)
    t1=(time.time())
    print("Done in {} minutes.".format((t1-t0)/60))
    print("Wrote {}.".format(outfileName))
