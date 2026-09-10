#ConsensusTool.py

# Usage: 
#   > ConsensusTool.py infile ORFstart ORFend 

from Bio import SeqIO as io
import numpy as np
import pandas as pd 
import os 
import sys

#function: readSeqs():
# Args:
#	infile: consensus fastas from anchovy 

def readSeqs(infile):
    print("Parsing Consensus Sequences...")
    parseFA=io.parse(infile,'fasta')
    depth=[float(k.description.split(" ")[2].split("coverage:")[1]) for k in parseFA]

    #parse and read length
    parseFA=io.parse(infile,'fasta')
    length=[int(k.description.split(" ")[3].split("length:")[1]) for k in parseFA]

    #parse and read all entries
    parseFA=io.parse(infile,'fasta')
    elements=[k for k in parseFA]
    return(elements,length,depth)

def selectSeqs(start,end,elements,length,depth,depthMin=1,lengthMin=1):
    print("Filtering Alignment...")
#    filtered=[s for s in elements if sum([s.seq[i]=="-" for i in range(start,end)])==0]
    filtered=[elements[x] for x in range(len(elements)) if ((depth[x]>depthMin) and (sum([elements[x].seq[i]=="-" for i in range(start,end)])<3))]
    print("Trimming Alignment...")
    trimmedSeqs=[filtered[i].seq[start:end] for i in range(len(filtered))]
    for i in range(len(filtered)):
        filtered[i].seq=trimmedSeqs[i]
    print("Total: "+str(len(filtered)))
    return(filtered)

def genotypeSummary(filtered,reference=""):
    if reference=="":
        clists=[[c for c in s.seq] for s in filtered]#for column in filtered alignment...
        tclists=np.transpose(clists)#transpose columns 
        #print(tclists)
        variantSites=[i for i in range(len(tclists)) if len(np.unique(tclists[i],return_counts=True)[0])>1]# sites with variation
        maxChars=[]
        for site in range(len(tclists)):
            maxChars=maxChars+[np.unique(tclists[site],return_counts=True)[0][np.argmax(np.unique(tclists[site],return_counts=True)[1])]]
        consensusSeq="".join(maxChars)
        print("Consensus Seq: "+consensusSeq[0:30]+"... ..."+consensusSeq[(len(consensusSeq)-30):len(consensusSeq)])
        output=[]    
        for s in filtered:
            s.dbxrefs="_".join([str(i+1)+s.seq[i] for i in variantSites if s.seq[i]!=maxChars[i]])
            #print([s.seq[(i/3-i%3):(i%3+3)] for i in variantSites if s.seq[i]!=maxChars[i]])
            #s.translation="_".join([str(i)+s.seq[i] for i in variantSites if s.seq[i]!=maxChars[i]])
            #print(s.dbxrefs)
            output=output+[s]
        return(output,consensusSeq)
    else:
        clists=[[c for c in s.seq] for s in filtered]#for column in filtered alignment...
        tclists=np.transpose(clists)#transpose columns 
        variantSites=[i for i in range(len(tclists)) if len(np.unique(tclists[i],return_counts=True)[0])>1]# sites with variation
        output=[]  
        for s in filtered:
            s.dbxrefs="_".join([str(i+1)+s.seq[i] for i in variantSites if s.seq[i]!=reference[i]])
            #print([s.seq[(i/3-i%3):(i%3+3)] for i in variantSites if s.seq[i]!=maxChars[i]])
            #s.translation="_".join([str(i)+s.seq[i] for i in variantSites if s.seq[i]!=maxChars[i]])
            #print(s.dbxrefs)
            output=output+[s]
        return (output, consensusSeq)

def main():
    #Arg inputs
    infile=sys.argv[1] #infile="/Volumes/lvd_qve/Projects/DENV_SEARCHLIGHT/Anchovy_2/DENV_6dpi_allConsensus.fasta"
    start=int(sys.argv[2]) #start of ORF or region of interest
    end=int(sys.argv[3]) #end of ORF or region of interest

    # Arg inputs
    elements,length,depth = readSeqs(infile)
    filtered = selectSeqs(start,end,elements,length,depth,depthMin=10,lengthMin=1)
    
    output,consensusSeq=genotypeSummary(filtered)
    
    with open(str.replace(infile,"_allConsensus.fasta","_consensus_reference.txt"), "w") as file:
        file.write(consensusSeq)
        print("... Wrote consensus to file.")
    
    io.write(output,handle=str.replace(infile,"_allConsensus.fasta","_filtConsensus.fasta"),format="fasta")
    print("Writing tables...")
    infoFrame=pd.DataFrame([[outputs.id,outputs.dbxrefs,"".join(outputs.seq),outputs.description] for outputs in output],columns=["CBC_ID","genotype","sequence","description"])
    infoFrame.to_csv(str.replace(infile,"_allConsensus.fasta","_filtConsensus.csv"))
    print("Done.")

if __name__ == "__main__":
    main()