#CBCtoFasta.py
#Patrick T. Dolan, Ph.D.
#Unit Chief Quantitative Virology and Evolution Unit
#Usage: python CBCtoFasta.py anchovyOutfile.csv
#
import pandas as pd
import numpy as np
import sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align

def produceFastas(anchout, inputDir):
    print("Converting Anchovy CBCs to fastas for mapping and assembly...")
    for j in np.unique(anchout.CBC):
        if len(anchout.CBC[anchout.CBC==j])>=25:#CHANGED COVERAGE 1/13/23 PD
            cellUMIs=np.unique(anchout.UMI[anchout.CBC==j])
            for i in cellUMIs:
                print(i)
                sequences=anchout.mappedSeq[(anchout.CBC==j)&(anchout.UMI==i)]
                seqObjList=[Seq(S) for S in sequences]
                Align.
            with open(inputDir+"/"+j.strip()+".fa",'w') as OF:
                OF.write("\n".join(["".join([">",str(i[9]),"\n",str(i[10])]) for i in anchout[anchout.CBC==j].itertuples()]))
    print("done.")

if __name__=="__main__":
    inputDir=sys.argv[1]
    FI=inputDir+"/merge_anchovy_v2.csv"
    anchOut=pd.read_csv(FI)
    produceFastas(anchOut, inputDir)
