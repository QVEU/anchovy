#CellConsensus.py
import pandas as pd
import numpy as np
import sys

def produceFastas(anchout, inputDir):
    print("Mapping cell consensus")
    for j in np.unique(anchout.CBC):
        print(len(anchout.CBC[anchout.CBC==j]))
        print("--------")
        print(j)
        if len(anchout.CBC[anchout.CBC==j])>=5:
            with open(inputDir+"/"+j+".fa",'w') as OF:
                OF.write("\n".join(["".join([">",str(i[9]),"\n",str(i[10])]) for i in anchout[anchout.CBC==j].itertuples()]))

if __name__=="__main__":
    inputDir=sys.argv[1]
    FI=inputDir+"/merge_anchovy_v1.csv"
    anchOut=pd.read_csv(FI)
    produceFastas(anchOut, inputDir)
