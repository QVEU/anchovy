#CBCtoFasta.py
#Patrick T. Dolan, Ph.D.
#Unit Chief Quantitative Virology and Evolution Unit
#Usage: python CBCtoFasta.py inputDirectory anchovyOutfile.csv
#
import pandas as pd
import numpy as np
import sys

def produceFastas(anchout, inputDir):
    print("Converting Anchovy CBCs to fastas for mapping and assembly...")
    for j in np.unique(anchout.CBC):
        #print(len(anchout.CBC[anchout.CBC==j]))
        #print("--------")
        #print(j)

        if len(anchout.CBC[anchout.CBC==j])>=5:
            with open(inputDir+"/"+j.strip()+".fa",'w') as OF:
                OF.write("\n".join(["".join([">",str(i[9]),"\n",str(i[10])]) for i in anchout[anchout.CBC==j].itertuples()]))
    print("done.")
if __name__=="__main__":
    inputDir=sys.argv[1]
    inputFile=sys.argv[2]
    FI=inputDir+"/"+inputFile
    anchOut=pd.read_csv(FI)
    produceFastas(anchOut, inputDir)
