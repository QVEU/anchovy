#CBCtoFasta.py
#Patrick T. Dolan, Ph.D.
#Unit Chief Quantitative Virology and Evolution Unit
#Usage: python CBCtoFasta.py anchovyOutfile.csv
#
import pandas as pd
import numpy asi np
import sys

def produceFastas(anchout, inputDir):
    print("Converting Anchovy CBCs to fastas for mapping and assembly...")
    for j in np.unique(anchout.CBC):
        print(len(anchout.CBC[anchout.CBC==j]))
        print("--------")
        print(j)
        if len(anchout.CBC[anchout.CBC==j])>=5:
            with open(inputDir+"/"+j+".fa",'w') as OF:
                OF.write("\n".join(["".join([">",str(i[9]),"\n",str(i[10])]) for i in anchout[anchout.CBC==j].itertuples()]))
    print("done.")
if __name__=="__main__":
    inputDir=sys.argv[1]
    FI=inputDir+"/merge_anchovy_v1.csv"
    anchOut=pd.read_csv(FI)
    produceFastas(anchOut, inputDir)
