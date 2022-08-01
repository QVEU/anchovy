#CellConsensus.py
import pandas as pd
import numpy as np
def mapCellConsensus(anchout):
    print("Mapping cell consensus")
    for j in np.unique(anchout.CBC):
        print(j)
        with open("reads.fa",'w') as OF:
            OF.write("\n".join(["".join([">",str(i[9]),"\n",str(i[10])]) for i in anchout[anchout.CBC==j].itertuples()]))


if __name__=="__main__":
    FI="/nethome/dolanpt/lab_share/Sequencing_Data/QVEU_Seq_0010_Minion_ndas10xlib2/no_sample/20220727_2153_MC-113212_FAT27957_6eb71326/fastq_pass/barcode01/merge_anchovy_v1.csv"
    anchOut=pd.read_csv(FI)
    print(anchOut.loc[8])
    mapCellConsensus(anchOut)
