#!/usr/bin/env python
#
#  Takes eval file and outputs various useful info...
#
import os
import sys
gitpath=os.path.expanduser("~/git/cafa4")
sys.path.append(gitpath)

import pandas as pd
import numpy as np

from fastcafa.fastcafa import *

if not len(sys.argv) > 1:
    print("need input file arg")
    sys.exit()
f1maxes = []
    
infiles = sys.argv[1:]
for infile in infiles:
    basename = os.path.basename(infile)
    # print(f"{infile}")
    try:
        df = pd.read_csv(infile, index_col=0, comment="#")
        
        totaltrue = len(df[df.correct == True])
        totalnum = df.shape[0] 
        accuracy = totaltrue / totalnum
        numcids = len(df.cid.unique() )
        f1max = df.groupby('cid')['f1max'].max().mean()  
        f1maxes.append(f1max)
                
        #print(f"num cafaids: {numcids}")
        print(f"{basename}\tf1max: {f1max}")
        #print(f"accuracy: {accuracy}")
    except FileNotFoundError:
        print(f"no such file {infile}")

ar = np.array(f1maxes)
meanf1max = np.mean(ar)

print(f">meanf1max(allfiles):\t\t{meanf1max}")
