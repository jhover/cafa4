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
    print("need input file arg(s)")
    sys.exit()
f1maxes = []
# cgid,cid,correct,f1max,goterm,nterms,score,f1maxtotal,pr
# topdf = pd.DataFrame(columns=['cid','goterm','score','cgid'])

topdf = None  
  
infiles = sys.argv[1:]
for infile in infiles:
    basename = os.path.basename(infile)
    # print(f"{infile}")
    try:
        df = pd.read_csv(infile, index_col=0, comment="#")
        if topdf is None:
    	    topdf = pd.DataFrame(columns=list(df.columns))   
        topdf = topdf.append(df, ignore_index=True, sort=True)
    
    except:
    	print('something went wrong')

topdf.to_csv(sys.stdout)


