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

  
infiles = sys.argv[1:]
for infile in infiles:
    basename = os.path.basename(infile)
    print(f"handling {infile}")
    try:
        df = pd.read_csv(infile, index_col=0, comment="#")
        df['method'] = 'phmmer'
        df['aspect'] = 'all'
        df.to_csv(infile)
    
    except:
        print('something went wrong')
