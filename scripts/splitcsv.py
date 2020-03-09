#!/usr/bin/env python
#
#  Takes filenames on input splits it by CID (field 3)  
#  Usage: splitcsv.py  <infile>
#    Outputs   <infile>.0 - <infile>.X   each with <sperf> sequences. 
#
import os
import sys

gitpath=os.path.expanduser("~/git/cafa4")
sys.path.append(gitpath)

usage=" Usage: splitcsv.py  <infile>"
sperf = 2000
fnum = 0
snum = 0

try:
    infile = sys.argv[1]
    print(f"infile={infile}")
except IndexError:
    print(usage)
    sys.exit()
    
f =  open(infile, 'r')
lines = f.readlines()
f.close()
print(f"{len(lines)} lines in file. ")

ofilename = f"{infile}.{fnum}"
of = open(ofilename, 'w')
currentcid = None

for line in lines:
    line=line.strip()
    (idx, cgid, cid, goterm, score) = line.split(',')
    if currentcid is None:
        #beginning line...
        currentcid = cid
        of.write(','.join([idx, cgid, cid, goterm,score]) + '\n')
        snum += 1
    
    elif cid == currentcid:
        # same cid, not reached limit, proceed unconditionally.
        of.write(','.join([idx, cgid, cid, goterm,score]) + '\n')
    
    elif cid != currentcid and snum < sperf:
        # not reached limit, proceed.
        of.write(','.join([idx, cgid, cid, goterm,score]) + '\n')
        snum += 1
        currentcid = cid        
    
    elif cid != currentcid and snum >= sperf :
        # switch file. close out old..
        print(f"{snum} >= {sperf}, new file...")
        of.close()
        fnum +=1
        # setup new...
        ofilename = f"{infile}.{fnum}"    
        of = open(ofilename, 'w')
        currentcid = cid
        of.write(','.join([idx, cgid, cid, goterm,score]) + '\n')
        snum = 1
    
    else:
        print("logic wrong!!")

of.close()

