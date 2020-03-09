#!/usr/bin/env python
#
#  Takes filenames on input splits it. 
#  Usage: fastasplit.py  <infile>
#    Outputs   <infile>.0 - <infile>.X   each with <sperf> sequences. 
#
import os
import sys

gitpath=os.path.expanduser("~/git/cafa4")
sys.path.append(gitpath)

usage=" Usage: fastasplit.py  <infile>"
sperf = 1800
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

for line in lines:
    line=line.strip()
    #print(f"line='{line}'")
    if line.startswith(">"):
        snum += 1
    if snum == sperf:
        of.close()
        fnum += 1
        # set up for next file...
        ofilename = f"{infile}.{fnum}"
        of = open(ofilename, 'w') 
        snum = 0
    of.write(f'{line}\n')

of.close()




 