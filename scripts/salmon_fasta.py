#!/usr/bin/env python
#
#  Takes filenames on input splits it by CID (field 3)  
#  Usage: splitcsv.py  <infile>
#    Outputs   <infile>.0 - <infile>.X   each with <sperf> sequences. 
#
import os
import sys
import logging

gitpath=os.path.expanduser("~/git/cafa4")
sys.path.append(gitpath)

print("salmon!")

import pandas as pd

#unisalmon='/Users/jhover/play/jones/small-uniprot-trembl-salmon.8030.fasta'
unisalmon='/Users/jhover/play/jones/uniprot-trembl-salmon.8030.fasta'
salmonhi='/Users/jhover/play/jones/salmon_hipriority_uniprotIDs.tsv'
sequencefile='/Users/jhover/play/jones/salmon_hiprio_gillis.tfa'


def parse_tfa_file(infile):
    """
    Reads .tfa file, determines species, target ids, geneids. 
    
    returns dataframe:
    pacc pid sequence
    
    """
    listoflists = []
    try:
        f = open(infile, 'r')
        currentheader = None
        currentseq = None
                
        for line in f:
            #  tr|B5X0T5|B5X0T5_SALSA Glutamate dehydrogenase OS=Salmo salar OX=8030 GN=DHE3 PE=2 SV=1
            if line.startswith(">"):
                # new entry, deal with old:
                if currentseq is not None:
                    logging.debug(f"header: {currentheader}")
                    logging.debug(f"sequence: {currentseq}")
                    currentheader.append(currentseq)
                    listoflists.append(currentheader)
                    currentseq = None
                # start new entry
                currentheader = line[1:].split()
                currentheader = currentheader[0].split('|')[1:]
            else:
                line = line.strip()
                if currentseq is not None:
                    currentseq = f"{currentseq}{line}"
                else:
                    currentseq = line
    except FileNotFoundError:
        logging.error(f"file not readable {infile} ")
        
    logging.debug(f"got {len(listoflists)} entries...")
    logging.debug(f"listoflists: {listoflists}") 
    df = pd.DataFrame(listoflists, columns=['pacc','pid','sequence'])
    #df = pd.DataFrame(priormatrix, index=ontobj.gotermlist, columns=['score'])
    #df.set_index('pacc', inplace=True)
    logging.debug(f"dimension:  {df.shape}")
    return df
    
def read_species(infile):
    sdf=pd.read_csv(infile, sep='\t')
    sdf.rename(columns={'Entry name': 'pid'}, inplace=True)
    sdf = sdf.drop(['Status', 'Organism','Protein names'], axis=1)
    return sdf
    
def write_tfa_file(df, outfile):
    s=""
    snum = 1
    header="G8030"
    x = 60
    
    for index, row in df.iterrows():
       
        s += f">G8030{snum:08} {row.pid}\n"
        chunklist = [ row.sequence[y-x:y] for y in range(x, len(row.sequence)+x, x) ] 
        #logging.debug(f"chunklist {chunklist}")
        for c in chunklist:
            s += f"{c}\n"
        snum += 1
    logging.debug(s)   
    
    try:
        f = open(outfile, 'w')
        f.write(s)
        logging.debug(f"Wrote TFA sequence to file {outfile}")
    except IOError:
        logging.error(f"could not write to file {outfile}")
        traceback.print_exc(file=sys.stdout) 
    finally:
        f.close()
            

if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.INFO)
    #logging.getLogger().setLevel(logging.DEBUG)
    
    tdf = parse_tfa_file(unisalmon)
    print(f"{tdf}\n")
    
    sdf = read_species(infile=salmonhi)
    print(sdf)
    #newdf = 
    
    # merge dataframes joining on pid == name
    mdf = sdf.merge(tdf, how='inner',on=['pid'])
    print(f"{mdf}\n")

    write_tfa_file(mdf, sequencefile)
    