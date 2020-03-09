#!/usr/bin/env python
#
#  Look up GO annotations from UniProtKB:entrys. 
#  E.g.  P35213  (1422B_RAT  gene: Ywhab
# 

import configparser
import logging
import requests
import sys
from io import  StringIO

import pandas as pd
import pdpipe as pdp
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import subprocess


from bioservices.uniprot import QuickGo # to query online REST interface
from Bio import SeqIO   # to parse uniprot.dat


class QuickGOPlugin(object):

    '''
GENEPRODUCTDB GENEPRODUCTID SYMBOL QUALIFIER   GOTERM     GOASPECT ECO ID      GOEVCODE  REFERENCE       WITH/FROM         TAXONID ASSIGNEDBY ANNOTATIONEXTENSIONDATE
UniProtKB     A4K2U9        YWHAB  involved_in GO:0045744 P        ECO:0000250 ISS       GO_REF:0000024  UniProtKB:P31946  9601    UniProt    20160330

  db          proteinid     gene   goqual      goterm     goaspect  ecoid      goevidence goref          withfrom          taxonid assignby  extdate
    
    '''
   
    def __init__(self, config):
        self.log = logging.getLogger('QuickGOPlugin')
        self.requestbase = "https://www.ebi.ac.uk/QuickGO/services/annotation"
        self.qg = QuickGo()
        self.config = config

    def get_df(self,dataframe):
        '''
        Takes 
        Returns Pandas DataFrame
        
        '''
        entries = dataframe['tacc'].unique().tolist()
        txt = self._query_entries(entries)


    def _query_entries(self, entrylist):
        self.log.debug("querying entry list: %s" % entrylist)
        entrystr=','.join(entrylist)
        self.log.debug("querying entry string: %s" % entrystr)
        requestURL = "%s/downloadSearch?geneProductId=%s" % (self.requestbase, entrystr )
        self.log.debug("RequestURL=%s"% requestURL )        
        r = requests.get(requestURL, headers={ "Accept" : "text/tsv"})
        if not r.ok:
            r.raise_for_status()
            #sys.exit()
        response = r.text
        self.log.debug("response=%s" % response)
        return response

    def _query_taxon(self, taxon):
        self.log.debug("querying taxon: %s" % taxon)      
        requestURL = "%s/downloadSearch?taxonId=%s" % (self.requestbase, taxon )
        self.log.debug("RequestURL=%s"% requestURL ) 
        r = requests.get(requestURL, headers={ "Accept" : "text/tsv"})

        if not r.ok:
            r.raise_for_status()
            #sys.exit()
        response = r.text
        self.log.debug("response=%s" % response)
        return response        


    def _query_bioservices_taxon(self, taxon):
        pass

if __name__ == '__main__':
    config = configparser.ConfigParser()
    
    qg = QuickGOPlugin(config)
    taxon='4577'
    #    entrylist = ['Q9CQV8', 'P35213', 'A4K2U9', 'P31946', 'Q4R572', 'P68250']
    #    out = qg._query_entries(entrylist)
    out = qg._query_taxon(taxon)
    sio = StringIO(out)
    df = pd.read_table(sio, 
                  names=['db','proteinid','gene','goqual','goterm','goaspect','ecoid',
                          'goevidence','goref','withfrom','taxonid','assignby','extdate' ],
                  skip_blank_lines=True,
                  comment='#',
                  skiprows=1,
                  index_col=False,
                  )
    df.to_csv('%s.quickgo.csv' % taxon)
    print(df)    






