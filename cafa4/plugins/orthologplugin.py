#!/usr/bin/env python
#
#
import os
import sys
gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

from cafa4.cafalib import CAFAPlugin
from cafa4 import uniprot
#from cafa4 import quickgo

class OrthologPlugin(CAFAPlugin):
    """
    Pipeline object. Takes Pandas dataframe
    looks up orthologs and GO codes.  
    Input:  Pandas DataFrame
    Output: Pandas DataFrame
    """
    REPR_ATTRS=['outdir','backend']

    def __init__(self, config):
        """
        
        """
        super(OrthologPlugin, self).__init__(config)
        self.configname = self.__class__.__name__.lower()
        self.backend = config.get(self.configname ,'backend').strip()

        # self.exc_ec_list = [i.strip() for i in config.get(self.configname, 'excluded_evidence_codes').split(',')] 
        
      
        
    def execute(self, dataframe):
        """
        within cafalib pipeline...
        
        for each row of dataframe, look up ortholog in uniprot and for each GO code
        add a new row with gene, goterm, goaspect
        
        iterate input df fully, putting new info in new df. 
        merge old + new df, return resulting dataframe
        
        """
        self.log.info("Looking up each ortholog")

        newdf = None
        
        if self.backend == 'uniprot':
            self.uniprot = uniprot.UniProt(self.config) 
            self.log.debug("Calling uniprot back end.")
            newdf = self.uniprot.cafa_execute(dataframe)
            
        if self.backend == 'quickgo':
            self.quickgo = quickgo.QuickGOPlugin(self.config)
            self.log.debug("Querying QuickGo for %d unique entries" % len(entries))
            udf = self.quickgo.get_df(entries)
            self.log.debug(qdf)
            qdf.to_csv("%s/quickgo.csv" % self.outdir)
        
        self.log.debug("\n%s" % str(newdf))
        return newdf
    