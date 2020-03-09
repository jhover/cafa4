#!/usr/bin/env python
#
#
#
import os
import sys
gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

import pickle

import pandas as pd
import pronto as pt

from cafa4.cafalib import CAFAPlugin
from cafa4.ontology import GeneOntology



class OntologyPlugin(CAFAPlugin):
    """
    Pipeline object. 
    Takes Pandas dataframe 
        1) Walks up GO tree, fanning out each line with different GOterm, *increasing* cafaprob estimate.  

    Input:  Pandas DataFrame
    Output: Pandas DataFrame
        Added rows
        Added columns: cafaprob

    """
    REPR_ATTRS=['outdir']

    def __init__(self, config):
        """
        
        """
        super(OntologyPlugin, self).__init__(config)
        self.log.debug("Creating ontology with pronto...")
        self.cachepath = f"{self.cachedir}/geneontology.pickle"
        #self.pgo = self.get_ontology()
        #self.pgo = pt.Ontology(os.path.expanduser(self.config.get(self.lkname,'obofile')))
        self.go = GeneOntology(config, config.get(self.lkname,'obofile'))
        self.initprob = float(self.config.get(self.lkname,'initial_probability'))
        self.probstep = float(self.config.get(self.lkname,'probability_step'))
        self.rankprob = self.config.getboolean(self.lkname,'rank_probability')
        self.log.debug("OntologyPlugin initialized.")


    def get_ontology(self):
        o = None
        if os.path.exists(self.cachepath):
            self.log.debug(f"{self.cachepath} exists. Loading...")
            f = open(self.cachepath, 'r')
            o = pickle.load(f)
            f.close()
            self.log.debug(f"Loaded ontology from {self.cachepath}")
        else:
            o = pt.Ontology(os.path.expanduser(self.config.get(self.lkname,'obofile')))
            f = open(self.cachepath, 'wb')
            pickle.dump(o,f)
            f.close()
            self.log.debug(f"Dumped ontology to {self.cachepath}")
        return o

    def execute(self, dataframe):
        return self.nativeexecute(dataframe)

    def nativeexecute(self, dataframe):
        # add column for probability
        dataframe['cafaprob'] = self.initprob

        # fan out for each go term up to root...
        newdfdict= {}
        ix = 0
        for row in dataframe.itertuples():
            self.log.debug("inbound row = %s" % str(row))
            (cafaid, evalue, score, bias, db, proteinacc, protein, species, 
            cafaprot, cafaspec, goterm, goaspect, goevidence, cafaprob) = row[1:]
            self.log.debug(f"Searching for parents of '{goterm}'...")        
            gt = self.go.get_term(goterm)
            
            superclasses = gt.superclasses()
            self.log.debug(f"Got generator of superclasses for {goterm}")
            
            # Add row for existing base goterm:
            newrow = [cafaid, evalue, score, bias, db, 
                         proteinacc, protein, species, cafaprot, 
                         cafaspec, goterm, goaspect, goevidence, cafaprob ]
            newdfdict[ix] = newrow
            ix += 1
            
            # Add row for each parent:
            for sclass in superclasses:
                #try:
                #    getattr(gt, 'cafaprob')
                #except AttributeError:
                #    gt.cafaprob = self.initprob
                #cafaprob = sclass.cafaprob + self.probstep 
                cafaprob = cafaprob + self.probstep
                if cafaprob >= 1.0:
                    cafaprob = 0.99
                newrow = [cafaid, evalue, score, bias, db, 
                          proteinacc, protein, species, cafaprot, 
                          cafaspec, sclass , goaspect, goevidence, cafaprob ]
                newdfdict[ix] = newrow
                ix += 1
                
        self.log.debug(f"Processed {ix} rows for new DF.")
        newdf = pd.DataFrame.from_dict(newdfdict, 
                                       orient='index', 
                                       columns = ['cafaid', 
                                                  'evalue', 
                                                  'score', 
                                                  'bias', 
                                                  'db', 
                                                  'proteinacc', 
                                                  'protein', 
                                                  'species', 
                                                  'cafaprot', 
                                                  'cafaspec',
                                                  'goterm',
                                                  'goaspect',
                                                  'goevidence',
                                                  'cafaprob',
                                                  ] )    
        if self.rankprob:
            pass
        
        self.log.debug(f"New dataframe with {len(newdf.index)} rows.")
        return newdf




    def prontoexecute(self, dataframe):
        # add column for probability
        dataframe['cafaprob'] = self.initprob

        # fan out for each go term up to root...
        newdfdict= {}
        ix = 0
        for row in dataframe.itertuples():
            self.log.debug("inbound row = %s" % str(row))
            (cafaid, evalue, score, bias, db, proteinacc, protein, species, 
            cafaprot, cafaspec, goterm, goaspect, goevidence, cafaprob) = row[1:]
            self.log.debug(f"Searching for parents of '{goterm}'...")        
            gt = self.pgo.get_term(goterm)
            #try:
            #    getattr(gt, 'cafaprob')
            #except AttributeError:
            #    gt.cafaprob = self.initprob
            
            superclasses = gt.superclasses()
            self.log.debug(f"Got generator of superclasses for {goterm}")
            
            # Add row for existing base goterm:
            newrow = [cafaid, evalue, score, bias, db, 
                         proteinacc, protein, species, cafaprot, 
                         cafaspec, goterm, goaspect, goevidence, cafaprob ]
            newdfdict[ix] = newrow
            ix += 1
            
            # Add row for each parent:
            for sclass in superclasses:
                #try:
                #    getattr(gt, 'cafaprob')
                #except AttributeError:
                #    gt.cafaprob = self.initprob
                #cafaprob = sclass.cafaprob + self.probstep 
                cafaprob = cafaprob + self.probstep
                if cafaprob >= 1.0:
                    cafaprob = 0.99
                newrow = [cafaid, evalue, score, bias, db, 
                          proteinacc, protein, species, cafaprot, 
                          cafaspec, sclass.id , goaspect, goevidence, cafaprob ]
                newdfdict[ix] = newrow
                ix += 1
                
        self.log.debug(f"Processed {ix} rows for new DF.")
        newdf = pd.DataFrame.from_dict(newdfdict, 
                                       orient='index', 
                                       columns = ['cafaid', 
                                                  'evalue', 
                                                  'score', 
                                                  'bias', 
                                                  'db', 
                                                  'proteinacc', 
                                                  'protein', 
                                                  'species', 
                                                  'cafaprot', 
                                                  'cafaspec',
                                                  'goterm',
                                                  'goaspect',
                                                  'goevidence',
                                                  'cafaprob',
                                                  ] )    
        if self.rankprob:
            pass
        
        self.log.debug(f"New dataframe with {len(newdf.index)} rows.")
        return newdf

 
    
    