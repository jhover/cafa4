#!/usr/bin/env python
#
#     Framework for building/automating CAFAx pipelines. 
#
#     Top level  Runplugins      Infoplugins
#     CAFA4Run
#                PhmmerPlugin                         
#                OrthologPlugin                       looks up ortholog info
#                                QuickGOPlugin        gets GO info
#                                UniprotGOPlugin      gets GO info
#                Expression
#
# 
#      ************************CANONICAL DATAFRAME COLUMNS*********************************
#   
# COLUMN        DESCRIPTION               MAPPINGS                             EXAMPLES
# cafaid        cafa4 target identifier   N/A                                  T2870000001
# cafaprot      cafa4 target protein id                                        1433B
# cafaspec      cafa4 target protien species                                   MOUSE
# proteinid     UniProtKB: entry/ID                                            1433B_RAT
# proteinacc    UniProtKB: accession  ?quickgo:gene product_id                 P63103
# protein       all caps name                                                  1433B
# gene          Free-text gene name.                                           Lrrk2  Ywahb
# geneid        Gene name+species.                                             LRRK2_MOUSE     
# taxonid       NCBI taxon id                                                  9606                 
# species       all caps code                                                  MOUSE   PONAB
# goterm        annotated term                                                 GO:0005634
# goaspect      biological process|molecular function|cellular component       bp       mf   cc
# goevidence    evidence codes for GO annotation.                              IEA 
# evalue        BLAST/HMMER/PHMMER expect statistic                                               1.000000e-126
# bias          Adjustement to score for char prevalence                       3.5
# score         BLAST/HMMER/PHMMER bit-score                                   400.3
# db            
# probest       Probability estimate for prediction.                           0.68  [.000001-1.0]  
#
#

__author__ = "John Hover"
__copyright__ = "2019 John Hover"
__credits__ = []
__license__ = "Apache 2.0"
__version__ = "0.99"
__maintainer__ = "John Hover"
__email__ = "hover@cshl.edu"
__status__ = "Testing"

import os
import sys
gitpath=os.path.expanduser("~/git/cafa4")
sys.path.append(gitpath)

import argparse
from configparser import ConfigParser
from importlib import import_module
import logging
import tempfile

import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import subprocess

    
def get_plugin(klassname):
    """
     load plugin by name from sub-directory

    """
    mname = "cafa4.plugins.%s" %  klassname.lower()
    logging.debug("attempting import of module %s" % mname)
    mod = import_module(mname)
    logging.debug("got module, getting klass=%s" % klassname)
    klass = getattr(mod, klassname)
    logging.debug("got klass, returning...")
    return klass


class CAFA4Run(object):
    
    def __init__(self, config, targetfile, name=None):
        """
        Embodies all the processing for a single run against all targets.
        Overall input is a set of Target sequence files. 
        Overall output is a properly-formatted CAFA4 prediction file.   
        
        """
        self.config = config
        self.targetfile = targetfile
        if name is None:
            self.name = config.get('global','name')
        else:
            self.name = name
            self.config.set('global','name',self.name)
        self.log = logging.getLogger('CAFA4Run')
        self.outdir = os.path.expanduser( config.get('global','outdir') )
        self.author = config.get('global','author')
        self.modelnumber = config.get('global','modelnumber')
        self.keywords = config.get('global','keywords')
        self.profile = config.get('global','profile')
        self.pipeline = [ x.strip() for x in config.get( self.profile,'pipeline').split(',')]
        self.outbase = config.get( self.profile,'outbase')
    
    
    
    def __repr__(self):
        s = "CAFA4Run:"
        for atr in ['name', 'outdir','targetfile','pipeline']:
            s += " %s=%s" % (atr, self.__getattribute__(atr))
        return s


    def cafafile(self, dataframe):
        """
        Produce properly-formated CAFA submission file:
        
        E.g. filename:    gillislab_1_10090_go.txt
                          gillislab_1_78_go.txt   
        
        AUTHOR          GILLIS_LAB
        MODEL          1
        KEYWORDS       orthologs, phmmer
        ACCURACY       1  PR=0.86; RC=0.30
        T100900000001  GO:0042203    .85
        T100900000002  GO:0003998    .03
        .
        .
        .
        
        """
        species = dataframe['cafaspec'][0]
        cafafile = "%s/gillislab_1_%s_go.txt" % (self.outdir, species)
        self.log.debug("Opening cafafile=%s" % cafafile)
        f = open( cafafile, 'w' )
        
        s = ""
        s += f"AUTHOR\t{self.author}\n" 
        s += "MODEL\t%s\n" % 1
        s += "KEYWORDS\tortholog, gene expression\n"
        s += "ACCURACY\t1\tPR=%f;\tRC=%f\n" % (1.00, 1.00) 
        self.log.debug("dataframe columns=%s" % dataframe.columns )
        #for row in dataframe.iterrows():
        #    target = row[1]['cafaid']
        #    goterm = row[1]['goterm']
            #probest = float(row[1]['probest'])  # 1.00 for call/no-call ->  precision/recall *point*. 
        #    probest = row[1]['cafaprob'] # rounding not needed. format does it correctly below. 
        #    s += "%s\t%s\t%.2f\n" % (target, goterm, probest)
        for row in dataframe.itertuples():
            target = row.cafaid
            goterm = row.goterm
            #probest = float(row[1]['probest'])  # 1.00 for call/no-call ->  precision/recall *point*. 
            probest = row.cafaprob # rounding not needed. format does it correctly below. 
            s += f"{target}\t{goterm}\t{probest:.2f}\n"        
        s+="END\n"
        
        f.write(s)
        f.close()
        self.log.info(f"Wrote cafafile with {len(dataframe.index)} entries. ")
        return s
    
    def execute(self):
        self.log.info("Begin run...")
        
        pk = get_plugin('PhmmerPlugin')
        
        phm = pk(self.config, self.targetfile)
        self.log.debug(phm)
        df = phm.execute()
        
        self.log.info("\n%s" % str(df))
        df.to_csv("%s/%s-%s-phmmer.csv" % (self.outdir, self.name,  self.outbase))
        
        ok = get_plugin('OrthologPlugin')
        ortho = ok(self.config)
        self.log.debug(ortho)
        df = ortho.execute(df)
        self.log.info("\n%s" % str(df))
        df.to_csv("%s/%s-%s-ortho.csv" % (self.outdir, self.name, self.outbase))

        gk = get_plugin('OntologyPlugin')
        go = gk(self.config)
        self.log.debug(go)
        df = go.execute(df)
        self.log.info(f"\n{str(df)}")
        df.to_csv("%s/%s-%s-go.csv" % (self.outdir, self.name, self.outbase))
        
        cfstr = self.cafafile(df)
        lines = cfstr.split('\n')
        summary = '\n'.join(lines[:10] + ['.', '.','.']  + lines[-10:]) 
        self.log.info("\n%s" % summary)
        self.log.info("Ending run...")


class CAFAPlugin(object):
    """
    A plugin that normally takes a Pandas dataframe as input, executes, and produces another DF.
    Returns a DF from disk cached if regen=False
    """
    
    REPR_ATTRS=['outdir','cachedir', 'usecache']
    
    def __init__(self, config):
        self.config = config
        self.kname = self.__class__.__name__
        self.lkname = self.kname.lower()
        self.log = logging.getLogger(self.kname)
        self.outdir = os.path.expanduser( config.get('global','outdir') )
        self.cachedir = os.path.expanduser(config.get(self.lkname ,'cachedir'))
        self.cachefile = "%s.csv" % self.lkname
        self.usecache = config.getboolean(self.lkname, 'usecache')


    def __repr__(self):
        s = "%s:" % self.kname
        for atr in self.__class__.REPR_ATTRS:
            s += " %s=%s" % (atr, self.__getattribute__(atr))
        return s 

    def execute(self, dataframe):
        self.log.warn("You have called a base Plugin. No change.")
        return dataframe

    def _df_from_cache(self):
        self.log.debug("Trying to load from cache: %s" % self.cachepath )  
        try:
            self.df = pd.read_csv(self.cachepath, index_col=0)
        except FileNotFoundError:
            self.df = None 
            
    def get_df(self):
        pass
        


   
if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    
    parser = argparse.ArgumentParser()
      
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')

    parser.add_argument('infile', 
                        metavar='infile', 
                        type=str, 
                        help='a .fasta sequence files')
    
    parser.add_argument('-c', '--config', 
                        action="store", 
                        dest='conffile', 
                        default='~/etc/cafa4.conf',
                        help='Config file path [~/etc/cafa4.conf]')

    parser.add_argument('-n', '--name', 
                        action="store", 
                        dest='runname', 
                        default='default',
                        help='Run-specific identifier to use in file output.')    
                    
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    cp = ConfigParser()
    cp.read(args.conffile)
           
    c4run = CAFA4Run(cp, args.infile, args.runname)
    logging.debug(c4run)
    c4run.execute()
