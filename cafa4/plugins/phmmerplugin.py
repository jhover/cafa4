#!/usr/bin/env python
#
#
import os
import sys
gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

import argparse
from configparser import ConfigParser
import logging
import tempfile
import traceback

import pandas as pd
import pdpipe as pdp
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import subprocess

from cafa4.cafalib import CAFAPlugin


class PhmmerPlugin(CAFAPlugin):
    """
    Pipeline object. Takes list of Fasta files to run on, returns pandas dataframe of 
    similar sequences with score. 
    Input:   List of FASTA files 
    Output:  Pandas DataFrame
    """
    
    REPR_ATTRS=['targetfile','database', 'score_threshold','cpus']

    def __init__(self, config, targetfile):
        super(PhmmerPlugin, self).__init__(config)
        self.targetfile = targetfile
        self.database = os.path.expanduser( config.get(self.lkname, 'database') )
        self.score_threshold = config.get(self.lkname,'score_threshold', fallback=None)
        self.eval_threshold = config.get(self.lkname, 'eval_threshold', fallback=None)
        self.cpus = config.get(self.lkname,'cpus')
        

    def execute(self):
        self._get_targetinfo()
        outfile = self.run_phmmer_file(self.targetfile , self.database)
        self.log.debug("phmmer outfile=%s" % outfile)
        df = self.read_phmmer_table(outfile)
        #self.log.debug("dflist is %s" % outdfs)
        #df = pd.concat(outdfs, ignore_index=True)
        self.log.debug(f"\n{df}")
        return df
            
    def _get_targetinfo(self):
        """
        Reads self.targetlist files, pulls out 
        target entry info 
        Returns dict of list.
           <targetid> : [<cafaprotein>, <cafaspecies> ]
        
        """
        self.targetinfo={}
        filehandle = None
        try:
            filehandle = open(self.targetfile, 'r')
            for line in filehandle:
                if line.startswith('>'):
                    (targetid , prot_spec) = line.split()
                    targetid = targetid[1:]
                    (cafaprot, cafaspec) = prot_spec.split('_')
                    self.targetinfo[targetid] = [cafaprot, cafaspec]
                else:
                    pass
               
        except Exception as e:
            traceback.print_exc(file=sys.stdout)                
        
        finally:
            if filehandle is not None:
                filehandle.close()
        
        self.log.debug("Parsed file(s) with %d targets" % len(self.targetinfo) )
        return self.targetinfo
    
    
    def run_phmmer_files(self, targetlist, database="/data/hover/data/uniprot/uniprot_sprot.fasta"):
    #
    #  time phmmer --tblout 7955.phmmer.2.txt --cpu 16  --noali 
    #              ~/data/cafa4/TargetFiles/sp_species.7955.tfa 
    #              ~/data/uniprot/uniprot_sprot.fasta 
    #              > 7955.phmmer.console.out 2>&1
        
        dbase = database
        outlist = []
        for file in targetlist:
            if not os.path.exists(file):
                raise FileNotFoundError()
            (cmd, outfile) = self._make_phmmer_cmdline(file)
            self.log.debug("Running cmd='%s' outfile=%s " % (cmd, outfile))
            cp = subprocess.run(cmd, 
                                shell=True, 
                                universal_newlines=True, 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE)
            
            outlist.append(outfile)
            self.log.debug("Ran cmd='%s' outfile=%s returncode=%s " % (cmd,outfile, cp.returncode))
        return outlist

    def run_phmmer_file(self, targetfile, database="/data/hover/data/uniprot/uniprot_sprot.fasta"):
      
        dbase = database
        if not os.path.exists(targetfile):
            raise FileNotFoundError()
        (cmd, outfile) = self._make_phmmer_cmdline(targetfile)
        self.log.info("Running cmd='%s' outfile=%s " % (cmd, outfile))
        cp = subprocess.run(cmd, 
                            shell=True, 
                            universal_newlines=True, 
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE)
        self.log.debug("Ran cmd='%s' outfile=%s returncode=%s " % (cmd,outfile, cp.returncode))
        return outfile

            
    def _make_phmmer_cmdline(self, filename):
        outpath = os.path.dirname(filename)
        filebase = os.path.splitext(os.path.basename(filename))[0]
        outfile = "%s/%s.phmmer.tbl.txt" % (self.outdir, filebase)
        #self.log.debug("outfile=%s" % outfile)
        cmdlist = ['time', 'phmmer']
        cmdlist.append( '--tblout  %s ' % outfile )
        cmdlist.append('--noali' )
        cmdlist.append('--cpu %s ' % self.cpus)
        if self.score_threshold is not None:
            cmdlist.append("-T %s " % self.score_threshold)
        elif self.eval_threshold is not None:
            cmdlist.append("-E %s " % self.eval_threshold)
        
        cmdlist.append(' %s ' % filename )
        cmdlist.append(' %s ' % self.database )
        cmd = ' '.join(cmdlist).strip()
        #self.log.debug("command is '%s'" % cmd)
        return (cmd, outfile)

    def read_phmmer_table(self, filename):
        self.log.info("Reading %s" % filename)
        df = pd.read_table(filename, 
                         names=['target','t-acc','cafaid','q-acc',
                                'evalue', 'score', 'bias', 'e-value-dom','score-dom', 'bias-dom', 
                                'exp', 'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc', 'description'],
                         skip_blank_lines=True,
                         comment='#',
                         index_col=False,
                         skiprows=3,
                         engine='python', 
                         sep='\s+')
        self.log.debug(str(df))
        self.log.debug("Dropping unneeded columns..")
        df = df.drop(['t-acc', 'q-acc','e-value-dom','score-dom', 'bias-dom', 'exp', 
                 'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc', 
                 'description'] , axis=1)
        self.log.debug("Parsing compound fields to define new columns...")
        self.log.debug("Splitting first field for db")
        df['db'] = df.apply(lambda row: row.target.split('|')[0], axis=1)
        self.log.debug("Splitting first field for target accession")
        df['proteinacc'] = df.apply(lambda row: row.target.split('|')[1], axis=1)
        self.log.debug("Splitting first field for prot_species")
        df['prot_spec'] = df.apply(lambda row: row.target.split('|')[2], axis=1)
        self.log.debug("Splitting protein_species to set protein")
        df['protein'] =   df.apply(lambda row: row.prot_spec.split('_')[0], axis=1)
        self.log.debug("Splitting protein_species to set species")
        df['species'] =   df.apply(lambda row: row.prot_spec.split('_')[1], axis=1)
        self.log.debug("Dropping split columns...")
        df.drop(columns=['target','prot_spec'], axis=1, inplace=True)
        self.log.debug("Adding target metadata from input files")
        df['cafaprot']= df.apply(lambda row:  self.targetinfo[row.cafaid][0], axis=1)      
        df['cafaspec']= df.apply(lambda row:  self.targetinfo[row.cafaid][1], axis=1)          
        self.df = df
        return df