#!/usr/bin/env python
#
#   ************************CANONICAL COLUMNS*********************************
#   
# COLUMN    DESCRIPTION               MAPPINGS                         EXAMPLES
# cid       cafa4 target identifier   N/A                              T2870000001
# cafaprot  cafa4 target protein id                                    1433B
# cafaspec  cafa4 target protien species                               MOUSE
# cgid      cafa4 target gene id                                       LRRK2_MOUSE 
# pid/gid   UniProtKB: entry/ID                                        LRRK2_MOUSE 
# pacc      UniProtKB: accession  ?quickgo:gene product_id             P63103
# protein   all caps name                                              1433B
# gene      Free-text gene name.  uppercased for consistency?          Lrrk2  Ywahb   
# taxonid   NCBI taxon id                                              9606                 
# species   all caps code                                              MOUSE   PONAB
# goterms   (in dictionaries, not DFs)
# goterm    Gene Ontology Annotation ID                                GO:0005634
# aspect    biological process|molecular function|cellular component   bp mf cc
# goev      evidence code for GO annotation.                           IEA 
# eval      BLAST/HMMER/PHMMER expect statistic                        1.000000e-126
# bias      Adjustment to score for char prevalence                   3.5
# pscore    PHMMER bit-score                                           400.3
# db        database against which orthology query is done             sp (swissprot)
# score     General numeric score of a prediction. [ any scalar ]
# pest      Probability estimate for prediction.   [.01-1.0]           0.68    
# seq       Raw (protein) sequence. 
# seqlen    Number of AA

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
from collections import defaultdict
from configparser import ConfigParser
import logging
import math
import pickle
import pprint as pp
import random
import subprocess
import tempfile
import traceback

import pandas as pd

import numpy as np
np.set_printoptions(threshold=400)
from scipy import sparse
from sklearn import metrics

import seaborn as sns
import matplotlib.pyplot as plt

import h5py

mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)


GOASPECTMAP= { 'biological_process' : 'bp',
               'cellular_component' : 'cc',
               'molecular_function' : 'mf',
               'external'           : 'ex'
            }

UPASPECTMAP = { 'C': 'cc',
                'F': 'mf',
                'P': 'bp'
              }


class Ontology(dict):
    """
    gomatrix:  goterm x goterm np.ndarray fully propagated relationships. 
    gotermidx: { <str> goterm : <int> index in gomatrix, ... } 
    gotermlist:
    
    NOTE: dict key order is now assumed stable as of Python 3.7. If this is run on 
    earlier version, unstable key order will cause chaos. 
    
    """
    
    instance = None
    
    def __init__(self, gomatrix, gotermidx, altidx):
        self.kname = self.__class__.__name__
        self.lkname = self.kname.lower()
        self.log = logging.getLogger(self.kname)
        
        # a list of np.arrays
        self.data = gomatrix
        # a dictionary of goterms -> array_index
        self.gotermidx = gotermidx
        # a dictionary of (alternate) goterms -> (real) goterms
        self.altidx = altidx
        # list of go terms
        # **depends on key sort order stability**
        self.gotermlist = list(gotermidx)
        
        
    def __getitem__(self, key):
        #
        val = None
        try:
            val = self.data[self.gotermidx[key]]
        except KeyError:
            #self.log.debug(f"bad primary goterm key {key}, looking up alt. ")
            try:
                realgt = self.altidx[key]
                realidx = self.gotermidx[realgt]
                val = self.data[realidx]
            except KeyError:
                if not math.isnan(key):
                    self.log.error(f"bad altidx key: '{key}' type:{type(key)} should not happen.")
                # otherwise a nan is sometimes expected for unannotated entries. 
        return val

    def __repr__(self):
        #return repr(self.__dict__)
        s = "{"
        for k in self.gotermidx.keys(): 
            s += f"'{k}' : {self.data[self.gotermidx[k]]} "
        s += "}"
        return s 

    def __len__(self):
        return len(self.data)

    def keys(self):
        return self.gotermidx.keys()
    
    def indexof(self, key):
        return self.gotermidx[key]

    def keyof(self, indexid ):
        return self.gotermlist[indexid]


class UniprotByPid(dict):
    """
    dictionary from proteinid -> propagated boolean go vector for that 
    proteinid (i.e. TSSK4_HUMAN, LRRK2_MOUSE). 
    
    Assumes given Ontology object above, provided on construction.
    
    Provides contains() to check if a given proteinid is annotated 
    with a given goterm.  

    """
    
    instance = None
    
    def __init__(self, dict,  ontology):
        super(UniprotByPid, self).__init__(dict)
        self.kname = self.__class__.__name__
        self.lkname = self.kname.lower()
        self.log = logging.getLogger(self.kname)
        self.ontobj = ontology
        
    def contains(self, gid, goterm):
        """
        Return boolean   True if goterm is annotated for that proteinid  
        """
        gv = self[gid]
        #logging.debug(f"govector fo {gid} is {gv}")
        gtidx = self.ontobj.indexof(goterm)
        #logging.debug(f"goterm index for {goterm} is {gtidx}")
        return gv[gtidx]

class UniprotByPacc(dict):
    """
    Same as UniprotByPid but indexed by accession id
    
    dictionary from pacc -> propagated boolean go vector for that 
    proteinid (i.e. TSSK4_HUMAN, LRRK2_MOUSE). 
    
    Assumes given Ontology object above, provided on construction.
    
    Provides contains() to check if a given proteinid is annotated 
    with a given goterm.  

    """
    
    instance = None
    
    def __init__(self, dict,  ontology):
        super(UniprotByPacc, self).__init__(dict)
        self.kname = self.__class__.__name__
        self.lkname = self.kname.lower()
        self.log = logging.getLogger(self.kname)
        self.ontobj = ontology
        
    def contains(self, pacc, goterm):
        """
        Return boolean   True if goterm is annotated for that pacc  
        """
        gv = self[pacc]
        #logging.debug(f"govector fo {gid} is {gv}")
        gtidx = self.ontobj.indexof(goterm)
        #logging.debug(f"goterm index for {goterm} is {gtidx}")
        return gv[gtidx]


class UniprotByGene(dict):
    """
    dictionary from a gene -> propagated boolean go vector for that 
    gene, regardless of species=. 
    
    Assumes given Ontology object above, provided on construction.
    
    Provides contains() to check if a given gene is annotated 
    with a given goterm.  

    """
    
    instance = None
    
    def __init__(self, dict,  ontology):
        super(UniprotByGene, self).__init__(dict)
        self.kname = self.__class__.__name__
        self.lkname = self.kname.lower()
        self.log = logging.getLogger(self.kname)
        self.ontobj = ontology
        
    def contains(self, gid, goterm):
        """
        Return boolean   True if goterm is annotated for that geneid  
        """
        gv = self[gid]
        #logging.debug(f"govector fo {gid} is {gv}")
        gtidx = self.ontobj.indexof(goterm)
        #logging.debug(f"goterm index for {goterm} is {gtidx}")
        return gv[gtidx]

class SpeciesMap(object):
    """
    Provides lookup between NCBI taxon id, short species code, and Latin name. 
    Derived from NCBI speclist.txt.

       data = [
         { taxid     : speccode, ...  },
         { speccode  : taxonid, ...},
         { linnean   : taxonid, ...},
         { taxid     : linnean, ...}
       ]
    
    """
    
    instance = None
    
    def __init__(self, listofdicts ):
        self.kname = self.__class__.__name__
        self.lkname = self.kname.lower()
        self.log = logging.getLogger(self.kname)
        self.data = listofdicts
        
    def get_speccode(self, taxid):
        """
            Get NCBItaxid -> short code.     
        """
        r = None
        try:
            r = self.data[0][taxid]
        except KeyError:
            self.log.warn("No code for that taxon ID.")
        return r


    def get_taxid_from_linnean(self, linnean):
        """
        
        """        
        r = None
        try:
            r = self.data[2][linnean] 
        except KeyError:
            self.log.warn(f"No taxon ID for '{linnean}'. Check case, space")
        return r        
        
        
    def get_taxonid(self, speccode):
        r = None
        try:
            r = self.data[1][speccode]
        except KeyError:
            self.log.warn("No code for that taxon ID.")        
        return r

    def get_linnean(self, taxonid):
        """
        
        """        
        r = None
        try:
            r = self.data[3][taxonid] 
        except KeyError:
            self.log.warn("No Linnean name for that taxon id.")
        return r   


class ExpressionSet(object):
    """
    Convenience container for aggregate gene expression datasets. 
    Assumes sets are indexed by pacc (UP protein accession id)
    
    Only loads sets requested.
    Previous failures saved as None

    Single instance
    
    """

    instance = None
    
    def __init__(self, config):
        self.config = config
        self.kname = self.__class__.__name__
        self.lkname = self.kname.lower()
        self.log = logging.getLogger(self.kname)        
        # indexed by species code 'HUMAN' 'PIG' 'CAEEL'
        self.datasets = {}
        ExpressionSet.instance = self
    
    def get_dataset(self, scode):
        ds = None
        try:
            ds = self.datasets[scode]
            logging.debug(f"Found existing dataset for {scode}")
        except KeyError:
            logging.debug(f"Loading dataset for {scode}")
            try:
                ds = get_expression_dataset(self.config, scode)
            except Exception:
                logging.warning('problem in get_expression_dataset(). leaving DS None')
            self.datasets[scode] = ds                    
        return ds 
        

def get_uniprot_by_object(config, usecache=True, by='pid', goaspect=None):
    goobj = get_ontology_object(config, usecache=True)
    ubxdict = get_uniprot_by(config,by=by, usecache=True)
    bystr = by.capitalize()
    klassname= f"UniprotBy{bystr}"
    kls = globals()[klassname]
    logging.debug(f"Creating {klassname} object... ")
    uobj = kls(ubxdict, goobj )
    kls.instance = uobj
    return uobj

        
def get_ontology_object(config, usecache=True):
    if Ontology.instance is None:
        build_ontology(config, usecache)
    return Ontology.instance


def get_uniprot_by(config, by='pid', usecache=True, version='2019'):
    """
    loads from cache, or triggers build...
    supports by pid (proteinid) or gene (gene name) or accession (pacc)
    
    """
    logging.debug(f"usecache={usecache}")
    cachedir = os.path.expanduser(config.get('uniprot','cachedir'))
    cachefile = f"{cachedir}/uniprotby{by}.{version}.pickle"    
    ubxdict = None
    if os.path.exists(cachefile) and usecache:
        logging.debug("Cache hit. Using existing info...")    
        try:
            cf = open(cachefile, 'rb')    
            ubxdict = pickle.load(cf)
        except Exception:
            logging.error(f'unable to load via pickle from {cachefile}')
            traceback.print_exc(file=sys.stdout)    
        finally:
            cf.close()       
    else:
        ubxdict = build_uniprot_by(config, by, version)
        logging.debug(f"saving dict: to {cachefile}")
        try:
            cf = open(cachefile, 'wb')    
            pickle.dump(ubxdict, cf )
            logging.debug(f"saved dict: to {cachefile}")
        except Exception as e:
            logging.error(f'unable to dump via pickle to {cachefile}')
            traceback.print_exc(file=sys.stdout)     
        finally:
            cf.close()        
    return ubxdict



def build_uniprot_by(config, by='pacc', version='2019', goaspect=None):
    """

    by = gene | pid | pacc

        pacc species      goterm  gasp goev          pid   gene
0     P0DJZ0   PAVHV  GO:0030430  bp   IDA    11K_PAVHV    11K
1     P32234   DROME  GO:0005525  mf   IDA  128UP_DROME  128UP
2     P83011   SCYCA  GO:0043231  cc   IDA  13KDA_SCYCA    NaN
7143  P24224   ECOLI  GO:0008897  mf   IDA   ACPS_ECOLI   ACPS
7144  P24224   ECOLI  GO:0018070  mf   IDA   ACPS_ECOLI   ACPS    

NOTES:  NaN expected for gene, omit...
        gene will not be unique. handle all...
    
        NaN now expected for goterm, in unannotated paccs. 
    builds dictionary gene -> propagated govector with vectors aggregated by 'by'
    
    Returns dict for later use..
    """
    logging.debug(f"building uniprot by={by} and version={version}")
    exponly_transfer= config.getboolean('uniprot','exponly_transfer')
    logging.debug(f"exponly_transfer={exponly_transfer}")
    ontobj = get_ontology_object(config, usecache=True)
    ubtdf = get_uniprot_byterm_df(config, usecache=True, version=version, goaspect=goaspect)
   
    if exponly_transfer:
        ubtdf = filter_goevidence(config, ubtdf, exponly=True)
    ubtdf = ubtdf[ubtdf[by].notna()]
    ubtdf.sort_values(by=by, inplace=True)
    ubtdf.reset_index(drop=True, inplace=True)
    ubtd = ubtdf.to_dict(orient = 'index')
    logging.debug(f"converted DF to dict: e.g. \n{ [ubtd[i] for i in range(0,3)] } ")
   
    byxdict = {}
    sumreport = 1
    suminterval = 10000
    repthresh = sumreport * suminterval
    gtlength = len(ontobj.gotermidx)
    logging.debug(f"gv length is {gtlength}")
    
    i = 0
    currentxid = None
    currentgv = np.zeros(gtlength, dtype=bool)
    while i < len(ubtd):
        row = ubtd[i]
        xid = row[by]
        #logging.debug(f"row {i} is pid {pid}")
        if currentxid is None:
            currentxid = xid
            ngv = ontobj[row['goterm']]
            if ngv is not None:
                currentgv = currentgv + ngv
            # otherwise ignore...
        elif currentxid == xid:
            # same xid, continue...
            ngv = ontobj[row['goterm']]
            if ngv is not None:
                currentgv = currentgv + ngv 
            # otherwise ignore...
        else:
            # new xid
            byxdict[currentxid] = currentgv 
            currentxid = xid
            currentgv = np.zeros(gtlength, dtype=bool)

        if len(byxdict) >= repthresh:
            logging.info(f"Processed {len(byxdict)} entries... ")
            sumreport +=1
            repthresh = sumreport * suminterval    
        i += 1
        
    samplekeys = list(byxdict.keys())[:3]
    logging.debug(f"Made dict by {by}: {[ (k, byxdict[k]) for k in samplekeys]} ")
    return byxdict

def build_goa_gomatrix_old(config, usecache=False, version='2019', infile=None, outfile=None):
    logging.debug(f"Getting GO ontology object... ")
    ontobj = get_ontology_object(config, usecache=True)
    logging.debug("Parsing GOA file...")
    lol = parse_goa_gaf(config, infile)
    logging.debug(f"List of lists: {len(lol)}")
    
    logging.debug(f"Building matrix...")
    genebygo, genelist = build_genematrix(lol, ontobj)
    logging.debug(f"Done. genebygo: t{type(genebygo)} shape {genebygo.shape} dtype {genebygo.dtype} ")
    #logging.debug("converting to sparse matrix.")
    #genebygo = sparse.lil_matrix(genebygo, dtype=bool)
    
    logging.debug(f"Done. genebygo: t{type(genebygo)} shape {genebygo.shape} dtype {genebygo.dtype} ")
  
    if outfile.endswith(".csv"):
        logging.debug(f"Saving matrix to {outfile}")
        genebygo = genebygo.astype(int)        
        df = pd.DataFrame(genebygo, index=genelist, columns=ontobj.gotermlist)
        logging.debug(f"")
        df.to_csv(outfile, index=True, header=True, sep=',')


    elif outfile.endswith(".tsv"):
        logging.debug(f"Saving matrix to {outfile}")
        genebygo = genebygo.astype(int)        
        df = pd.DataFrame(genebygo, index=genelist, columns=ontobj.gotermlist)
        logging.debug(f"")
        df.to_csv(outfile, index=True, header=True, sep='\t')
            
    logging.debug("Done.")       
    return df
    

def build_genematrix_old(goadata, ontobj):
    '''
    Takes goadata [ [ <gene>, <goterm> ], [ <gene>, <goterm> ] ...]
    produces boolean matrix of genes by goterm vector. 
    
        47k
    
    G1    <vector>
    G2    <vector>
    G3
    
    '''
    logging.debug("In build_genematrix...")
    gotermlist = ontobj.gotermlist  # columnlabels
    genelist = []   # rowlabels
    govectors = []  # data
    # govectors = sparse.lil_matrix( (0,47417),'bool')
    # govectors = []
    
    currentg = None
    currentv = None
    
    for e in goadata:
        (gene, goterm ) = e
        #govect = sparse.lil_matrix( ontobj[goterm])
        govect = ontobj[goterm]
        
        if currentg is None:
            logging.debug(f"First entry. gene is {gene} govect is {govect}")
            currentv = govect
            currentg = gene
            
        elif gene == currentg:
            currentv = currentv + govect
            
        else:
            #logging.debug(f"gene: {gene} != currentgene: {currentg} ")
            genelist.append(gene)
            #if govectors is None:
            #    govectors = govect
            
            #else:
            #    govectors.append(govect)
                #govectors = np.vstack( ( govectors, govect ) )
            govectors.append(currentv)
            currentg = gene
            currentv = govect
    #logging.debug(f"genelist= {genelist}")
    logging.debug(f"collected {len(govectors)} govectors. {len(genelist)} genes.")
    logging.debug("Done building structures. Creating matrix.")
    
    logging.debug(f"govectors[0] = {govectors[0]}")
    
    #m = np.array(govectors)
    m = np.vstack(govectors )
    #return govectors
    return m, genelist

def build_genematrix(goadata, ontobj):
    '''
    Takes goadata [ [ <gene>, <goterm> ], [ <gene>, <goterm> ] ...]
    produces boolean matrix of genes by goterm vector. 
    
          47k
    G1    <vector>
    G2    <vector>
    G3
    
    '''
    logging.debug("In build_genematrix...")
    gotermlist = ontobj.gotermlist  # columnlabels
    govectors = {}  # genename -> govector
    # govectors = sparse.lil_matrix( (0,47417),'bool')
   
    for e in goadata:
        (gene, goterm ) = e
        #govect = sparse.lil_matrix( ontobj[goterm])
        govect = ontobj[goterm]
        
        try:
            gv = govectors[gene]
            govectors[gene] = gv + govect
            
        except KeyError:
            govectors[gene] = govect

    veclist = []  
    genelist = []   # rowlabels    

    for k in govectors.keys():
        veclist.append(govectors[k])
        genelist.append(k)

    #logging.debug(f"genelist= {genelist}")
    logging.debug(f"collected {len(govectors)} govectors. {len(genelist)} genes.")
    logging.debug("Done building structures. Creating matrix.")
    
    #logging.debug(f"govectors[0] = {govectors[0]}")
    
    #m = np.array(govectors)
    m = np.vstack(veclist )
    #return govectors
    return m, genelist
    

def parse_goa_gaf(config, infile=None):
    '''
    create list of lists of GOA database. 
    [ <gene>, <goterm> ]
    '''
    if infile is not None:
        filepath = infile
    else:
        filepath = os.path.expanduser(config.get('goa','datafile'))
    
    try:
        logging.debug(f" attempting to open '{filepath}'")
        filehandle = open(filepath, 'r')
        allentries = []
        current = None
        sumreport = 1
        suminterval = 10000
        repthresh = sumreport * suminterval
        
        while True:
            line = filehandle.readline()
            if line == '':
                break
    
            if line.startswith("!"):
                pass
    
            else:
                fields = line.split('\t')
                fields = fields[:7]
                geneterm = [ fields[2], fields[4] ]
                allentries.append(geneterm)

    except FileNotFoundError:
        logging.error(f"No such file {filepath}")   

    except Exception as e:
        traceback.print_exc(file=sys.stdout)                
    
    if filehandle is not None:
        filehandle.close()      
    
    logging.info(f"Parsed file with {len(allentries)} entries" )
    logging.debug(f"Some entries:  {allentries[1000:1005]}")
    return allentries
        

def do_build_prior(config, usecache=True ):    

    logging.info("calculating prior..")
    priormatrix = calc_prior(config, usecache, species=None)
    logging.info(f"priormatrix = {matrix_info(priormatrix)}")
    priordf = get_prior_df(config, usecache)
    logging.info(f"priordf:\n{priordf}")    


def do_orthoexpression(config, infile, outfile, usecache=True, version='2019', goasp=None):
    """
    Peform phmmer search for hits and pass output to method to perform 
    expression lookup. 
    
    """
    logging.info(f"Running phmmer on {infile}")
    pdf = get_phmmer_df(config, infile)
    logging.info(f"got phmmer df:\n{pdf}")

    if pdf is not None:
        logging.info(f"Making expression prediction against version '{version}'")
        df = calc_orthoexpression_prediction(config, pdf, usecache, version)
        if df is not None:
            logging.debug(f"prediction=\n{df}")#
            df = filter_goaspect(config, df, goasp)
            if df is not None:
                logging.info(f"writing to outfile {outfile}")
                df.to_csv(outfile)
            else:
                logging.debug(f"no data for prediction.")
        else:
            logging.debug(f"no data for prediction.")
    else:
        logging.info(f"No phmmer hits for infile {infile} so can't run expression.")    


def apply_pid2pacc(row, idmap):
    a = np.nan
    try:
        a = idmap[row.cgid]
    except KeyError:
        logging.warning(f"KeyError in idmap: {row.cgid}") 
    return a

def apply_go2aspect(row, godict ):
    """ 
      'GO:0001654': {'is_a': ['GO:0007423'],
      'part_of': ['GO:0150063'],
      'goterm': 'GO:0001654',
      'goname': 'eye development',
      'goasp': 'bp'},
    """
    a = np.nan
    try:
        a = godict[row['goterm']]['goasp']
        #logging.debug(f"got aspect {a} from row:\n{row}")
    except KeyError:
        logging.warning(f"KeyError in godict: {row.goterm}") 
    return a


def do_expression(config, infile, outfile, usecache=True, version='2019',goasp=None ):
    """
    Infile can be FASTA with second ID field <gene>_<species>
    E.g. YEFM_SALTY, CA194_HUMAN
    
    """
    logging.info(f"Doing expression on file {infile}...")
    idmaps = get_idmaps(config, usecache=usecache, version=version)
    p2a = idmaps['pid2pacc']
    
    infile = os.path.expanduser(infile)
    if os.path.exists(infile):
        indf = parse_tfa_file(infile)
        logging.debug(f"indf=\n{indf}")
        # add column for accession no map: pacc 
        indf['pacc'] = indf.apply(apply_pid2pacc, axis=1, idmap=p2a)
        indf['score'] = 1.0
        indf['pid'] = indf['cgid']

        df = calc_expression_prediction(config, indf, usecache, version)
        logging.debug(f"prediction=\n{df}")
        if df is not None:
            logging.debug(f"Adjusting prediction limit to {goasp}")
            df = filter_goaspect(config, df, goasp)            
            if df is not None:
                logging.info(f"writing to outfile {outfile}")
                df.to_csv(outfile)
            else:
                logging.info("No content. Not writing outfile.")
        else:
            logging.info("No content. Not writing outfile.")
        logging.info("done.")
    else:
        logging.error(f"No such infile: {infile}")



def do_phmmer(config, infile, outfile, usecache=True, version='current', goasp=None):
    """
    Perform phmmer on infile sequences. 
    Output prediction to outfile.
    
    Also cache predictions....
    
    """
    logging.info(f"running phmmer version={version} goasp={goasp}")
    infile = os.path.expanduser(infile)
    pcachedir = os.path.expanduser(config.get('phmmer','cachedir')) 
    filename = os.path.basename(infile)
    (filebase, e) = os.path.splitext(filename)
    predcachefile = f"{pcachedir}/{filebase}.phmmer.{version}.{goasp}.pred.csv"
    logging.debug(f"predcachefile={predcachefile}")
    
    if os.path.exists(infile):
        df = None
        if os.path.exists(predcachefile):
            logging.info("Predication cache hit. loading..")
            df = pd.read_csv(predcachefile, index_col=0)
        else:      
            pdf = get_phmmer_df(config, infile, version)
            if pdf is not None:
                logging.info(f"got phmmer df:\n{pdf}")
                logging.info("making phmmer prediction...")
                df = calc_phmmer_prediction(config, pdf, usecache, version)
                logging.debug(f"prediction=\n{df}")
                if df is not None:
                    logging.debug(f"Adjusting prediction limit to {goasp}")
                    df = filter_goaspect(config, df, goasp)
                    logging.info(f"Caching to {predcachefile}")
                    df['method'] = 'phmmer'
                    df.to_csv(predcachefile)
                    logging.info(f"writing to outfile {outfile}")
                    df.to_csv(outfile)
                    logging.info("done.")
                else:
                    logging.warning("No phmmer hits. No output.")
            else:
                logging.warning("No phmmer hits. No output.")
    else:
        logging.error(f"No such infile: {infile}")



def parse_tfa_file( infile):
    """
    Reads .tfa file, determines species, target ids, geneids. 
    
    returns dataframe:
    cid  cgid
    
    """
    listoflists = []
       
    try:
        f = open(infile, 'r')
    except FileNotFoundError:
        logging.error(f"file not readable {filename} ")
    for line in f:
        # >T100900000004 1433G_MOUSE
        if line.startswith(">"):
            fields = line[1:].split()
            cid = fields[0].strip()
            cgid = fields[1].strip()
            listoflists.append( [cid, cgid] )
    logging.debug(f"got {len(listoflists)} cids with geneids.") 
    df = pd.DataFrame(listoflists, columns=['cid','cgid']) 
    return df    

def parse_tfa_file_lol( infile):
    """
    Reads .tfa file, determines species, target ids, geneids. 
    
    returns dataframe:
    cid  cgid
    
    """
    listoflists = []
       
    try:
        f = open(infile, 'r')
        for line in f:
            # >T100900000004 1433G_MOUSE
            if line.startswith(">"):
                fields = line[1:].split()
                cid = fields[0].strip()
                cgid = fields[1].strip()
                listoflists.append( [cid, cgid] )
            logging.debug(f"got {len(listoflists)} cids with geneids.") 
    except FileNotFoundError:
        logging.error(f"file not readable {infile} ")
    
    return listoflists



def do_prior(config, infile, outfile, usecache=True, species=None, version='current',goasp=None):
    """
    Apply prior likelihood (global or species) to all infile sequences. 
    Output prediction to outfile for later eval. 
    
    """
    logging.info("making prior prediction...")
    df = make_prior_prediction(config, infile, species)
    logging.debug(f"prediction=\n{df}")
    logging.debug(f"Adjusting prediction limit to {goasp}")
    df = filter_goaspect(config, df, goasp)        
    if len(df) > 0:
        logging.info(f"writing to outfile {outfile}")
        df.to_csv(outfile)
    else:
        logging.warning(f"No data. Not writing to outfile.")
    logging.info("done.")


  

def run_combine_weighted_avg(config, predict1, predict2, outpred, weight=1.0 ):
    """
    Weighted average rank aggregration. 
    
    Takes two predictions, creates combined prediction, weighting the goterm scores
    Where weight is not 1.0 (average between p1 and p2), the input weight is the second 
    with respect to the first, i.e. the first 'perturbed' by the second. 

    Must be done *per cid*. Missing cids from one or the other are passed through unchanged.  

        ,cgid,    cid,        goterm,    score
    0,AACC3_PSEAI,T2870000002,GO:0003674,182.0
    1,AACC3_PSEAI,T2870000002,GO:0003824,182.0
    2,AACC3_PSEAI,T2870000002,GO:0008080,182.0

  
   
    """
    obj = get_ontology_object(config)
    
    pdf1 = pd.read_csv(predict1, index_col=0)
    pdf2 = pd.read_csv(predict2, index_col=0)
    
    for idf in [pdf1, pdf2]:
        logging.debug(f"df=\n{idf}")
        
    # sanity check
    p1cids = set(pdf1.cid.unique())
    logging.debug(f"input 1 has {len(p1cids)}")
    p2cids = set(pdf1.cid.unique())
    logging.debug(f"input 2 has {len(p2cids)}")
        
    # Dataframe to collect all calculated values. 
    topdf = pd.DataFrame(columns=['cgid', 'cid', 'goterm','score'])

    allcids = list(p1cids | p2cids)
    allcids.sort()
    logging.debug(f"got {len(allcids)} total cids, sorted, e.g. {allcids[0:3]} ")
    

    #logging.debug("Starting to handle each CID...")
    for cid in allcids:     
        #select
        cdf1 = pdf1[pdf1.cid == cid]
        cdf2 = pdf2[pdf2.cid == cid]
        # copy
        #cdf1 = cdf1.copy()
        #cdf2 = cdf2.copy()
        # reset index        
        cdf1.reset_index(drop=True)
        cdf2.reset_index(drop=True)  
        
        # sort       
        cdf1.sort_values(by='score',inplace=True, ascending=False)
        cdf2.sort_values(by='score',inplace=True, ascending=False)

        # drop=True?
        cdf1 = cdf1.reset_index(drop=True)
        cdf2 = cdf2.reset_index(drop=True) 
        
        logging.debug(f"Before ranks: cdf1=\n{cdf1}\ncdf2=\n{cdf2} ")
      
        # set ranks    
        cdf1['rank'] = cdf1.index
        cdf2['rank'] = cdf2.index        

        logging.debug(f"before merge:\ncdf1=\n{cdf1}\ncdf2=\n{cdf2}")       
        cdf = pd.merge(cdf1, cdf2 , how='outer', on=['cid','cgid','goterm'] )
        
        # set NaN to maxrank, minscore
        logging.debug(f"after merge:\n{cdf}")
        xmaxrank = cdf.rank_x.max()
        xminscore = cdf.score_x.min()
        ymaxrank = cdf.rank_y.max()
        yminscore = cdf.score_y.min()
        
        values = { 'rank_x': xmaxrank, 'score_x': xminscore, 'rank_y': ymaxrank , 'score_y' : yminscore }
        cdf.fillna(value=values, inplace=True)
        logging.debug(f"after fillna:\n{cdf}")        
        
        # correct ranks
        cdf.rank_x = cdf.rank_x + 1
        cdf.rank_y = cdf.rank_y + 1    
        logging.debug(f"cdf is:\n{cdf}")
        #cdf['rank'] = ( cdf.rank_x + cdf.rank_y / 2 )
        # Take weighted mean of ranks:
        #
        #     ( 1 * rank_x) + ( weight * rank_y )
        #     -------------------------------
        #                1 + weight
        #
        cdf['rank'] = ( cdf.rank_x + ( cdf.rank_y * weight) ) /  ( 1 + weight ) 

        logging.debug(f"after calc cdf is:\n{cdf}")        
        dropcol = ['score_x', 'rank_x','score_y','rank_y' ] 
        cdf.drop(columns=dropcol, inplace=True)
        logging.debug(f"after drop is:\n{cdf}")
        cdf.sort_values(by='rank', inplace=True, ascending=False)        
        logging.debug(f"after sort cdf is:\n{cdf}")    
        cdf.reset_index(drop=True, inplace=True)
        cdf['score'] = cdf.index + 1
        #logging.debug(f"after reset cdf is:\n{cdf}")
        cdf.sort_values(by='score', inplace=True, ascending=False)
        #logging.debug(f"after resort cdf is:\n{cdf}")
        cdf.reset_index(drop=True, inplace=True)
        #logging.debug(f"after reset cdf is:\n{cdf}")
        cdf.drop(columns='rank', inplace=True)
        logging.debug(f"after reset cdf is:\n{cdf}")
        topdf = topdf.append(cdf, ignore_index=True)

    logging.debug(f"topdf after all=\n{topdf}")
    topdf.to_csv(outpred)
    return topdf




def run_combine_rr(config, predict1, predict2, outpred):
    """
    
    Round-robin rank aggregation. 
    
    Takes two predictions, creates combined prediction, popping the top ranked prediction
    off each input prediction to create the output.  

    Must be done *per cid*. Missing cids from one or the other are passed through unchanged.  

        ,cgid,    cid,        goterm,    score
    0,AACC3_PSEAI,T2870000002,GO:0003674,182.0
    1,AACC3_PSEAI,T2870000002,GO:0003824,182.0
    2,AACC3_PSEAI,T2870000002,GO:0008080,182.0
   
    """
    obj = get_ontology_object(config)
    
    pdf1 = pd.read_csv(predict1, index_col=0)
    pdf2 = pd.read_csv(predict2, index_col=0)
    
    for idf in [pdf1, pdf2]:
        logging.debug(f"df=\n{idf}")
        
    # sanity check
    p1cids = set(pdf1.cid.unique())
    logging.debug(f"input 1 has {len(p1cids)}")
    p2cids = set(pdf1.cid.unique())
    logging.debug(f"input 2 has {len(p2cids)}")
        
    # Dataframe to collect all calculated values. 
    topdf = pd.DataFrame(columns=['cgid', 'cid', 'goterm','score'])

    allcids = list(p1cids | p2cids)
    allcids.sort()
    logging.debug(f"got {len(allcids)} total cids, sorted, e.g. {allcids[0:3]} ")
    

    #logging.debug("Starting to handle each CID...")
    for cid in allcids:     
        #select
        cdf1 = pdf1[pdf1.cid == cid]
        cdf2 = pdf2[pdf2.cid == cid]
        # copy
        cdf1 = cdf1.copy()
        cdf2 = cdf2.copy()
        # reset index        
        cdf1.reset_index(drop=True)
        cdf2.reset_index(drop=True)  
        
        # sort       
        cdf1.sort_values(by='score',inplace=True, ascending=False)
        cdf2.sort_values(by='score',inplace=True, ascending=False)

        # drop=True?
        cdf1 = cdf1.reset_index(drop=True)
        cdf2 = cdf2.reset_index(drop=True) 
        
        logging.debug(f"Before ranks: cdf1=\n{cdf1}\ncdf2=\n{cdf2} ")
       
        print(cdf1)
        print(cdf2) 
        
        cdf = mergerankrr(cdf1, cdf2)
       
        topdf = topdf.append(cdf, ignore_index=True)

    logging.debug(f"topdf after all=\n{topdf}")
    topdf.to_csv(outpred)
    return topdf


def mergerankrr(df1, df2):
    '''
    performs round robin rank aggregation. 
    
           cgid              cid      goterm    score
0     ODPB_MYCGE  G24327300000000  GO:0005575  64753.3
1     ODPB_MYCGE  G24327300000000  GO:0110165  64317.1
2     ODPB_MYCGE  G24327300000000  GO:0008150  61633.3

           cgid              cid      goterm  score
0    ODPB_MYCGE  G24327300000000  GO:0008150  535.2
1    ODPB_MYCGE  G24327300000000  GO:0008152  535.2
2    ODPB_MYCGE  G24327300000000  GO:0032787  356.8    

    creates new df ordered by rank. 
    
    '''
    logging.debug("In mergerankrr...")
    
    cgid = df1['cgid'].iloc[0]
    cid = df1['cid'].iloc[0]
    logging.debug(f"cgid={cgid}")
    logging.debug(f"cid={cid}")
    
    
    l1 = df1['goterm'].tolist()
    l1.reverse()
    l2 = df2['goterm'].tolist()
    l2.reverse()
    logging.debug(f"len l1: {len(l1)} len l2: {len(l2)}")
    
    newlist = list()
    
    while (len(l1) + len(l2)) > 0:
        if len(l1)>0:
            gt = l1.pop()
            if gt not in newlist:
                newlist.append(gt)
            
        if len(l2)>0:
            gt = l2.pop()
            if gt not in newlist:
                newlist.append(gt)        
    
    logging.debug(f'newlist len: {len(newlist)} ')        
    
    cidf = pd.DataFrame(newlist, columns=['goterm'])
    cidf['cgid'] = cgid
    cidf['cid'] = cid
    # when merging distinct predictions, scores are no longer comparable 
    cidf['score'] = 1.0
    
    logging.debug(f"new cidf:\n{cidf}")    
    return cidf
    


def run_evaluate(config, predictfile, outfile, goaspect=None):
    """
    Consume a prediction.csv file, and score based on accuracy. 
    X.eval.csv
   
    """
    df = pd.read_csv(os.path.expanduser(predictfile), index_col=0)
    logging.debug(f"got predictdf types:\n{df.dtypes}\n{df}")
    edf = do_evaluate(config, df, goaspect)

    #max_goterms = config.get('global','max_goterms')
    eval_threshold = config.get('phmmer','eval_threshold')
    topx_threshold = config.get('phmmer','topx_threshold')
    score_method = config.get('phmmer','score_method')
    #logging.info(f"hyperparams:\nmax_goterms={max_goterms}\neval_threshold={eval_threshold}\ntopx_threshold={topx_threshold}\nscore_method={score_method}  ")
    logging.info(f"hyperparams:\neval_threshold={eval_threshold}\ntopx_threshold={topx_threshold}\nscore_method={score_method}  ")
    #logging.debug(f"got evaluation df:\n{edf}")
    fh = open(outfile, 'w')
    fh.write(f"# eval_threshold={eval_threshold}\n# topx_threshold={topx_threshold}\n# score_method={score_method}\n")
    edf.to_csv(fh)
    fh.close()    
    print(edf)



def is_correct_apply(row):
    """
    Function to check prediction for use with .apply and Pandas DF. 

    """
    return UniprotByPid.instance.contains(row.cgid, row.goterm)


def do_evaluate(config, predictdf, goaspect):
    """
    Calculate evalutation numbers based on known terms. 
    
    i    cid           goterm       score    cgid
    0    G960600000001 GO:0086041   53.0   Q9Y3Q4_HUMAN
    1    G960600000001 GO:0086089   49.0   Q9Y3Q4_HUMAN
    2    G960600000001 GO:0086090   49.0   Q9Y3Q4_HUMAN
    
    Return:
    
    eval_threshold=1.0e-120
    topx_threshold=200
    score_method=phmmer_score_weighted  
              cgid             cid  correct      goterm      pauc      pest     score     pr
        CHIA_MOUSE  G1009000000001     True  GO:0008150  0.759799  0.990000  0.046959  0.230651
        CHIA_MOUSE  G1009000000001    False  GO:0005575  0.759799  0.442188  0.020743  0.230651
        CHIA_MOUSE  G1009000000001    False  GO:0110165  0.759799  0.423416  0.019845  0.230651

    """
    #logging.debug(f"got predictdf:\n{predictdf}")
    ubp = get_uniprot_by_object(config, by='pid', usecache=True, goaspect=goaspect)
    ontobj = get_ontology_object(config, usecache=True)
    logging.debug(f"got known uniprot and ontology object.")  
    
    outdf = pd.DataFrame(columns = ['cid','goterm','score','cgid', 'correct', 'nterms'])

    cidlist = list(predictdf.cid.unique())
    logging.debug(f"cid list: {cidlist}")

    ntermsum = 0
    for cid in cidlist:
        cdf = predictdf[predictdf.cid == cid].copy()

        # get gene id. 
        cgid = cdf.cgid.unique()[0]
        logging.debug(f"cgid is {cgid}")
        try:
            nterms = ubp[cgid].sum()
            if nterms > 0:
                ntermsum = ntermsum + nterms
                #logging.debug(f"there are {nterms} goterms associated with {cgid}")
                #logging.debug(f"geneid for this target is is {cgid}")
                cdf['correct'] = cdf.apply(is_correct_apply, axis=1)
                cdf['nterms'] = nterms
                #cdf.reset_index(drop=True, inplace=True) 
                logging.debug(f"cdf after assessment:\n{cdf.dtypes}\n{cdf}")
                #logging.debug(f"cdf is:\n{cdf}")
                #logging.debug(f"appending: outdf.columns={outdf.columns} cdf.columns={cdf.columns}")
                #cdf = calc_f1_max(cdf)
                cdf = calc_f1(cdf)
                logging.debug(f"cid {cid}:\n{cdf}")
                outdf = outdf.append(cdf, ignore_index=True)
            else:
                logging.warning(f"nterms=0 ignoring.??")
        except KeyError:
            logging.warning(f"No entry found for {cgid}...")    
        
    outdf[['gene','species']] = outdf.cgid.str.split('_',expand=True)
    return outdf


def calc_f1_max(dataframe):
    """
    F1 max calculation by protein. 
    Dataframe assumes one cid represented. 
    Calculates f1max for set. 
    Sets all rows 'f1max' to that value.
    
    nterms is total number of correct terms
    
    In:
             cgid             cid  correct      goterm nterms      pest     score  ntermsum
0       CHIA_MOUSE  G1009000000001     True  GO:0008150     85  0.990000  0.046959       764
1       CHIA_MOUSE  G1009000000001    False  GO:0005575     85  0.442188  0.020743       764
2       CHIA_MOUSE  G1009000000001    False  GO:0110165     85  0.423416  0.019845       764
3       CHIA_MOUSE  G1009000000001     True  GO:0009987     85  0.387863  0.018143       764 
    ...
    
    F1 at each index i is:
            pr = num-true  / i 
            rc = num-true / (nterms - num-true)
            2  *   (   pr * rc / pr + rc )
    
    """
    logging.debug(f"inbound frame is:\n{dataframe}")
    nterms = dataframe['nterms'].max()
    logging.debug(f"annotated terms is {nterms}")
    dataframe.sort_values(by='score', ascending=False, inplace=True)
    dataframe.reset_index(drop=True, inplace=True)
    logging.debug(f"frame after sorting:\n{dataframe}")    
    
    numtrue = 0
    i = 1
    f1max = 0
    logging.debug(f"initial. i={i} numtrue={numtrue} f1max={f1max}")
    for row in dataframe.itertuples():
        if row.correct == True:
            numtrue += 1
            logstr = "correct"
        else:
            logstr = " not correct"
        try:
            pr = (numtrue / i)
            #logging.debug(f"pr[{i}] = {pr}")
            rc = (numtrue / nterms )
            #logging.debug(f"rc[{i}] = {rc}")
            # from Cafa tool:   f = (2*pr*rc)/(pr+rc)  SAME VALUE
            f1 = 2 * ( (pr * rc) / (pr + rc ) )
            logging.debug(f"{row.cid}: {row.goterm} {logstr}. f1={f1}")
        except ZeroDivisionError:
            f1 = 0
            logging.warning("Division by zero error during f1 calculation.")
        if f1 > f1max:
            f1max = f1
        i += 1
    logging.debug(f"final. i={i} numtrue={numtrue} f1max={f1max}")
    dataframe['f1max'] = f1max
    return dataframe


def calc_f1(dataframe):
    """
    F1 calculation for each goterm prediction.  
    Dataframe assumes one cid represented. 
        
    nterms is total number of correct terms
    
    In:
             cgid             cid  correct      goterm nterms      pest     score  ntermsum
0       CHIA_MOUSE  G1009000000001     True  GO:0008150     85  0.990000  0.046959       764
1       CHIA_MOUSE  G1009000000001    False  GO:0005575     85  0.442188  0.020743       764
2       CHIA_MOUSE  G1009000000001    False  GO:0110165     85  0.423416  0.019845       764
3       CHIA_MOUSE  G1009000000001     True  GO:0009987     85  0.387863  0.018143       764 
    ...
   
    F1 at each index i is:
            pr = num-true  / i 
            rc = num-true / (nterms - num-true)
            2  *   (   pr * rc / pr + rc )
    
    """
    logging.debug(f"inbound frame is:\n{dataframe} length={dataframe.shape[0]}")
    nterms = dataframe['nterms'].max()
    logging.debug(f"annotated terms is {nterms}")
    dataframe.sort_values(by='score', ascending=False, inplace=True)
    dataframe.reset_index(drop=True, inplace=True)
    logging.debug(f"frame after sorting:\n{dataframe}")    
    
    numtrue = 0
    i = 1
    #f1max = 0
    f1= 0.0 
    logging.debug(f"initial. i={i} numtrue={numtrue} f1={f1}")
    #dataframe['f1'] = 0
    f1list = []
        
    for row in dataframe.itertuples():
        if row.correct == True:
            numtrue += 1
            logstr = "correct"
        else:
            logstr = " not correct"
        try:
            pr = (numtrue / i)
            #logging.debug(f"pr[{i}] = {pr}")
            rc = (numtrue / nterms )
            #logging.debug(f"rc[{i}] = {rc}")
            # from Cafa tool:   f = (2*pr*rc)/(pr+rc)  SAME VALUE
            f1 = 2 * ( (pr * rc) / (pr + rc ) )
            logging.debug(f"{row.cid}: {row.goterm} {logstr}. f1={f1}")
        except ZeroDivisionError:
            f1 = 0
            logging.warning("Division by zero error during f1 calculation.")
        f1list.append(f1)
        #if f1 > f1max:
        #    f1max = f1
        i += 1
    dataframe['f1'] = pd.Series(f1list)
    logging.debug(f"final. i={i} numtrue={numtrue}")
    #dataframe['f1max'] = f1max
    return dataframe



def calc_f1_max_overall(dataframe):
    """
    In:
             cgid             cid  correct      goterm nterms      pest     score  ntermsum
0       CHIA_MOUSE  G1009000000001     True  GO:0008150     85  0.990000  0.046959       764
1       CHIA_MOUSE  G1009000000001    False  GO:0005575     85  0.442188  0.020743       764
2       CHIA_MOUSE  G1009000000001    False  GO:0110165     85  0.423416  0.019845       764
3       CHIA_MOUSE  G1009000000001     True  GO:0009987     85  0.387863  0.018143       764 
    ...
    
    F1 at each index i is:
            pr = num-true  / i 
            rc = num-true / (numtermsum - num-true)
           2  *   (   pr * rc / pr + rc )

    
    """
    logging.debug(f"inbound frame is:\n{dataframe}")
    totalterms = dataframe.groupby('cid')['nterms'].max().sum()
    logging.debug(f"total annotated terms is {totalterms}")
    dataframe.sort_values(by='score', ascending=False, inplace=True)
    dataframe.reset_index(drop=True, inplace=True)
    logging.debug(f"frame after sorting:\n{dataframe}")    
    
    
    
    numtrue = 0
    i = 1
    f1max = 0
    logging.debug(f"initial. i={i} numtrue={numtrue} f1max={f1max}")
    for row in dataframe.itertuples():
        if row.correct == True:
            numtrue += 1
        #logging.debug(f"numtrue is {numtrue}")
        pr = (numtrue / i)
        #logging.debug(f"pr[{i}] = {pr}")
        rc = (numtrue / totalterms )
        #logging.debug(f"rc[{i}] = {rc}")
        f1 = 2 * ( (pr * rc) / (pr + rc ) )
        #logging.debug(f"f1 is {f1}")
        if f1 > f1max:
            f1max = f1
        i += 1
    logging.debug(f"final. i={i} numtrue={numtrue} f1max={f1max}")
    dataframe['f1maxtotal'] = f1max
    return dataframe



def calc_precision_recall(posidxlist, totalnum):
    """
        Calculates precision recall    
    """
    # precision-recall -> pr
    n = len(posidxlist)
    sum = 0
    i = 1
    for posi in posidxlist:
        sum = sum + ( i / ( posi + 1 ) )
        #logging.debug(f"sum is {sum}")
        i += 1
    if n == 0:
        pr = 0.0
    else:
        pr = (1 / n) * sum 
    return pr


def get_evaluate_df(config, predictdf, goaspect=None,  threshold=None ):
    """
    lol:
         cid            numpredict   numcorrect  numannotated
    [['G960600000001',  4317,        539,        1637], 
     ['G960600000002',  10782,       0,          2], 
     ['G960600000003',  5634,        371,        555], 
    ]
    
    """
    
    df = do_evaluate_auroc(config, predictdf, goaspect)
    return df



def do_summarize(config, evaldf):
    """
    Calculate final evaluation one row per target. 
    Input is evaluation dataframe (one row per goterm). 
    
    In:
      i, cid,           goterm,     score, cgid,       correct, nterms, method, aspect, f1
      0, G722700000021, GO:0005575, 177.4, F172A_DROME,False,   48,     phmmer, cc,     0.0
    
    Out:
        cgid             cid            f1max   nterms   ncorrect  pr
        CHIA_MOUSE  G1009000000001              4.3


    """
    #outdf = pd.DataFrame(columns = ['cid','score','cgid', 'ncorrect', 'nterms', 'method','aspect','f1max'])
    
    
    gb = evaldf.groupby('cid')
    
    outdf = gb.max().reset_index()
    droplist = ['goterm','correct','score']   
    outdf.drop(droplist, inplace=True, axis=1)
    
    #outdf['ncorrect'] = evaldf.groupby('cid')['correct'].transform('sum')
        
    ncdf = gb.apply(lambda x: (x.correct == True).sum()).reset_index(name='ncorrect')
    outdf['ncorrect'] = ncdf.ncorrect
    outdf.rename(columns={'f1': 'f1max'}, inplace=True)
    logging.debug(f"outdf:\n{outdf}") 
    return outdf


def run_phmmer(config, filename, version='current'):
    logging.debug(f"running phmmer on {filename}")
    (outfile, exclude_list, cidcgidmap) = execute_phmmer(config, filename, version)
    logging.debug(f"got outfile={outfile}")
    logging.debug(f"cidcgidmap is {cidcgidmap}")
    logging.debug(f"running parse_phmmer on {outfile} exclude_list={exclude_list}")
    phdict = parse_phmmer(config, outfile, exclude_list, cidcgidmap )
    logging.debug(f"got phdict...")
    return phdict

def get_phmmer_dict(config, filepath, version='current'):
    logging.info("running phmmer")
    phdict = run_phmmer(config, filepath, version)
    logging.debug(f"got phmmer dict length: {len(phdict)}")
    return phdict

def get_phmmer_df(config, filepath, version='current'):
    """
    orders by target,  evalue ascending (best fit first).  
    cache output for later usage....
    
    """
    logging.debug(f"Handling input file {filepath}")
    infile = os.path.expanduser(filepath)
    pcachedir = config.get('phmmer','cachedir') 

    filename = os.path.basename(os.path.expanduser(filepath))
    (filebase, e) = os.path.splitext(filename)
    pcachefile =f"{pcachedir}/{filebase}.phmmerdf.{version}.csv"
    
    df = None
    try:
        df = pd.read_csv(pcachefile, index_col=0)
        logging.info("Cached phmmer run found. ")
    except FileNotFoundError:
        logging.info("No cached phmmer data. Running...")      
        phd = get_phmmer_dict(config, filepath, version)
        if len(phd) > 0:
            df = pd.DataFrame.from_dict(phd, orient='index')
            df = df.sort_values(by=['cid','eval'], ascending=[True, True])
        if df is not None:
            logging.debug(f"Caching phmmer output to {pcachefile}")
            df.to_csv(pcachefile)
    return df


def execute_phmmer(config, filename, version='current'):
    """
    cpus = 8
    eval_threshold = 1.0e-120
    database=~/data/uniprot/uniprot_sprot.fasta
    remove_self_hits = True
    
    *Excludes geneids ( <protein>_<species> ) of sample from phmmer hit results*
    So the inbound sequence files *must* contain correct geneids (or some other
    method must be used to exclude self-hits). 

    """
    logging.debug(f"executing with filename={filename} version={version} ")
    exclude_list = [] 
    cidgidmap = get_tfa_geneids(filename)
    for c in cidgidmap.keys():
        exclude_list.append(c)     
    outdir = os.path.expanduser(config.get('phmmer','cachedir'))
    filename =os.path.expanduser(filename)
    outpath = os.path.dirname(filename)
    filebase = os.path.splitext(os.path.basename(filename))[0]
    outfile = "%s/%s.phmmer.tbl.txt" % (outdir, filebase)
    logging.debug(f"outfile={outfile}")
    cpus = config.get('phmmer','cpus')
    eval_threshold = config.get('phmmer','eval_threshold')
    database = os.path.expanduser(config.get('phmmer', 'database'))
    logging.debug(f"Using non-current version of uniprot for phmmer database: {database}")
    
    cmd = ["/usr/bin/which","phmmer"]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    res = p.stdout.readlines()
    logging.debug(f'which phmmer gave {res}')
    
    # cmd = f"phmmer --tblout {outfile} --noali --cpu {cpus} -E {eval_threshold} {filename} {database}"
    # args = f" --tblout {outfile} --noali --cpu {cpus} -E {eval_threshold} {filename} {database} "
    # cmd = ['phmmer', args ]
    cmd = [ 'phmmer',
           '--tblout', outfile , 
           '--noali',
           '--cpu', cpus,
           '-E', eval_threshold,
           filename,
           database 
           ]
    
    logging.debug(f"Running: {cmd}")
    cp = subprocess.run(cmd)
    #cp = subprocess.run(cmd, 
    #                    shell=False, 
    #                    universal_newlines=True, 
    #                    capture_output=True)
    
    logging.debug(f"Ran cmd='{cmd}' outfile={outfile} returncode={cp.returncode} " )
    #logging.debug(f"stdout='{cp.stdout}")
    #logging.debug(f"stderr='{cp.stderr}")
    logging.debug(f"returning outfile={outfile} exclude+list={exclude_list}")
    
    return (outfile, exclude_list, cidgidmap)


def get_tfa_geneids(filename):
    """
    Get the geneids of the TargetFiles in order to exclude these results from queries, e.g. 
    phmmer runs. This wouldn't necessarily cause a problem with the real files, because they 
    will be un-annotated. But for validation testing, the targets *will* be annotated (and
    need to be excluded). 
    
    Also need the mappings between cid and geneids for validation, where we need
    to look up the correct, annotated geneid to check predicitons. 
        
    """
    cidgids = {}
    
    try:
        f = open(filename, 'r')
    except FileNotFoundError:
        logging.error(f"file not readable {filename} ")
    for line in f:
        # >T100900000004 1433G_MOUSE
        if line.startswith(">"):
            fields = line[1:].split()
            cid = fields[0].strip()
            geneid = fields[1].strip()
            cidgids[cid] = geneid
    logging.debug(f"got {len(cidgids)} cids with geneids to exclude from results.")
    return cidgids
            

def parse_phmmer(config, filename, excludelist, cidcgidmap):
    """
    Read phmmer tbl out. Return dict  
    
    excludelist = [ '<geneids of inbound target proteins.>'] e.g. 2A5D_ARATH
     
    remove_self_hits = True|False
        shouldn't be necessary in standard case. 
    
    topx_threshold=<int>
        only include top X hits 
    
    """
    logging.info(f"Reading {filename}")
    df = pd.read_table(filename, 
                     names=['target','t-acc','cid','q-acc',
                            'eval', 'pscore', 'bias', 'e-value-dom','score-dom', 'bias-dom', 
                            'exp', 'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc', 'description'],
                     skip_blank_lines=True,
                     comment='#',
                     index_col=False,
                     skiprows=3,
                     engine='python', 
                     sep='\s+')
    logging.debug(f"got phmmer df:\n{df}")
    
    logging.debug("Dropping unneeded columns..")
    df = df.drop(['t-acc', 'q-acc','e-value-dom','score-dom', 'bias-dom', 'exp', 
             'reg', 'clu',  'ov', 'env', 'dom', 'rep', 'inc', 
             'description'] , axis=1)
           
    topx = config.getint('phmmer','topx_threshold')
    logging.debug(f"topx ={topx}")
    
    if topx is not None:
        df = df.groupby('cid').head(topx).reset_index(drop=True) 
    
    dict = df.to_dict(orient='index')
    logging.debug(f"dict is {dict}")
    idxtodel = []
    
    for idx in dict.keys():      
        (db,pacc,pid) = dict[idx]['target'].split('|')
        logging.debug(f"split |-separated target field...") 
        cid = dict[idx]['cid']
        # only exclude hits for the exact target protein...
        if pid == cidcgidmap[cid]:
            logging.debug(f"Found the pid {pid} excluded to be excluded for this cid {cid}")
            idxtodel.append(idx)
        else:
            (protein, species) = pid.split('_')
            dict[idx]['pacc'] = pacc
            dict[idx]['pid'] = pid
            dict[idx]['cgid'] = cidcgidmap[dict[idx]['cid']]
            del dict[idx]['target']
    for idx in idxtodel:
        del dict[idx]
    
    if logging.root.level >= logging.DEBUG: 
        s ="{"
        for k in dict.keys():
            s+=f"{k} : {dict[k]}\n"
        s += "}"
        logging.debug(f"{s}")
    """
    { 0 : {'cid': 'T100900000001', 'eval': 1.0999999999999997e-156, 'score': 523.6, 'bias': 8.5, 'pacc': 'Q9CQV8', 'pid': '1433B_MOUSE'}
      1 : {'cid': 'T100900000001', 'eval': 4.099999999999999e-155, 'score': 518.4, 'bias': 7.7, 'pacc': 'P35213', 'pid': '1433B_RAT'}
    }
    
    """
    logging.debug(f"returning OK. dict={dict}")
    return dict


def parse_expression_hd5(filename):
    """
    Loads data in file to dataframe.
    """
    with h5py.File(filename, 'r') as f:
        logging.debug("reading matrix...")
        matrix = f['agg'][:]
        logging.debug("reading rows. converting to unicode.")
        rows = [ s.decode() for s in  f['row'][:] ]
        logging.debug("reading columns. converting to unicode")
        columns = [ s.decode() for s in  f['col'][:] ]
        logging.debug("making dataframe...")
        df = pd.DataFrame(matrix,  index=rows, columns = columns )
    return df    

def get_expressionset_object(config):

    """
    Get handler to lazy load gene expression datasets...
    
    Sets are already implicitly cached in their datafiles, so no usecache needed. 
    
    """
    eobj = None
    if ExpressionSet.instance is not None:
        logging.debug("Found existing object.")
        eobj= ExpressionSet.instance
    else:
        logging.debug("No object, creating...")
        eobj = ExpressionSet(config)
        ExpressionSet.instance = eobj
    return eobj


def get_expression_dataset(config, species='YEAST'):
    """
    Assumes expression data filename format <taxid>_prioAggNet.hdf5
    Assumes mapping filename format <taxid>_gene_info.tab
    
    Current available:  
    ARATH BOVIN CAEEL CHICK DANRE DROME HUMAN MAIZE PIG RAT SOYBN YEAST   

    """
    df = None
    smo = get_specmap_object(config)
    basedir = os.path.expanduser(config.get('expression','datadir'))
    datasuffix = os.path.expanduser(config.get('expression', 'datasuffix'))
    prefix = smo.get_taxonid(species)
    
    filename = f"{basedir}/{prefix}{datasuffix}"
    logging.debug(f"Looking for expression data at: {filename}")
    if os.path.exists(filename):
        df = parse_expression_hd5( filename )
    else:
        df = None
    return df

    
def get_expression_genemap(config, species='YEAST'):
    """
    Assumes mapping filename format <taxid>_gene_info.tab
    
    Current available:  
    ARATH BOVIN CAEEL CHICK DANRE DROME HUMAN MAIZE PIG RAT SOYBN YEAST   

    """
    smo = get_specmap_object(config)
    basedir = os.path.expanduser(config.get('expression','datadir'))
    mapsuffix = os.path.expanduser(config.get('expression', 'mapsuffix'))
    prefix = smo.get_taxonid(species)
    
    filename = f"{basedir}/{prefix}{mapsuffix}"    
    logging.debug(f"Reading map table file {filename}")
    df = pd.read_table(filename)
    return df
    
    

def calc_orthoexpression_prediction(config, dataframe, usecache, version='2019'):
    """
    Takes phmmer df PDF: 
                cid           eval  pscore     bias   pacc          pid       cgid
    1     T100900000001  4.100000e-155  518.4   7.7  P35213    1433B_RAT    1A1L1_MOUSE
    2     T100900000001  5.400000e-155  518.0   7.2  A4K2U9    1433B_PONAB
    
    For each ortholog:
        parses species tag from pid:
        loads gene expression dataset for that species code (if it exists)
            gets ordered list of coexpressed UP paccs for ortholog.pacc
                take top <topx_threshold> [10] unless score is < <expscore_threshold> [.8 ]
                                
                sums the govectors (if they exist) for those UP paccs. 
                assigns those to the ortholog prediction. 
                score is the pscore/eval of the ortholog hit. 

    Outputs prediction:
   
                   cid      goterm         score         cgid
    0     G1009000000001  GO:0008150  7.164505e+63   CHIA_MOUSE
    1     G1009000000001  GO:0003674  2.687680e+63   CHIA_MOUSE
    2     G1009000000001  GO:0008152  2.687680e+63   CHIA_MOUSE
        
    
    """
    logging.debug(f"getting uniprot_by pacc.. version={version}")
    upbpacc = get_uniprot_by(config, by='pacc', version=version)
    ontobj = get_ontology_object(config, usecache)
    expdataobj = get_expressionset_object(config)
        
    gtlength = len(ontobj.gotermidx)
    topx_threshold = config.getint('orthoexpression', 'topx_threshold')
    expscore_threshold = config.getfloat('orthoexpression','expscore_threshold')
    score_method = config.get('orthoexpression' ,'score_method')
    max_goterms = config.getint('global','max_goterms')

    pdf = dataframe
    cidlist = list(pdf.cid.unique())
    logging.debug(f"cid list: {cidlist}")    
    gtarray = np.array(list(ontobj.gotermidx))
    
    # Dataframe to collect all calculated values. 
    topdf = pd.DataFrame(columns=['cid','goterm','score','cgid'])
    
    # Handle each cafa target in input list. 
    logging.debug("Starting to handle each CID...")
    for cid in cidlist:     
        cdf = pdf[pdf.cid == cid]
        cgid = cdf.reset_index().iloc[0].cgid       
        logging.debug(f"cgid for cid {cid} is {cgid}")
        logging.debug(f"one cid df:\n{cdf}")
        #               cid           eval  pscore  bias    pacc         pid        cgid
        # 1  G1009000000001  5.700000e-178   594.8   2.0  P04218    OX2G_RAT  OX2G_MOUSE
        # 2  G1009000000001  1.100000e-144   485.6   2.4  P41217  OX2G_HUMAN  OX2G_MOUSE
        
        # handle each ortholog protein, lookup expression data...
        # calculate score for inferred goterm set...       
        gvlist = []
        
        for (i, row) in cdf.iterrows():
            gv = np.zeros(gtlength)
            pscore = row.pscore
            (p, scode) = row.pid.split('_')
            logging.debug(f"got {p} and speciescode: {scode} ...")
            
            eds = expdataobj.get_dataset(scode)
            if eds is not None:
                logging.debug(f"Got dataset for species {scode}")   
                try:
                    expser = eds[row.pacc]
                    logging.debug(f"expser type is {type(expser)}")
                    logging.debug(f"Got expression info for pacc={row.pacc}")
                    if (len(expser.shape) > 1 ) and ( expser.shape[1] > 1 ):
                        logging.warning(f"Expression DF has too many columns!! sp={scode} col={row.pacc}" )
                        expser = expser.iloc[:,0]
                    
                    expser = expser.nlargest(n=topx_threshold)
                    epacclist = list(expser.index)
                    # list of accession codes: ['A0A024R410', 'P49411', ...]
                    for epacc in epacclist:
                        #gv = gv + upbpacc[prow.goterm].astype(np.int64)
                        try:
                            gv = gv + upbpacc[epacc]
                            logging.debug(f"Added GO vector for {epacc}")
                        except KeyError:
                            logging.debug(f"No GO vector for {epacc}")
                             
                except KeyError:
                    logging.debug(f"No expression info for pacc={row.pacc}")
            else:
                logging.debug(f"No expression data for species {scode}")
            
            if gv.sum() > 0:
                logging.debug(f"pscore={pscore}, applying...")
                if score_method == 'phmmer_score':
                    ones = np.ones(gtlength)
                    gv = np.minimum(gv, ones)
                    gv = gv * pscore
            
                if score_method == 'phmmer_score_weighted':
                    #logging.debug(f"gv.dtype={gv.dtype} max={gv.max()} pscore={pscore}")
                    gv = gv * pscore
                    #logging.debug(f"after gv.dtype={gv.dtype}") 
                gvlist.append(gv)
            
        if len(gvlist) > 0:    
            gv = np.zeros(gtlength)
            for v in gvlist:
                gv = gv + v            
            logging.debug(f"gv: {matrix_info(gv)}")
            gvnz = gv > 0
            gotermar = gtarray[gvnz]
            govalar = gv[gvnz]
            cidar = np.full( len(govalar), fill_value=cid)        
            cgidar = np.full(len(govalar), fill_value=cgid)
            df = pd.DataFrame({ 'cid': cidar, 
                               'goterm': gotermar, 
                               'score' : govalar, 
                               'cgid' : cgidar })
            
            logging.debug(f"dataframe is {df}")
            # This pulls out values, sorted by whatever 'score' is...
            df = df.nlargest(max_goterms, 'score')
            logging.debug(f"made dataframe for cid {cid}:\n{df}")
            topdf = topdf.append(df, ignore_index=True, sort=True)
        else:
            logging.info(f"No expression inference for {cid}. No predictions...")

    if len(topdf) > 0:
        logging.debug(f"made dataframe for all:\n{topdf}")
        topdf['method'] = 'orthoexpression'
        return topdf
    else:
        logging.debug(f"dataframe empty. returning None")
        return None

def calc_phmmer_prediction(config, dataframe, usecache, version='current'):
    """
    Takes phmmer df PDF: 
                cid           eval  pscore  bias    pacc          pid   cgid
    1     T100900000001  4.100000e-155  518.4   7.7  P35213    1433B_RAT   1A1L1_MOUSE
    2     T100900000001  5.400000e-155  518.0   7.2  A4K2U9  1433B_PONAB
    
    Gets uniprot_byterm_df UBTDF:
         pacc     species      goterm   goev    pid
    0        Q6GZX4   FRG3G     GO:0046782   IEA  001R_FRG3G
    1        Q6GZX3   FRG3G     GO:0033644   IEA  002L_FRG3G
    
    Algorithm:
        Determines predictive score by building goterm vector (term + parents) of all orthologs
        Score assigned is the score phmmer gave to hit. 
    
    Outputs prediction:
   
                   cid      goterm         score         cgid
    0     G1009000000001  GO:0008150  7.164505e+63   CHIA_MOUSE
    1     G1009000000001  GO:0003674  2.687680e+63   CHIA_MOUSE
    2     G1009000000001  GO:0008152  2.687680e+63   CHIA_MOUSE
   
    """
    logging.debug(f"inbound dataframe:\n{dataframe}")
    logging.debug(f"getting uniprot_byterm_df, version {version}")
    ubtdf = get_uniprot_byterm_df(config, usecache=usecache, version=version)
    ontobj = get_ontology_object(config, usecache)
    gtlength = len(ontobj.gotermidx)
    max_goterms = config.getint('global','max_goterms')
    score_method = config.get('phmmer','score_method')

    pdf = dataframe

    cidlist = list(pdf.cid.unique())
    logging.debug(f"cid list: {cidlist}")    
    gtarray = np.array(list(ontobj.gotermidx))
    
    # Dataframe to collect all calculated values. 
    topdf = pd.DataFrame(columns=['cid','goterm','score','cgid'])
    
    # Handle each cafa target in input list. 
    for cid in cidlist:
        #gv = np.zeros(gtlength)       
        
        cdf = pdf[pdf.cid == cid]
        cgid = cdf.reset_index().iloc[0].cgid       
        logging.debug(f"cgid for cid {cid} is {cgid}")
        logging.debug(f"one cid df:\n{cdf}")
        #               cid           eval  pscore  bias    pacc         pid        cgid
        # 1  G1009000000001  5.700000e-178   594.8   2.0  P04218    OX2G_RAT  OX2G_MOUSE
        # 2  G1009000000001  1.100000e-144   485.6   2.4  P41217  OX2G_HUMAN  OX2G_MOUSE
        
        # handle each ortholog protein, calculate score for inferred goterm set...
        for (i, row) in cdf.iterrows():
            gv = np.zeros(gtlength)       
            udf = ubtdf[ubtdf.pacc == row.pacc]
            pscore = row.pscore
            #           pacc   species   goterm     goev      pid
            # 2500682  O55201   MOUSE    GO:0032044  ISS  SPT5H_MOUSE
            # 2500683  O55201   MOUSE    GO:0005654  ISO  SPT5H_MOUSE
            # handle each goterm for this ortholog
            for (j, prow) in udf.iterrows():
                #try:
                newgv = ontobj[prow.goterm]
                #logging.debug(f"goterm for row is {prow.goterm} adding...")
                if newgv is not None:
                    logging.debug(f"Got GO vector for {prow.goterm}")
                    gv = gv + newgv.astype(np.int64)
                else:
                    logging.warning(f"No GO vector for {prow.goterm}")
               
            # we now have a govector with all goterms indicated by this ortholog.
            # each entry is sum of number of times that goterm was indicated (as annotated or
            # parent of an annotation term).
            # gv = array([123,   0,   3,   7, 345, 0])
        
        logging.debug(f"pscore={pscore} udf=\n{udf}")
        
        if score_method == 'phmmer_score':
            ones = np.ones(gtlength)
            gv = np.minimum(gv, ones)
            gv = gv * pscore
            #logging.debug(f"cdf is {cdf}")               
            #  cdf:
            #     cid           eval             pscore  bias    pacc          pid        cgid
            #74   G1009000000009  3.500000e-256  854.2  13.3  P33534   OPRK_MOUSE  OPRK_MOUSE
            #75   G1009000000009  5.900000e-254  846.9  13.3  P34975     OPRK_RAT  OPRK_MOUSE

        if score_method == 'phmmer_score_weighted':
            #logging.debug(f"gv.dtype={gv.dtype} max={gv.max()} pscore={pscore}")
            gv = gv * pscore

        gvnz = gv > 0
        #logging.debug(f"gvnz: {matrix_info(gvnz)}")
        #logging.debug(f"gtarray:length: {len(gtarray)} type:{type(gtarray)}")
        gotermar = gtarray[gvnz]
        govalar = gv[gvnz]
        cidar = np.full( len(govalar), fill_value=cid)        
        cgidar = np.full(len(govalar), fill_value=cgid)
        df = pd.DataFrame({'cid': cidar, 
                           'goterm': gotermar, 
                           'score' : govalar, 
                           'cgid' : cgidar })
        if len(df) > 0:
            logging.debug(f"dataframe is:\n{df}")
            # This pulls out values, sorted by whatever 'score' is...
            df = df.nlargest(max_goterms, 'score')
            logging.debug(f"made dataframe for cid {cid}:\n{df}")
            topdf = topdf.append(df, ignore_index=True, sort=True)
        else:
            logging.warning(f"No data for cid: {cid}")
        
    if len(topdf) > 0:
        topdf['method'] = 'phmmer'
        logging.debug(f"made dataframe for all:\n{topdf}")
        return topdf
    else:
        logging.warning("No data for input.")
        return None


def calc_expression_prediction(config, dataframe, usecache, version='2019'):
    """
    Takes generic df PDF: 
             cid       cgid    pacc
0  G982300000000   CD59_PIG  O62680
1  G982300000001  APOA1_PIG  P18648
2  G982300000002  ANKR1_PIG  Q865U8
    
    *NOTE* Assumes one and only one unique CGID per CID
    
    For each protein
        parses species tag from pid:
        loads gene expression dataset for that species code (if it exists)
        gets ordered list of coexpressed UP paccs for ortholog.pacc
        take top <topx_threshold> [100], filter by < <expscore_threshold> [.5 ]          
        for each coexpressed pacc:
            gets govector (if it exists)
            sums govector (value 1 per goterm)   


    Outputs prediction:
   
                   cid      goterm         score         cgid
    0     G1009000000001  GO:0008150  <num-seen>     CHIA_MOUSE
    1     G1009000000001  GO:0003674  <num-seen>    CHIA_MOUSE
    2     G1009000000001  GO:0008152  <num-seen>    CHIA_MOUSE
        
    
    """
    logging.debug(f"getting uniprot_by pacc.. version={version}")
    topx_threshold = config.getint('expression', 'topx_threshold')
    expscore_threshold = config.getfloat('expression','expscore_threshold')
    score_method = config.get('expression' ,'score_method')
    max_goterms = config.getint('global','max_goterms')

    upbpacc = get_uniprot_by(config, by='pacc', version=version)
    ontobj = get_ontology_object(config, usecache)
    gtlength = len(ontobj.gotermidx)
    gtarray = np.array(list(ontobj.gotermidx)) 
           
    expdataobj = get_expressionset_object(config)
    pdf = dataframe
    cidlist = list(pdf.cid.unique())
    logging.debug(f"cid list: {cidlist}")    

    # Dataframe to collect all calculated values. 
    topdf = pd.DataFrame(columns=['cid','goterm','score','cgid'])
    
    # Handle each cafa target in input list. In this case only one per row...
    #         cid         cgid     pacc
    # 0  G982300000000   CD59_PIG  O62680
    # 1  G982300000001  APOA1_PIG  P18648
    # 2  G982300000002  ANKR1_PIG  Q865U8
    #
    
    logging.debug("Starting to handle each CID...")
    for cid in cidlist:     
        cdf = pdf[pdf.cid == cid]
        cgid = cdf.reset_index().iloc[0].cgid       
        pacc = cdf.reset_index().iloc[0].pacc
        logging.debug(f"cgid for cid {cid} is {cgid}") 
        (p, scode) = cgid.split('_')
        logging.debug(f"got {p} and speciescode: {scode} ...")        
        eds = expdataobj.get_dataset(scode)
        
        gvlist = []
        gv = np.zeros(gtlength)
          
        if eds is not None:
            logging.debug(f"Got dataset for species {scode}")   
            try:
                expser = eds[pacc]
                logging.debug(f"expser type is {type(expser)}")
                logging.debug(f"Got expression info for pacc={pacc}")
                if (len(expser.shape) > 1 ) and ( expser.shape[1] > 1 ):
                    logging.warning(f"Expression DF has too many columns!! sp={scode} col={row.pacc}" )
                    expser = expser.iloc[:,0]
                
                expser = expser.nlargest(n=topx_threshold)
                logging.debug(f"expression series: expser:\n{expser}")
                #
                #  2020-04-07 11:58:08,054 (UTC) [ DEBUG ] fastcafa.py:1756 root.calc_expression_prediction(): expression series: expser:
                #  Q865U8        1.000000
                #  I3LAY6        0.642747
                #  A0A287A909    0.633572
                #  A0A287A8C2    0.630955
                #  A0A286ZT07    0.615943
                logging.debug(f"Removing values less than {expscore_threshold}")
                expser = expser[expser > expscore_threshold ]
                
                epacclist = list(expser.index)
                # list of accession codes: ['A0A024R410', 'P49411', ...]
                for epacc in epacclist:
                    try:                        
                        newgv = upbpacc[epacc].astype(int)
                        logging.debug(f"newgv is {newgv}")
                        gv = gv + newgv
                        logging.debug(f"Added GO vector for {epacc}")
                    
                    except KeyError:
                        logging.debug(f"No GO vector for {epacc}")

                logging.debug(f"{cid}: gv.sum={gv.sum()} gv.max={gv.max()}")
                
                gvnz = gv > 0
                gotermar = gtarray[gvnz]
                govalar = gv[gvnz]
                cidar = np.full( len(govalar), fill_value=cid)        
                cgidar = np.full(len(govalar), fill_value=cgid)
                df = pd.DataFrame({ 'cid': cidar, 
                                    'goterm': gotermar, 
                                    'score' : govalar, 
                                    'cgid' : cgidar })
                logging.debug(f"dataframe is {df}")
                df = df.nlargest(max_goterms, 'score')
                df = df.sort_values(by='score', ascending=False)
                topdf = topdf.append(df, ignore_index=True, sort=True)
                
            except KeyError:
                logging.debug(f"No expression info for pacc={pacc}")
        
        else:
            logging.debug(f"No expression data for species {scode}")
            
    topdf['method'] = 'expression'
    if len(topdf) > 0:
        logging.debug(f"made dataframe for all:\n{topdf}")
        return topdf
    else:
        logging.warning("No data for input...")
        return None


def check_filename_for_taxids(config, filename):
    """
    If a filename contains a valid CAFA taxid (e.g. 44689) between dots,
    returns species code, e.g. 'DICDI'
    else None 
    
    """
    logging.debug("Checking filename for species taxid...")
    t2s = get_cafaspecies(config)
    infile = os.path.expanduser(filename)
    filename = os.path.basename(infile)
    (filebase, ext) = os.path.splitext(filename)
    fields = filebase.split('.')
    ctids = list(t2s.keys())
    logging.debug(f"Looking for tids: {ctids} in {fields} ...")
    taxid = None
    found = False
    for ctid in ctids:
        for f in fields:
            if ctid == f:
                logging.debug(f"Found taxid {ctid}")
                return t2s[ctid]
    return None
    

def make_prior_prediction(config, infile, species=None):
    """
    Same as calc_phmmer_prediction, but assigns prior likelihoods as score
    Automatically detects species names/codes in filename, uses. 
    
    """
    tfalol = parse_tfa_file_lol(infile)
    logging.debug(f"Got list of lists length: {len(tfalol)}") 
    
    fnspecies = check_filename_for_taxids(config, infile)
    if species is None and fnspecies is not None:
        species = fnspecies
        logging.debug(f"Auto-identified species {species} from filename taxonid.")
    pdf = get_prior_df(config, True, species)
    logging.debug(f"Got prior frame:\n{pdf}")     
    
    plol = pdf.values.tolist()
    outlist = []
    
    for (cid, cgid) in tfalol:
        for (goterm, score) in plol:
            outlist.append([ cgid, cid, goterm, score])
    logging.debug(f"make list of lists length={len(outlist)}")
    outdf = pd.DataFrame(outlist, columns=['cgid','cid', 'goterm','score'] )
    outdf['method'] = 'prior'  
    return outdf


def get_goprior(config, usecache, species=None):
    gp = calc_prior(config, usecache, species)
    return gp

def get_altiddict(config, usecache):
    if ALTIDDICT is not None:
        return ALTIDDICT
    else:
        build_ontology(config, usecache)    
        return ALTIDDICT


def build_ontology(config, usecache):
    """
    obofile=~/data/go/go.obo
    cachedir = ~/play/cafa4      
    
    from parse_obo():
    { 'GO:2001315': 
         {'is_a': ['GO:0009226', 'GO:0046349', 'GO:2001313'], 
         'part_of': [], 
         'goterm': 'GO:2001315', 
         'goname': 'UDP-4-deoxy-4-formamido-beta-L-arabinopyranose biosynthetic process', 
         'goasp': 'bp', 
       ...  
    }
 
    result:  Numpy boolean matrix of all relations in ontology, ordered by sorted goterm name.  
       
    """
    #global GOMATRIX
    #global GOTERMIDX
    #global ALTIDDICT
    #global GOTERMLIST
    # def __init__(self, gomatrix, gotermidx, altidx):
    
    logging.debug(f"usecache={usecache}")
    cachedir = os.path.expanduser(config.get('ontology','cachedir'))
    ontologyfile = f"{cachedir}/ontology.npy"
    termidxfile = f"{cachedir}/gotermindex.pickle"
    altiddictfile = f"{cachedir}/altiddict.pickle"
    include_partof = config.getboolean('ontology','include_partof')
    
    gomatrix = None
    
    if os.path.exists(ontologyfile) and usecache:
        logging.debug("Cache hit. Using existing matrix...")
        gomatrix = np.load(ontologyfile)
        logging.debug(f"loaded matrix: {matrix_info(gomatrix)} from {ontologyfile}")
        
        f = open(termidxfile, 'rb')
        gotermidx = pickle.load(f)
        f.close()

        f = open(altiddictfile, 'rb')
        altiddict = pickle.load(f)
        f.close()
                
        logging.debug(f"goterm index, e.g. : \n{list(gotermidx)[0]} :  {gotermidx[list(gotermidx)[0]]} ")
    
    else:
        (godict, altiddict) = parse_obo(config)
        
        # get keys from dict
        gotermlist = list(godict)
        logging.debug(f"parsed obo with {len(gotermlist)} entries. ")
        logging.debug(f"example entry:\n{gotermlist[0]}")
        logging.debug("sorting goterms")
        gotermlist.sort()
        logging.debug(f"sorted: e.g. {gotermlist[0:5]} ")
        logging.debug("creating goterm index dict.")
        #
        i = 0
        gotermidx = {}
        for gt in gotermlist:
            gotermidx[gt] = i
            i = i + 1
              
        logging.debug(f"creating zero matrix of dimension {len(gotermlist)}")
        shape = (len(gotermlist), len(gotermlist))
        gomatrix = np.zeros( shape, dtype=bool )
        logging.debug(f"filling in parent matrix for all goterms...")
        for gt in godict.keys():
            for parent in godict[gt]['is_a']:
                    gomatrix[gotermidx[gt]][gotermidx[parent]] = True
        if include_partof:
            logging.debug("Including part_of relationships as is_a")
            for gt in godict.keys():
                for parent in godict[gt]['part_of']:
                        gomatrix[gotermidx[gt]][gotermidx[parent]] = True
        
        #logging.debug(f"initial matrix:\n{print_square(gomatrix, GOTERMLIST)}")
        logging.debug("Calculating sparsity...")
        logging.debug(f"sparsity = { 1.0 - np.count_nonzero(gomatrix) / gomatrix.size }")
        logging.debug("converting to sparse matrix.")
        gomatrix = sparse.lil_matrix(gomatrix, dtype=bool)
        logging.debug(f"converging matrix: {matrix_info(gomatrix)}")
        gomatrix = converge_sparse(gomatrix)
        logging.info(f"got converged matrix:\n{matrix_info(gomatrix)} ")
        logging.debug(f"converged matrix sum={gomatrix.sum()}")
        #logging.debug("Calculating sparsity...")
        #sparsity = 1.0 - np.count_nonzero(gomatrix) / gomatrix.size
        #logging.debug(f"sparsity = { 1.0 - np.count_nonzero(gomatrix) / gomatrix.size }")        
        gomatrix = gomatrix.todense()
        gomatrix = np.asarray(gomatrix, dtype='bool')
            
        logging.debug(f"Caching all values/indexes...")
        logging.debug(f"Saving matrix: {matrix_info(gomatrix)} to {ontologyfile}")
        np.save(ontologyfile, gomatrix)

        logging.debug(f"Saving gotermidx {len(gotermidx)} items to {termidxfile}")
        f = open(termidxfile, 'wb')   
        pickle.dump(gotermidx, f )
        f.close()
        
        logging.debug(f"Saving altiddict {len(altiddict)} items to {altiddictfile}.")
        f = open(altiddictfile, 'wb')
        pickle.dump(altiddict, f)
        f.close()
        
        logging.debug("Done constructing input for Ontology().")

    ontobj = Ontology(gomatrix, gotermidx, altiddict)
    # set global instance
    Ontology.instance = ontobj  
    logging.debug("Done creating Ontology object.")



def print_square(matrix, labels):
    """
    pretty prints square matrix with labels for debugging.
    """
    df = pd.DataFrame(matrix, columns=labels, index=labels)
    logging.debug(f"\n{df}")
    


def converge_sparse(matrix):
    logging.debug(f"starting matrix: \n{matrix_info(matrix)}")
    #logging.debug(f"{print_square(matrix.todense(), GOTERMLIST)}")
    oldval = 0
    logging.debug("Summing inbound matrix...")
    newval = matrix.sum()
    logging.debug("Beginning convergence loop.")
    icount = 0
    while oldval != newval:
        #logging.debug(f"Inbound matrix:\n{matrix_info(matrix)}")
        #logging.debug(f"oldval={oldval} newval={newval}")
        oldval = newval
        if not isinstance(matrix,  sparse.lil.lil_matrix): 
            logging.debug(f"{type(matrix)} is not scipy.sparse.lil.lil_matrix, converting... ")
            matrix = sparse.lil_matrix(matrix, dtype=np.bool)
        else:
            pass
            #logging.debug("matrix already lil_matrix...")
        #logging.debug("Multiplying...")
        mult = matrix @ matrix
        #logging.debug("Adding back original...")
        matrix = mult + matrix
        #logging.debug("Getting new sum...")
        newval = matrix.sum()
        #logging.debug(f"New matrix {icount}:\n{matrix.todense()}")
        #logging.debug(f"{print_square(matrix.todense(), GOTERMLIST)}")
        icount += 1
    logging.debug(f"Convergence complete. {icount} iterations. matrix:\n{matrix_info(matrix)}")
    return matrix





def get_prior_df(config, usecache=True, species = None):

    ontobj = get_ontology_object(config, usecache=True)
    priormatrix = calc_prior(config, usecache, species)    
    priormatrix = priormatrix.transpose()
    logging.debug(f"priormatrix shape: {priormatrix.shape}")
    df = pd.DataFrame(priormatrix, index=ontobj.gotermlist, columns=['score'])
    df.reset_index(inplace=True)
    map = {'index' : 'goterm'}
    df.rename(axis='columns', mapper=map, inplace=True)
    df = df.sort_values(by='score', ascending=False)
    df.reset_index(drop=True, inplace=True)
    prevlength = df.shape[0]
    df = df[df.score != 0]
    newlength = df.shape[0]
    if newlength != prevlength:
        logging.debug(f"Reduced prior DF length by removing { prevlength - newlength } zero values.")
    return df
   
    
def calc_prior(config, usecache, species=None, version='current'):
    """
    take ontology matrix. goterms x goterms
    make total_termvector[ 17k ]
    
       0 proteinacc   1 species  2 goterm       3 goevidence
    0  '3AHDP'       'RUMGV'    'GO:0016491'    'IEA'
    1  '3AHDP'       'RUMGV'    'GO:0006694'    'TAS'

    for each protein:
        for each protein.goterm:
            gtvector = get_vector(goterm)
            total_termvector += gtvector
    
    if species is specified (as species code, e.g CAEEL), calculate prior for species only, 
    otherwise global
    
    ->
    
    vector of prob (.0001 - .99999)  [ 0.001, .0020, ... ] indexed by sorted gotermlist. 

    """ 

    logging.debug(f"Called with species={species}")
    freqarray = None
    fspec = None  # species code for filename.
    
    smo = get_specmap_object(config)
    
    if species is None:
        fspec = 'all'
    else:
        fspec = smo.get_taxonid(species)

    cachedir = os.path.expanduser(config.get('uniprot','cachedir'))
    cachefile = f"{cachedir}/uniprot.goprior.{fspec}.npy"       
    
    if usecache and os.path.exists(cachefile):
        freqarray = np.load(cachefile)
        logging.debug(f"Loaded prior freqarray from file: {cachefile}")
    
    else:
        ontobj = get_ontology_object(config, usecache=True)
        gtlength = len(ontobj.gotermlist)
        logging.debug(f"building sprot. species={species}")
        sprot = get_uniprot_byterm(config, usecache=True, version=version)
        sdf = pd.DataFrame(sprot,columns=['pacc', 'pid', 'protein', 
                                          'species', 'goterm','goasp','goev', 
                                          'seqlen', 'seq', 'gene'                                        
                                          ])
        logging.debug(f"Built dataframe:\n{sdf}")    
        
        # evaluation is done only on experimentally validated annotations. 
        # so calculate prior solely using experimentally validated annotations. 
        sdf = filter_goevidence(config, sdf)
        sprot = sdf.values.tolist()
        
        if species is not None: 
            logging.debug(f"species {species} specified. converted to df with {sdf.shape[0]} rows:\n{sdf}")
            sdf = sdf[sdf.species == species ]
            logging.debug(f"removed other species. {sdf.shape[0]} rows left.")
                             
        sumarray = np.zeros(gtlength, dtype=np.int)
        logging.debug(f"made array {sumarray}")
        
        i = 0
        elist = []
        logging.debug("summing over all terms and parents...")
        for r in sprot:
            gt = r[4]
            gv = ontobj[gt]
            if gv is not None:
                sumarray = sumarray + ontobj[gt]
                elist.append(gt)
            i += 1
        eset = set(elist)
        #logging.debug(f"missing keys: {list(eset)}")
        logging.debug(f"added {i} gomatrix lines, with {len(eset)} distinct missing keys. ")
        logging.debug(f"sumarray: {sumarray} max={sumarray.max()} min={sumarray.min()} dtype={sumarray.dtype}")    
        divisor = sumarray.sum()
        logging.debug(f"divisor={divisor}")
        freqarray = sumarray / divisor
        np.save(cachefile, freqarray)
        logging.debug(f"Saved freqarray to {cachefile}")

    logging.debug(f"freqarray: {freqarray.dtype} {freqarray} max={freqarray.max()} min={freqarray.min()}")    
    return freqarray


def get_uniprot_byterm_df(config, usecache=True, version='2019', goaspect=None):
    logging.debug(f"called with usecache={usecache} version={version}")
    lol = get_uniprot_byterm(config, usecache=usecache, version=version, goaspect=goaspect)
    df = pd.DataFrame(lol,columns=['pacc', 'pid', 'protein', 'species', 
                                      'goterm','goasp','goev','seqlen','seq','gene'])  
    logging.debug(f"Built dataframe:\n{df}")
    return df 


def get_uniprot_byterm(config, usecache, version='2019', goaspect=None):
    """
    [ {'proteinid': '001R_FRG3G', 
       'protein': '001R', 
       'species': 'FRG3G', 
       'proteinacc': 'Q6GZX4', 
       'taxonid': '654924', 
       'goterms': {'GO:0046782': ['cc','IEA'}, 
       'seqlength': 256, 
       'sequence': 'MAFSAEDVL......LYDDSFRKIYTDLGWKFTPL'},
       'gene' : '128UP',
       .
       .
       .
    ]
    
    Generate non-redundant uniprot with a row for every goterm:
        pacc species      goterm goasp goev          pid   gene
0     P0DJZ0   PAVHV  GO:0030430  bp    IDA    11K_PAVHV    11K
1     P32234   DROME  GO:0005525  mf    IDA  128UP_DROME  128UP
2     P83011   SCYCA  GO:0043231  cc    IDA  13KDA_SCYCA     {} 
   
       for p in lod:
            newgts = {}
            for gt in p['goterms'].keys():
                evcode = p['goterms'][gt]
                item = [ p['proteinacc'],
                         p['protein'],
                         p['species'],
                         gt, 
                         evcode,
                         p['seqlength'],
                         p['sequence'] 
                      ]
                lodt.append(item)
    """    
    logging.debug(f"calling build_uniprot usecache={usecache} version={version}")
    lod = build_uniprot(config, usecache=True, version=version)
    ubt = []
    for p in lod:
        #
        # many entries do *not* have a gene symbol assigned, and the defaultdict 
        # implementation gives back and empty dict. Fix this....
        e = p['gene']
        if len(e) == 0:
            p['gene'] = np.nan
        # and proceed...
        if len(p['goterms'].keys()) > 0:           
            for gt in p['goterms'].keys():
                evcode = p['goterms'][gt][1]
                aspcode = p['goterms'][gt][0]
                
                if goaspect is None or aspcode == goaspect:
                    #logging.debug(f"goaspect={goaspect} aspcode={aspcode} adding...")
                    item = [ p['proteinacc'],
                             p['proteinid'],
                             p['protein'],
                             p['species'],
                             gt,
                             aspcode,
                             evcode, 
                             p['seqlength'],
                             p['sequence'], 
                             p['gene'],
                         ]
                    ubt.append(item)
                    
        else:
            item = [ p['proteinacc'],
                         p['proteinid'],
                         p['protein'],
                         p['species'],
                         np.nan,
                         np.nan,
                         np.nan, 
                         p['seqlength'],
                         p['sequence'], 
                         p['gene'],
                     ]
            ubt.append(item)
    logging.debug(f"created uniprot_byterm with {len(ubt)} entries.")          
    return ubt  
        

def get_uniprot_testset(config, usecache, species, evidence=[ 'EXP', 'IDA', 'IMP', 'IGI', 'IEP' ]):
    """
    species  = 'MOUSE' 'HUMAN' -> taxonid matched. 
    evidence = list of OK codes | None means any
        Experimental: [ 'EXP', 'IDA', 'IMP', 'IGI', 'IEP' ] q
    
    Generate DataFrame with data for export as CAFA test file:
    G<taxonid>_<datetime>_<5-digit-number>    <proteinid>   <sequence>
     
 ['Q6GZX4', 'FRG3G', 'GO:0046782', 'IEA', 'IEA', '001R_FRG3G', 256, 'MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPSEKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLDAKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHLEKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDDSFRKIYTDLGWKFTPL', nan]
      
    """
    #lodt = build_uniprot_test(config, usecache)
    #lodt = build_uniprot(config, usecache)
    lodt = get_uniprot_byterm(config, usecache)
    logging.debug(f"lodt: \n{lodt[0]}")

    #tdf = pd.DataFrame(lodt, columns=['pacc','protein','species', 'goterm','goev','seqlen','seq'])
    tdf = pd.DataFrame(lodt, columns=['pacc', 'pid', 'protein', 'species', 
                                      'goterm','goev','seqlen','seq','gene'])
    logging.debug(f"testing tdf is:\n{tdf}")
    if evidence is not None:
        tdf = tdf[tdf['goev'].isin(evidence)] 
    tdf = tdf[tdf['species'] == species]
    
    tdf.reset_index(inplace=True)
    un = tdf.protein.unique()
    evcodes = tdf.goev.unique()
    logging.debug(f"generated {len(tdf)} item DF with {len(un)} proteins. evcode={evcodes}")
    return tdf



def do_testset_old(config, numseq, species, outfile):
    """
    Creates random target set from *annotated* swissprot proteins for given species. 
    Creates FASTA files exactly like CAFA TargetFiles
        
    :param    type      name:   desc
    :param    str       species          Species to generate: 
    :return   testfile  Path to test file generated in fasta format.
    :rtype
    :raise:
    
    """
    outfile = os.path.expanduser(outfile)
    numseq = int(numseq)
    logging.debug(f"numseq={numseq} species={species} outfile={outfile}")
    tdf = get_uniprot_testset(config, usecache=True, species=species)

    logging.info(f"got testset dataframe {len(tdf)} entries. ")
    
    #  [  { taxid     : speccode, ...  },
    #     { speccode  : taxonid, ...},
    #     { linnean   : taxonid, ...}   ]
    specmaps = get_specmaps(config)
    taxonid = specmaps[1][species.strip()]
    logging.debug(f"taxonid is {taxonid} for {species}")    
    up = tdf.protein.unique()
    upl = up.tolist()

    # throws error if 
    spl = random.sample(upl, numseq)
    snum = 1
    x = 60
    s = ""
    for p in spl:
        r = tdf[tdf.protein == p].reset_index().iloc[0]
        s += f">G{taxonid}{snum:08} {r.protein}_{r.species}\n"
        chunklist = [ r.seq[y-x:y] for y in range(x, len(r.seq)+x, x) ] 
        #logging.debug(f"chunklist {chunklist}")
        for c in chunklist:
            s += f"{c}\n"

        snum += 1
    #logging.debug(s)
    
    try:
        f = open(outfile, 'w')
        f.write(s)
        f.close()
        logging.debug(f"Wrote test data to file {outfile}")
    except IOError:
        logging.error(f"could not write to file {outfile}")
        traceback.print_exc(file=sys.stdout) 

    return outfile



def build_prior(config, usecache, species=None, outfile=None, version='current'):
    """
    
    Additionally, if outfile is specified, writes CSV to file of ranked 
    goterms with frequency as pest, similar to prior prediction but without input. 
    
    """
    
    logging.debug(f"running calc_prior with species {species}")
    out = calc_prior(config, usecache, species, version)
    if outfile is not None:
        outfile = os.path.expanduser(outfile)
        ontobj = get_ontology_object(config, usecache=True)
        df = pd.DataFrame(out, columns=['score'])
        df['goterm'] = ontobj.gotermlist
        df.sort_values(by='score', ascending=False, inplace=True)
        logging.debug(f"prior has {df.shape[0]} rows. ")
        df = df[df.score != 0.0]
        logging.debug(f"prior has {df.shape[0]} non-zero rows. ")
        df.reset_index(drop=True, inplace=True)
        df.to_csv(outfile)
        logging.debug(f"Wrote prior df to outfile {outfile}:\n{df}")
    return out


def get_uniprot_df(config, usecache=True, version='2019'):
    logging.debug(f"called with usecache={usecache}")
    lod = build_uniprot(config, usecache, version)
    df = pd.DataFrame(lod)  
    logging.debug(f"Built dataframe:\n{df}")    
    return df  


def build_uniprot(config, usecache, version='2019', goaspect=None):
    """
    Lowest-level uniprot build. calls parsing code .dat file. 
    Builds list of dictionaries, each element is item in uniprot/sprot
    
    """
    logging.debug(f"getting uniprot usecache={usecache} version={version}")
    cachedir = os.path.expanduser(config.get('uniprot','cachedir'))    
    evcode_flags = config.getboolean('uniprot','evcode_flags')
    cachefile = f"{cachedir}/uniprot.{version}.pickle"    
    listofdicts = None
    if os.path.exists(cachefile) and usecache:
        logging.debug("Cache hit. Using existing info...")    
        try:
            cf = open(cachefile, 'rb')    
            listofdicts = pickle.load(cf)
        except Exception:
            logging.error(f'unable to load via pickle from {cachefile}')
            traceback.print_exc(file=sys.stdout)    
        finally:
            cf.close()       
    else:
        listofdicts = parse_uniprot_dat(config, version )
        logging.debug(f"evcode_flags is {evcode_flags}")
        if evcode_flags:
            logging.debug(f"adding boolean evidence code flags to each dict entry.")
            listofdicts = add_evidence_flags(config, listofdicts)
            logging.debug(f"got adjusted list of dicts...")
        
        logging.debug(f"saving listofdicts: to {cachefile}")
        try:
            cf = open(cachefile, 'wb')    
            pickle.dump(listofdicts, cf )
            logging.debug(f"saved listofdicts: to {cachefile}")
        except Exception as e:
            logging.error(f'unable to dump via pickle to {cachefile}')
            traceback.print_exc(file=sys.stdout)     
        finally:
            cf.close()        
    return listofdicts

def add_evidence_flags(config, lod):
    """
    Adds boolean flag entries for evidence code categories:
    exp_codes = EXP, IDA, IMP, IGI, IEP, IPI 
    ele_codes = IEA
    cur_codes = ISS, NAS, TAS, IC, IRD, IGC, IBD, IBA, ISO, ISA, ISM, RCA, IKR
    
    is_exp   has experimental annotation
    is_ele   has IEA/electronic annotation
    is_cur   has curated non-experimental annotation
    is_ann   has any annotation. 
    
    'goterms': {'GO:0000334': 'IDA',
               'GO:0055114': 'IEA',
               'GO:0019363': 'IEA'},

    """
    logging.debug("adding evidence code boolean flags to uniprot data...")
    exp_codes = [ x.strip() for x in config.get('uniprot','exp_codes').split(',')]  
    ele_codes = [ x.strip() for x in config.get('uniprot','ele_codes').split(',')]
    cur_codes = [ x.strip() for x in config.get('uniprot','cur_codes').split(',')]
    
    for d in lod:
        #logging.debug(f"dict: {d}")
        try:
            gts = d['goterms']
            if len(gts) > 0:
                d['is_ann'] = True
                d['is_exp'] = False
                d['is_ele'] = False
                d['is_cur'] = False
                for gt in gts.keys():
                    ec = gts[gt][1]
                    if ec in exp_codes:
                        d['is_exp'] = True
                    elif ec in ele_codes:
                        d['is_ele'] = True
                    elif ec in cur_codes:
                        d['is_cur'] = True
                    else:
                        logging.warning(f"annotated protein but ecode='{ec}' not account for")
            else:
                d['is_ann'] = False                    
        except KeyError:
            d['is_ann'] = False

    return lod



def do_testset(config, numseq=10, species=None, outfile=None, limited=False, previous='2017',current='2019'):
    logging.debug(f"numseq={numseq} species={species} outfile={outfile} limited={limited}")
    logging.debug("Determining newly annotated proteins...")
    newannot = get_newannotated_df(config, limited, previous, current)
    newannot.drop_duplicates(['pacc'], inplace=True, keep='first')
    newannot.reset_index(drop=True, inplace=True)
    logging.debug(f"got {len(newannot)} unique newly-annotated entries with limited={limited}")
    bm_species = [ x.strip() for x in config.get('testset','cafa_species').split(',')]
    outdir = os.path.expanduser(config.get('global','outdir'))
    #logging.debug(f"species are {bm_species}")
    do_species = None
    if species is None:
        do_species = bm_species
    else:
        do_species = [species]    
    logging.debug(f"species are {do_species}")
    
    if limited:
        evtag = 'limited'
    else:
        evtag = 'noknow'
    
    for sp in do_species:
        logging.debug(f"Handling {sp}")
        specmap = get_specmap_object(config)
        taxonid = specmap.get_taxonid(sp)
        logging.debug(f"taxonid is {taxonid} for {sp}")    
        # keep just one row per accession num. no need for dupes. 
        naspec = newannot[newannot.species == sp] 
        upl = naspec.pacc.tolist()
        numavail = len(upl)
        logging.debug(f"there are {numavail} sequences available. ")
        if numavail < numseq:
            spl = upl
        else:
            spl = random.sample(upl, numseq)
        if numavail > 0:
            logging.debug(f"There are {numavail} available.")
            snum = 0
            x = 60
            s = ""
            for p in spl:
                r = naspec[naspec.pacc == p].reset_index().iloc[0]
                s += f">G{taxonid}{snum:08} {r.protein}_{r.species}\n"
                chunklist = [ r.seq[y-x:y] for y in range(x, len(r.seq)+x, x) ] 
                for c in chunklist:
                    s += f"{c}\n"
                snum += 1
            #logging.debug(s)
            if outfile is None:
                writefile = f"{outdir}/sp_species.{sp}.bmtest.{evtag}.{snum}.tfa"
            else:
                writefile = outfile
            write_tfa(writefile, s)
        else:
            logging.info(f"No newly annotated sequences for species {sp}")


def write_tfa(outfile, text):        
        try:
            f = open(outfile, 'w')
            f.write(text)
            logging.debug(f"Wrote TFA sequence to file {outfile}")
        except IOError:
            logging.error(f"could not write to file {outfile}")
            traceback.print_exc(file=sys.stdout) 
        finally:
            f.close()
            
  
    
    
def get_newannotated_df(config , limited=False, previous='2017', current='2019' ):
    """
    pulls out only entries newly annotated between 'previous' and 'current' versions
    of uniprot (as defined by config). 
    
    limited=True means previous may be annotated with non-experimental evidence. 
    limited=False means previous must be completely un-annotated.     
    
    In both cases current must have experimental evidence. 
    
    """
    nonexp_goev = ['IEA', 'ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'RCA', 'IBA', 'IBD', 'IKR', 'IRD']
    exp_goev=[ 'EXP', 'IDA', 'IMP', 'IGI', 'IEP' ]
    
    # get previous uniprot any evidence type 
    uprev = get_uniprot_byterm_df(config, usecache=True, version=previous)
    uprevunannot = uprev[uprev.goterm.isna()]
    uprevunannotp = uprevunannot['pacc'].unique()
    logging.debug(f"{len(uprevunannotp)} unique un-annotated proteins.")
    
    # start with annotated...
    uprevannot = uprev[uprev.goterm.notna()]
    if limited == False:
        # just use completely unannotated (no-knowledge)
        uprev = uprevunannot
        logging.debug(f"previously completely un-annotated:\n{uprev}")
    else:
        # electronic OK, must remove all proteins with *any* experimental goev
        # have unannotated. need non-experimental
        # get unannotated + non-experimental annotated (limited-knowledge)
        expannot = uprev[uprev['goev'].isin(exp_goev)] 
        expannotpset = set(expannot['pacc'].unique())
        logging.debug(f" got {len(expannotpset)} unique electronically annotated proteins.")
        
        eleannot = uprev[uprev['goev'].isin(nonexp_goev)]
        eleannotpset = set(eleannot['pacc'].unique())
        logging.debug(f" got {len(eleannotpset)} unique electronically annotated proteins.")

        eleconlypset = eleannotpset.difference(expannotpset)
        eleconlyp = list(eleconlypset)
        logging.debug(f" got {len(eleconlyp)} unique electronic-only annotated proteins.")
        eleconlydf = pd.DataFrame(eleconlyp, columns=['pacc'])
        logging.debug(f"made 1-column dataframe:\n{eleconlydf}")
        elecprev = pd.merge(eleconlydf, uprev, how='inner',on=['pacc','pacc'])
        elecprev.reset_index(drop=True, inplace=True)
        logging.debug(f"got filtered electronic-only previous DF:\n{elecprev}")
        uprev = elecprev.append(uprevunannot)
        logging.debug(f"got filtered previous DF:\n{uprev}")    
        #dataframe.reset_index(drop=True, inplace=True) 
        
    # get current uniprot, experimental only.
    unow = get_uniprot_byterm_df(config, usecache=True, version=current) 
    unow = unow[unow.goterm.notna()]
    logging.debug(f"current experimentally annotated:\n{unow}")
    
    unannot = uprev['pacc'].sort_values()
    unannot = pd.Series(unannot.unique())
    unannot.reset_index(drop=True, inplace=True) 
    
    nowannot = unow['pacc'].sort_values()      
    nowannot = pd.Series(nowannot.unique())
    nowannot.reset_index(drop=True, inplace=True)

    # get intersection of two series as DF
    newlydf = pd.DataFrame(pd.Series(np.intersect1d(unannot, nowannot)), columns=['pacc'])
    # provide full dataframe for just intersection...
    newlyannot = pd.merge(newlydf, unow , how='inner', on=['pacc','pacc'])
    return newlyannot



def build_specmaps(config, usecache):
    """
    builds three maps in form of list of dicts:
       [
         { taxid     : speccode, ...  },
         { speccode  : taxonid, ...},
         { linnean   : taxonid, ...}
         { taxonid   : linnean, ...}
       ]
    caches as pickle object. 
    
    """       
    logging.debug(f"usecache={usecache}")
    cachedir = os.path.expanduser(config.get('uniprot','cachedir'))
    cachefile = f"{cachedir}/specfile.pickle"
    specfile = os.path.expanduser(config.get('uniprot','specfile'))
    #     tax2spec  spec2tax  lin2tax
    map = [ {}, {}, {}, {} ] 
    if os.path.exists(cachefile) and usecache:
        logging.debug("Cache hit. Using existing info...")
        try:
            cf = open(cachefile, 'rb')    
            map = pickle.load(cf)
        except Exception:
            logging.error(f'unable to load via pickle from {cachefile}')
            traceback.print_exc(file=sys.stdout)    
        finally:
            cf.close()
        
    else:
        entries = parse_speclist(config, specfile )
        # [ ('CALCT', 'E', '227173', 'Calidris canutus', 'Red knot'), ... ]
        for entry in entries:
            map[0][entry[2]] = entry[0]
            map[1][entry[0]] = entry[2]
            map[2][entry[3]] = entry[2]
            map[3][entry[2]] = entry[3]
                    
        logging.debug(f"saving map: to {cachefile}")
        try:
            cf = open(cachefile, 'wb')    
            pickle.dump(map, cf )
            
        except Exception as e:
            logging.error(f'unable to dump via pickle to {cachefile}')
            traceback.print_exc(file=sys.stdout)     
        finally:
            cf.close()    
    logging.debug(f"saving map: to {cachefile}")
    return map


def build_idmaps(config, version='2019'):
    """
    Builds O(1) lookup maps of identifiers:
    
    pacc   pid
    idmaps = { 
               'pid2pacc' : { },  
               'pacc2pid' : {}}
            }

    """
    lod = build_uniprot(config, usecache=True, version=version )   
    idmaps ={}
    pid2pacc = {}
    pacc2pid = {}
    
    for item in lod:
        pid2pacc[item['proteinid']] = item['proteinacc']
        pacc2pid[item['proteinacc']] = item['proteinid']

    idmaps['pid2pacc'] = pid2pacc    
    idmaps['pacc2pid'] = pacc2pid
    
    return idmaps

def get_idmaps(config, usecache=True, version='2019'):
    
    
    cachedir = os.path.expanduser(config.get('uniprot','cachedir'))
    cachefile = f"{cachedir}/idmaps.{version}.pickle"

    if os.path.exists(cachefile) and usecache:
        logging.debug("Cache hit. Using existing info...")    
        try:
            cf = open(cachefile, 'rb')    
            idmaps = pickle.load(cf)
        except Exception:
            logging.error(f'unable to load via pickle from {cachefile}')
            traceback.print_exc(file=sys.stdout)    
        finally:
            cf.close()       
    else:
        idmaps = build_idmaps(config, version)
        logging.debug(f"saving dict: to {cachefile}")
        try:
            cf = open(cachefile, 'wb')    
            pickle.dump(idmaps, cf )
            logging.debug(f"saved dict: to {cachefile}")
        except Exception as e:
            logging.error(f'unable to dump via pickle to {cachefile}')
            traceback.print_exc(file=sys.stdout)     
        finally:
            cf.close()        
    return idmaps
    
    

def get_cafaspecies(config):
    """
    Creates dict for cafa species:
    taxid -> scode
    
    """
    smo = get_specmap_object(config)
    speclist = [x.strip() for x in config.get('testset','cafa_species').split(',') ]
    t2s = {}
    for scode in speclist:
        t = smo.get_taxonid(scode)
        t2s[t] = scode
    return t2s


def get_specmaps(config):
    lod = build_specmaps(config, usecache=True)
    return lod    

def get_specmap_object(config):
    sm = SpeciesMap.instance
    if SpeciesMap.instance is None:
        lod = build_specmaps(config, usecache=True)
        sm = SpeciesMap(lod)
        SpeciesMap.instance = sm
    return sm


def parse_uniprot_dat(config, version='2019'):
        """
        Parses uniprot/sprot DAT file, returns as list of dicts, with sub-dicts...
        [ {'proteinid': '4CLL9_ARATH', 
           'protein': '4CLL9', 
           'species': 'ARATH', 
           'proteinacc': 'Q84P23', 
           'taxonid': '3702', 
           'goterms': {'GO:0005777': ['cc', 'IDA'], 
                       'GO:0005524': ['mf', 'IEA'], 
                       'GO:0004321': ['mf', 'IDA'], 
                       'GO:0016874': ['mf', 'IEA'], 
                       'GO:0009695': ['bp', 'IDA'], 
                       'GO:0031408': ['bp', 'IEA']}
           },
           .
           .
           .
    
        """
        logging.debug(f"using config with sections: {config.sections()}...")
        filepath = os.path.expanduser(config.get('uniprot', version))
        logging.debug(f"opening datfile={filepath}")
        try:
            logging.debug(f" attempting to open '{filepath}'")
            filehandle = open(filepath, 'r')
        except FileNotFoundError:
            logging.error(f"No such file {filepath}")                
        
        allentries = []
        current = None
        sumreport = 1
        suminterval = 10000
        repthresh = sumreport * suminterval
        try:
            while True:
                line = filehandle.readline()
                if line == '':
                    break

                if line.startswith("ID "):
                    # ID   001R_FRG3G              Reviewed;         256 AA.
                    #      <prot_name>_<prot_spec>
                    proteinid = line[5:16].strip()
                    current = defaultdict(dict)
                    current['proteinid'] = proteinid
                    (protein, species) = proteinid.split('_')
                    current['protein'] = protein
                    current['species'] = species
                    #logging.debug("Handling ID. New entry.")                
                
                elif line.startswith("AC "):
                    # AC   Q6GZX4;
                    # AC   A0A023GPJ0; 
                    # AC   Q91896; O57469;
                    #logging.debug("Handling AC.")
                    rest = line[5:]
                    accession = rest.split(';')[0]
                    #accession = line[5:11].strip()
                    current['proteinacc'] = accession

                elif line.startswith("OX   "):
                    #OX   NCBI_TaxID=654924;
                    #logging.debug("Handling OX.")
                    taxonid = ""
                    val = line[5:]
                    fields = val.split('=')
                    if fields[0] == 'NCBI_TaxID':
                        taxonid = fields[1].strip().replace(';','')
                    current['taxonid'] = taxonid
                    
                elif line.startswith("DR   GO;"):
                    # DR   GO; GO:0046782; P:regulation of viral transcription; IEA:InterPro.
                    # P biological process, C cellular component, F molecular function.  
                    #logging.debug("Handling DR.")
                    fields = line.split(';')
                    goterm = fields[1].strip()
                    goinfo = fields[2]
                    aspcode = goinfo.split(':')[0].strip()
                    goaspect = UPASPECTMAP[aspcode]
                    goevsrc = fields[3]
                    (goevidence, evsrc) = goevsrc.split(':') 
                    goevidence = goevidence.strip()
                    current['goterms'][goterm] = [ goaspect,  goevidence ]

                elif line.startswith("SQ   SEQUENCE"):
                    #logging.debug("Handling SQ:  XXX")
                    # line = filehandle.readline()
                    current['seqlength'] = int(line.split()[2])
                    current['sequence'] = ""
                    seqlen = current['seqlength']
                    aaread = 0
                    while aaread < seqlen:
                        line = filehandle.readline()
                        lineseq = line.strip().replace(" ","")
                        current['sequence'] = "%s%s" % (current['sequence'], lineseq)
                        aaread += len(lineseq) 

                elif line.startswith("GN   "):
                    # Examples:
                    #  GN   ABL1 {ECO:0000303|PubMed:21546455},
                    #  GN   Name=BRCA1; Synonyms=RNF53;
                    #  GN   ORFNames=T13E15.24/T13E15.23, T14P1.25/T14P1.24;
                    # logging.debug(f"Handling GN. {line}")
                    val = line[5:]
                    if val.startswith("Name="):
                        fields = val.split()   # by whitespace
                        (n, gname) = fields[0].split("=")
                        gname = gname.upper()
                        #logging.debug(f"Gene name is {gname} ")
                        current['gene'] = gname.replace(';','') 
            
                elif line.startswith("//"):
                    #logging.debug("End of entry.")                  
                    allentries.append(current)
                    #logging.debug(f"All entries list now {len(allentries)} items... ")
                    if len(allentries) >= repthresh:
                        logging.info(f"Processed {len(allentries)} entries... ")
                        sumreport +=1
                        repthresh = sumreport * suminterval
                    current = None
                
        except Exception as e:
            traceback.print_exc(file=sys.stdout)                
        
        if filehandle is not None:
            filehandle.close()      
        
        logging.info(f"Parsed file with {len(allentries)} entries" )
        logging.debug(f"Some entries:  {allentries[1000:1005]}")
        return allentries


def parse_obo(config):
    """
    creates dict of dicts. key is goterm, contents is dict of 
       
       goterm  ""
       goname  ""
       goasp   ""
       godef   ""
       goisa   [] of goterms
       gohasa  [] of goterms
    
    
    """
    obofile = os.path.expanduser(config.get('ontology','obofile'))
    filehandle = open(obofile)
    godict = {}
    altids = {}
    current = None
    logging.info(f"Parsing file {obofile}")
    try:
        for line in filehandle:
            if line.startswith("[Typedef]"):
                godict[current['goterm']]= current
                break
            elif line.startswith("[Term]"):     
                if current is not None:
                    godict[current['goterm']]= current
                # create new item...
                current = {}
                current['is_a'] = []
                current['part_of'] = []
                
            elif line.startswith("id: "):
                current['goterm'] = line[4:].strip()
                
            elif line.startswith("name: "):
                current['goname'] = line[6:].strip()
            
            elif line.startswith("namespace: "):
                asp = line[11:].strip()
                
                current['goasp'] = GOASPECTMAP[asp]

            # must create a separate mapping from alt_ids that maps to
            # primary, so it can be added to gomatrix properly
            elif line.startswith("alt_id: "):
                #current['alt_id'] = line[8:18].strip()
                #current['goasp'] = GOASPECTMAP[asp]            
                altid = line[8:18].strip()
                altids[altid] = current['goterm'] 
            
            #elif line.startswith("def: "):
            #    current['godef'] = line[5:].strip()

            #elif line.startswith("synonym: "):
            #    current.synonym.append(line[9:].strip())

            elif line.startswith("is_a: "):
                current['is_a'].append(line[6:16].strip())
            
            elif line.startswith("relationship"):
                if "part_of" in line:
                    current['part_of'].append(line[22:32])
                                 
    except Exception as e:
        traceback.print_exc(file=sys.stdout)                
    
    logging.info(f"Parsed file with {len(godict)} terms and {len(altids)} alt terms.")    
    return (godict, altids)


def parse_speclist(config, filepath):
    '''
    Parses uniprot speclist.txt    https://www.uniprot.org/docs/speclist.txt
    to local .CSV
    
    taxonid   species   lineanname       commonname
    72259      ABANI    Abaeis nicippe   Sleepy orange butterfly
                                         
    OXYMO E  475340: N=Oxytenis modestia
                     C=Costa Rica leaf moth
                     S=Dead-leaf moth

    returns 
    
    '''
    logging.debug("Opening species map file %s" % filepath)

    
    try:
        fh = open(filepath, 'r')
    except FileNotFoundError:
        logging.error("No such file %s" % filename)                
   
    species = None
    kingdom = None
    taxonid = None
    lineanname = None
    commonname = None
    
    columnnames = ['species','kingdom','taxonid','lineanname','commonname']
    datalist = []
    # list of tuples
          
    try:
        for line in fh:
            #logging.debug("handling line %s" % line)
            if 'N=' in line and not line.startswith('Code')  :
                #logging.debug("handling N= line. taxonid is %s" % taxonid)
                if species is not None:
                    tup = (species, kingdom, taxonid, lineanname, commonname)
                    #logging.debug("Adding tuple: %s" % str(tup))
                    datalist.append( tup )
                    # reset all varaiables
                    species = kingdom = taxonid = lineanname = commonname = None
                species = line[:5].strip()
                kingdom = line[6]
                taxonid = line[7:15].strip()
                lineanname = line[19:].strip()
                #logging.debug("handling N= line. taxonid is %s" % taxonid)         
            elif 'C=' in line :
                commonname = line[19:].strip()
            elif 'S=' in line :
                 pass
    except Exception as e:
        traceback.print_exc(file=sys.stdout)                
    finally:
        fh.close()    
    logging.debug("Parsed file with %d terms" % len(datalist) )
    return datalist


def get_default_config():
    cp = ConfigParser()
    cp.read(os.path.expanduser("~/git/cafa4/etc/fastcafa.conf"))
    return cp


def run_tocafa(config, infile, outfile=None, modelnum=1):
    """
    Convert prediction .csv into valid CAFA submission file. 
    If no outfile name is provided, put outfile in infile directory, naming 
    according to standard. 
        
    Calculate probability estimate (pest) if needed. 
    Truncate to max_goterms. 
    teamID_modelNo_taxonID_go.txt
    
    <teamname>_<modelnum>_<species>_go.txt
    AUTHOR <teamname>
    MODEL <modelnum>
    KEYWORDS orthlog, expression
    ACCURACY   1  PR=0.75; RC=0.31
    ACCURACY   1  PR=0.65; RC=0.50
    ACCURACY   1  PR=0.55; RC=0.41
    <cid1>    GO:0008550    .99
    .
   <cid2>    GO:0008550    .99
    .
    END

    """
    max_goterms = int(config.get('global','max_goterms'))
    author = config.get('global','author')
    infile = os.path.expanduser(infile)
    df = pd.read_csv(infile, index_col=0)
    logging.debug(f"inbound df=\n{df}")

    s = ""
    s += f"AUTHOR\t{author}\n" 
    s += f"MODEL\t{modelnum}\n"
    s += "KEYWORDS\tortholog, gene expression\n"
    #s += "ACCURACY\t1\tPR=1.00; RC=1.00\n"
    logging.debug(f"dataframe columns={df.columns}" )

    # Dataframe to collect all calculated values. 
    topdf = pd.DataFrame(columns=['cid','goterm','pest'])
    cidlist = list(df.cid.unique())
    logging.debug(f"cid list: {cidlist}")    
    for cid in cidlist:
        cdf = df[df.cid == cid]
        cdf = cdf.nlargest(max_goterms, 'score')
        # Normalize all estimates from score to .01 - .99
        cmax = cdf.score.max()
        cmin = cdf.score.min()
        cdf['pest'] = np.interp(cdf['score'], (cmin, cmax ), (.01,.99))    
        logging.debug(f"cdf=\n{cdf}")
        topdf = topdf.append(cdf, ignore_index=True, sort=False)
    logging.debug(f"topdf=\n{topdf}")
    for (i, row) in topdf.iterrows():
        s += "{}\t{}\t{:.2f}\n".format( row.cid, row.goterm, row.pest)
    s+="END\n"
    try:
        f = open(outfile, 'w')  
        f.write(s)
    except IOError:
        logging.warning(f"error writing to {outfile}")
    finally:
        f.close()
    
    logging.info(f"Wrote cafafile with {len(topdf.index)} entries. ")
    return s

def run_summarize_old(config, infile, outfile):
    """
    Consume and merge a set of eval files, assuming provided parameters with which they 
    were created. 
    Produce f1max distribution graphs appropriately. 

        run_summarize(cp,
                      args.infiles # list
                      args.outfile, # outfile. 
                      args.method,     # phmmer, expression, prior, orthoexpression
                      args.knowledge,  # noknow, limited
                      args.aspects,    # all, bp, mf, cc 
                      args.evcodes,    # exp, iea
)
    """
    logging.debug(f'handling infiles= {infiles}')
    
    topdf = None

    # exp_codes = [ x.strip() for x in config.get('uniprot','exp_codes').split(',')]
    species_order =[x.strip() for x in config.get('summarize', 'species_order').split(',')]
    

    for infile in infiles:
        basename = os.path.basename(infile)
        logging.debug(f"{infile}")
        try:
            df = pd.read_csv(infile, index_col=0, comment="#")            
            df.drop(['f1maxtotal','pr','score'], inplace=True, axis=1)
            df['method'] = method
            df['knowledge'] = knowledge
            df['expcodes'] = evcode
            df['aspect'] = aspect        
            df[['gene','species']] = df.cgid.str.split('_',expand=True)
            #logging.debug("added split gene/species")
            #means = df.drop(['correct'], axis=1)
            #meansdf = means.groupby('cid').mean()
            #sums = df.drop(['nterms','f1max','goterm'], axis=1)
            #sumsdf =  sums.groupby('cid').sum()
            #topdf = pd.merge(df, sumsdf , how='outer', on=['cid'] )
            if topdf is None:
                logging.debug("topdf is none, creating...")
                topdf = pd.DataFrame(columns=list(df.columns))
            logging.debug(f"merging df: {df}") 
            #topdf.append(df, ignore_index=True, sort=True)
            topdf = pd.concat([topdf, df])
            logging.debug(f"topdf after merge: {topdf}")

        except FileNotFoundError:
            logging.error(f"no such file {infile}")

        except ValueError:
            logging.warning(f"Issue reading file {infile}. empty?")

    logging.debug(f"merged infiles: df:\n{topdf}")
    sns.set()
    chart = sns.catplot(x='species', y='f1max', kind='boxen', data=topdf, order=species_order)

    chart.set(title = f"method={method}\nknowledge={knowledge} aspect={aspect} evcode={evcode}  ")
    for ax in chart.axes.flat:
        ax.set(ylim=(0.0,1.0))
        for label in ax.get_xticklabels():
            label.set_rotation(45)
    
    logging.debug(f"chart obj is {chart}")

    logging.debug(f"Saving figure to {outfile}.png ")
    chart.savefig(f"{outfile}.png")
    logging.debug(f"writing merged eval to {outfile}")
    topdf.to_csv(outfile)



def run_summarize(config, infile, outfile):
    '''
    Consume one-row-per-goterm file and produce a one-row-per-target output, calculating f1max.
    Other per-protein statistics? 

    '''
    logging.info(f'handling infile: {infile}')
    f1maxes = []
    try:
        evaldf = pd.read_csv(infile, index_col=0, comment="#")
        logging.debug("Running do_summarize()")
        outdf = do_summarize(config, evaldf )
        logging.debug(f"writing summarize DF to {outfile}")
        outdf.to_csv(outfile)
        logging.debug("done.")
        
    except FileNotFoundError:
        print(f"no such file {infile}")
    

################################ utility functions ###################################

def matrix_info(matrix):
    """
    Returns string of inexpensive info on large matrix for logging. 
    """
    return f"type: {type(matrix)} shape: {matrix.shape} dtype: {matrix.dtype}"

def matrix_debug(matrix):
    """
    Returns string of somewhat more expensive info on large matrix for logging.
    
    """
    return f"type: {type(matrix)} shape: {matrix.shape} dtype: {matrix.dtype} sum: {matrix.sum()} min: {matrix.min()} max {matrix.max()}\n{matrix}"

def filter_goaspect(config, dataframe, goasp):
    """
    Adds aspect column. 
    Removes goterms from incorrect aspect. Re-indexes. 

            cgid             cid      goterm     score  
0      CP270_RAT  G1011600000000  GO:0008150  252622.3
1      CP270_RAT  G1011600000000  GO:0005575  215323.7
2      CP270_RAT  G1011600000000  GO:0110165  210381.7    
       
    """
    df = dataframe
    logging.debug(f"filtering to goaspect {goasp} initial df:\n{df}")
    (godict, altids) = parse_obo(config)   
    df['aspect'] = df.apply(apply_go2aspect, axis=1, godict=godict)
    logging.debug(f"df len={len(df)} with aspect:\n{df}")
    if goasp is not None:
        df = df[df.aspect == goasp]
        df.reset_index(drop=True, inplace=True)
    logging.debug(f"df len={len(df)} after filtering goaspect:\n{df}")
    if len(df) > 0:
        return df
    else:
        logging.debug(f"no data with aspect {goasp}")
        return None


def filter_goevidence(config, dataframe, exponly=True  ):
    """
    Remove any rows of a dataframe that do *not* represent experimentally confirmed 
    evidence.
    XXX  Leave un-annotated (NaN) and experimental goev codes.  
    
    exp_codes = EXP, IDA, IMP, IGI, IEP, IPI, HTP, HDA, HMP, HGI, HEP 
    ele_codes = IEA
    cur_codes = ISS, NAS, TAS, IC, IRD, IGC, IBD, IBA, ISO, ISA, ISM, RCA, IKR
    
    """
    exp_codes = [ x.strip() for x in config.get('uniprot','exp_codes').split(',')]  
    ele_codes = [ x.strip() for x in config.get('uniprot','ele_codes').split(',')]
    cur_codes = [ x.strip() for x in config.get('uniprot','cur_codes').split(',')]
    
    dataframe = dataframe[dataframe['goev'].isin(exp_codes)]
    #exp_goev=[ 'EXP', 'IDA', 'IMP', 'IGI', 'IEP' ]
    #dataframe = dataframe[~dataframe['goev'].isin(nonexp_goev)] 
    dataframe.reset_index(drop=True, inplace=True)
    return dataframe



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

    parser.add_argument('-C', '--usecache',
                        action='store_true', 
                        dest='nocache',
                        default=False, 
                        help='Use cached information.' )
    
    
    subparsers = parser.add_subparsers( dest='subcommand',
                                        help='sub-command help.')

################################ prior prediction ########################################

    parser_prior = subparsers.add_parser('prior',
                                          help='output prior as prediction for input TFAs')
    
    parser_prior.add_argument('-i','--infile', 
                               metavar='infile', 
                               type=str, 
                               help='a .fasta sequence files')

    parser_prior.add_argument('-o','--outfile', 
                               metavar='outfile', 
                               type=str, 
                               help='a DF .csv prediction file')

    parser_prior.add_argument('-s','--species', 
                               metavar='species', 
                               type=str,
                               default=None, 
                               help='species-specific prior, otherwise global')

    parser_prior.add_argument('-V','--version', 
                               metavar='version', 
                               type=str, 
                               default='current',
                               help='uniprot version to use')
    
    parser_prior.add_argument('-g','--goaspect', 
                               metavar='goaspect',
                               choices=['bp','mf','cc','all'], 
                               type=str, 
                               default=None,
                               help='Limit prediction to goaspect [bp|mf|cc|all]' )


################################ phmmer prediction ########################################
   
    parser_phmmer = subparsers.add_parser('phmmer',
                                          help='run phmmer and output prediction')
    
    parser_phmmer.add_argument('-i','--infile', 
                               metavar='infile', 
                               type=str, 
                               help='a .fasta sequence files')

    parser_phmmer.add_argument('-o','--outfile', 
                               metavar='outfile', 
                               type=str, 
                               help='a DF .csv prediction file')

    parser_phmmer.add_argument('-V','--version', 
                               metavar='version', 
                               type=str, 
                               default='current',
                               help='a DF .csv prediction file')

    parser_phmmer.add_argument('-g','--goaspect', 
                               metavar='goaspect', 
                               choices=['bp','mf','cc','all'],
                               type=str, 
                               default=None,
                               help='Limit prediction to goaspect [bp|mf|cc]' )



################################ direct expression ########################################
   
    parser_expr = subparsers.add_parser('expression',
                                          help='infer via expression and output prediction')
    
    parser_expr.add_argument('-i','--infile', 
                               metavar='infile', 
                               type=str, 
                               help='a .fasta sequence files')

    parser_expr.add_argument('-o','--outfile', 
                               metavar='outfile', 
                               type=str, 
                               help='a DF .csv prediction file')

    parser_expr.add_argument('-V','--version', 
                               metavar='version', 
                               type=str, 
                               default='current',
                               help='version of uniprot to use.')

    parser_expr.add_argument('-g','--goaspect', 
                               metavar='goaspect',
                               choices=['bp','mf','cc','all'], 
                               type=str, 
                               default=None,
                               help='Limit prediction to goaspect [bp|mf|cc]' )


################################ orthoexpression ########################################
   
    parser_oexpr = subparsers.add_parser('orthoexpression',
                                          help='infer via expression and output prediction')
    
    parser_oexpr.add_argument('-i','--infile', 
                               metavar='infile', 
                               type=str, 
                               help='a .fasta sequence files')

    parser_oexpr.add_argument('-o','--outfile', 
                               metavar='outfile', 
                               type=str, 
                               help='a DF .csv prediction file')

    parser_oexpr.add_argument('-V','--version', 
                               metavar='version', 
                               type=str, 
                               default='current',
                               help='version of uniprot to use.')

    parser_oexpr.add_argument('-g','--goaspect', 
                               metavar='goaspect', 
                               choices=['bp','mf','cc','all'],
                               type=str, 
                               default=None,
                               help='Limit prediction to goaspect [bp|mf|cc]' )


######################### build and cache ontology ####################################

    parser_buildontology = subparsers.add_parser('build_ontology',
                                          help='build and cache GO ontology')    

    parser_builduniprot = subparsers.add_parser('build_uniprot',
                                          help='build and cache uniprot info')

    parser_builduniprot.add_argument('-V','--version', 
                               metavar='version', 
                               type=str, 
                               default='2019',
                               help='version of uniprot to use.')

######################### build and cache goa gomatrix ####################################

    parser_buildgoa_gomatrix = subparsers.add_parser('build_goa_gomatrix',
                                          help='build and cache GO ontology')    

    parser_buildgoa_gomatrix.add_argument('-V','--version', 
                               metavar='version', 
                               type=str, 
                               default='2019',
                               help='version of GOA to use.')
    
    parser_buildgoa_gomatrix.add_argument('-i', '--infile', 
                               metavar='infile', 
                               type=str, 
                               help='a .gaf GOA file')        

    parser_buildgoa_gomatrix.add_argument('-o', '--outfile',  
                               metavar='outfile',
                               dest='outfile',
                               required=False,
                               default=None, 
                               type=str, 
                               help='a govector matrix file?')


######################### build and cache prior info ####################################
    
    parser_buildprior = subparsers.add_parser('build_prior',
                                          help='build and cache uniprot prior')
    
    parser_buildprior.add_argument('-s', '--species',
                        action="store", 
                        dest='species', 
                        default=None,
                        help='Limit to species where relevant.')

    parser_buildprior.add_argument('-o', '--outfile',
                        action="store", 
                        dest='outfile', 
                        default=None,
                        help='CSV output file of ranked goterms for prior. ')

    parser_buildprior.add_argument('-V','--version', 
                               metavar='version', 
                               type=str, 
                               default='2019',
                               help='version of uniprot to use.')

######################### build and cache species maps ####################################

    parser_buildspecies = subparsers.add_parser('build_species',
                                          help='build and cache NCBI species maps')             

######################### evaluate prediction for known ####################################    
        
    parser_evaluate = subparsers.add_parser('evaluate',
                                          help='evaluate prediction against known. output stats.')
    
    parser_evaluate.add_argument('-i', '--infile', 
                               metavar='infile', 
                               type=str, 
                               help='a .csv prediction file')        
          
    parser_evaluate.add_argument('-o', '--outfile', 
                               metavar='outfile', 
                               type=str, 
                               help='a .csv output file with stats')    

    parser_evaluate.add_argument('-g','--goaspect', 
                               metavar='goaspect', 
                               type=str, 
                               default=None,
                               help='Limit prediction to goaspect [bp|mf|cc|all]' )


    parser_builduniprot_test = subparsers.add_parser('build_uniprot_test',
                                          help='build and cache uniprot test source info')

######################### generate test targets from previous uniprot ####################################

    parser_testset = subparsers.add_parser('testset',
                                          help='generate a set of input .TFA files from species')

    parser_testset.add_argument('-n','--numseq', 
                               metavar='numseq',
                               dest='numseq', 
                               default=10,
                               type=int, 
                               help='Number of sequences to select.')   
    
    parser_testset.add_argument('-s', '--species', 
                               metavar='species',                                
                               dest='species', 
                               required=False,
                               default=None,
                               type=str, 
                               help='a .fasta sequence file with exp. annotated proteins')   
    
    parser_testset.add_argument('-o', '--out',  
                               metavar='outfile',
                               dest='outfile',
                               required=False,
                               default=None, 
                               type=str, 
                               help='a Fasta testfile. ')
    
    parser_testset.add_argument('-L', '--limited',
                        action='store_true', 
                        dest='limited',
                        default=False, 
                        help='Allow test targets that were electronically annotated.')

    parser_testset.add_argument('-O', '--previous',
                        metavar='previous', 
                        dest='previous',
                        default='2017',
                        type=str, 
                        help='Older Swissprot version to use as previous.')

    parser_testset.add_argument('-N', '--current',
                        metavar='current', 
                        dest='current',
                        default='2019', 
                        type=str,
                        help='Newer Swissprot version to use as current.')

################################ combine predictions ########################################
   
    parser_combine = subparsers.add_parser('combine',
                                          help='combine two predictions..')
    
    parser_combine.add_argument('-f','--infile1', 
                               metavar='infile1', 
                               type=str, 
                               help='a .csv prediction')

    parser_combine.add_argument('-s','--infile2', 
                               metavar='infile2', 
                               type=str, 
                               help='a .csv prediction')
    
    parser_combine.add_argument('-o', '--out',  
                               metavar='outfile',
                               dest='outfile',
                               required=False,
                               default=None, 
                               type=str, 
                               help='a DF .csv prediction file')

    parser_combine.add_argument('-w','--weight', 
                               metavar='weight', 
                               type=float,
                               default=1.0, 
                               help='weight of second prediction')

    parser_combine.add_argument('-m', '--method',  
                               metavar='method',
                               dest='method',
                               default='rank_average', 
                               type=str, 
                               help='How predictions should be combined.')
    
######################### write CAFA formatted file for submission ####################################    
        
    parser_tocafa = subparsers.add_parser('tocafa',
                                          help='write CAFA formatted file for submission.')
    
    parser_tocafa.add_argument('-i', '--infile', 
                               metavar='infile', 
                               type=str, 
                               help='a .csv prediction file')        
          
    parser_tocafa.add_argument('-o', '--outfile', 
                               metavar='outfile', 
                               type=str, 
                               help='a .txt output file for submission.')   

    parser_tocafa.add_argument('-m', '--model', 
                               metavar='modelnum', 
                               type=str,
                               default='1', 
                               help='Numeric label for this method/model.')

######################### summarize  ####################################    
    
        
    parser_summarize = subparsers.add_parser('summarize',
                                          help='create target summary from goterm eval. one row per target.  ')
    
    parser_summarize.add_argument('-i', '--infile', 
                               metavar='infile', 
                               type=str, 
                               help='a .csv prediction file')          
          
    parser_summarize.add_argument('-o', '--outfile', 
                               metavar='outfile', 
                               type=str, 
                               default=None,
                               help='a .csv output file for submission.') 



##############################################################################

    args= parser.parse_args()
    
    # default to INFO
    logging.getLogger().setLevel(logging.INFO)
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    cp = ConfigParser()
    cp.read(args.conffile)
    
    logging.debug(f"args: {args}")

    ################## create data resources #######################

    if args.subcommand == 'build_ontology':
        build_ontology(cp, usecache=False)

    if args.subcommand == 'build_uniprot':
        build_uniprot(cp, usecache=False, version=args.version)

    if args.subcommand == 'build_species':
        build_specmaps(cp,  usecache=False)

    if args.subcommand == 'build_prior':
        build_prior(cp, usecache=False, species=args.species, outfile=args.outfile, version=args.version )

    if args.subcommand == 'build_uniprot_test':
        build_uniprot_test(cp, usecache=False )

    if args.subcommand == 'build_goa_gomatrix':
        build_goa_gomatrix(cp, usecache=False, version=args.version, infile=args.infile, outfile=args.outfile)

    ################## generate prediction file #######################
    
    if args.subcommand == 'prior':
        if args.goaspect == 'all':
            args.goaspect = None
        do_prior(cp, args.infile, 
                 args.outfile, 
                 usecache=True, 
                 version=args.version, 
                 goasp=args.goaspect )
       
    if args.subcommand == 'phmmer':
        # make bulk processing easier...
        if args.goaspect == 'all':
            args.goaspect = None
        
        do_phmmer(cp, args.infile, 
                  args.outfile, 
                  usecache=True, 
                  version=args.version,
                  goasp=args.goaspect )
    
    if args.subcommand == 'orthoexpression':
        # make bulk processing easier...
        if args.goaspect == 'all':
            args.goaspect = None
        do_orthoexpression(cp, args.infile, 
                           args.outfile, 
                           usecache=True, 
                           version=args.version,
                           goasp=args.goaspect )

    if args.subcommand == 'expression':
        # make bulk processing easier...
        if args.goaspect == 'all':
            args.goaspect = None
        do_expression(cp, args.infile, 
                      args.outfile, 
                      usecache=True, 
                      version=args.version,
                      goasp=args.goaspect )

    if args.subcommand == 'combine':
        if args.method == 'round_robin':
            run_combine_rr(cp, args.infile1, args.infile2, args.outfile)
        elif args.method == 'rank_average':
            run_combine_weighted(cp, args.infile1, args.infile2, args.outfile, args.weight)

    ################## test performance  #######################

    if args.subcommand == 'testset':
        do_testset(cp, args.numseq, args.species, args.outfile, args.limited, args.previous, args.current )

    if args.subcommand == 'evaluate':
        # make bulk processing easier...
        if args.goaspect == 'all':
            args.goaspect = None
        run_evaluate(cp, args.infile, args.outfile, args.goaspect)

    if args.subcommand == 'summarize':
        run_summarize(cp, args.infile, args.outfile)


    ################## generate submission file  #######################
    
    if args.subcommand == 'tocafa':
        run_tocafa(cp, args.infile, args.outfile, args.model)
        


    