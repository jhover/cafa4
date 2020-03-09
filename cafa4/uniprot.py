#!/usr/bin/env python
#
# Encapsulates bioservices UniProt API and local files for use with CAFA
#

import argparse
from collections import defaultdict
from configparser import ConfigParser
import logging
import os
import pprint
import subprocess
import sys
import traceback

from Bio import SeqIO  # to parse uniprot.dat
from bioservices.uniprot import UniProt  # to query online REST interface
import numpy as np
import pandas as pd
import pdpipe as pdp
import pyarrow as pa
import pyarrow.parquet as pq


#  Cf  https://pypi.org/project/PyUniProt/
class UniProtRecord(object):
    """
    Created from a record object from the bioservices API.
    
    """
    
    def __init__(self, record):
        '''
        
        ID: Q6GZW9
        Name: 006R_FRG3G
        Description: RecName: Full=Uncharacterized protein 006R;
        Database cross-references: EMBL:AY548484, RefSeq:YP_031584.1, SMR:Q6GZW9, GeneID:2947778, KEGG:vg:2947778, Proteomes:UP000008770
        Number of features: 1
        /accessions=['Q6GZW9']
        
        record:
            id:
            name:
            description: 
            annotations: -> dict
            dbxrefs: -> list     
            features: -> SeqFeature

        'letter_annotations', 
        'lower', 
        'name', 'reverse_complement', 'seq', 'translate', 'upper'
        
        E.g 
        annotations -> dict
            'accessions': ['Q6GZW9'], 
            'protein_existence': 4, 
            'date': '28-JUN-2011', 
            'sequence_version': 1, 
            'date_last_sequence_update': '19-JUL-2004', 
            'date_last_annotation_update': '05-JUN-2019', 
            'entry_version': 28, 
            'gene_name': 'ORFNames=FV3-006R;', 
            'organism': 'Frog virus 3 (isolate Goorha) (FV-3)', 
            'taxonomy': ['Viruses', 'Iridoviridae', 'Alphairidovirinae', 'Ranavirus'], 
            'ncbi_taxid': ['654924'], 
            'organism_host': ['Ambystoma (mole salamanders)', 
                               'Dryophytes versicolor (chameleon treefrog)', 
                               'Lithobates pipiens (Northern leopard frog) (Rana pipiens)', 
                               'Notophthalmus viridescens (Eastern newt) (Triturus viridescens)',
                               'Rana sylvatica (Wood frog)'], 
                               'host_ncbi_taxid': ['8295', 
                                                    '30343', 
                                                    '8404', 
                                                    '8316', 
                                                    '45438'], 
                                'references': --> Reference 
        
        dbxrefs -> list
            ['EMBL:AY548484', 'RefSeq:YP_031579.1', 'SwissPalm:Q6GZX4', 'GeneID:2947773', 
            'KEGG:vg:2947773', 'Proteomes:UP000008770', 
            'GO:GO:0046782', 
            'InterPro:IPR007031', 'Pfam:PF04947']
    
        [Reference(title='Comparative genomic analyses of frog virus 3, type species of the genus Ranavirus (family Iridoviridae).',
                     ...)], 
                  'keywords': ['Complete proteome', 
                              'Reference proteome']
        '''
        self.proteinid = record.id
        self.protein = record.name
        self.goterms = []
        for xf in record.dbxrefs:
            if xf.startswith("GO:"):
                gt = xf[3:]
                self.goterms.append(gt)
        self.accessions = record.annotations['accessions']
        self.taxonid = record.annotations['ncbi_taxid'][0]
        #self.species = 
        #if 'accessions' in record.annotations.keys():
        #    self.accessions = record.annotations['accessions']
            
            
    def __repr__(self):
        REPR_ATTRS=['proteinid','protein','accessions', 'taxonid','goterms']
        s="UniProtRecord: "
        for atr in REPR_ATTRS:
            s+="%s=%s " % (atr, 
                           self.__getattribute__(atr))
        return s


class UniProt(object):
    '''
    Aux info plugin.
    Takes dataframe, extracts entry_ids, adds info from uniprot.  
    Returns modified dataframe. 

    '''
    
    ASPECTMAP = { 'C': 'cc',
                  'F': 'mf',
                  'P': 'bp'
                }
    

    def __init__(self, config):
        self.log = logging.getLogger(self.__class__.__name__)
        self.config = config
        self.uniprotapi = None
        self.outdir = os.path.expanduser( config.get('global','outdir') )
        self.taxid_mapfile = os.path.expanduser( config.get('global','taxid_mapfile'))
        self.sprotdatfile = os.path.expanduser( config.get('ontologyplugin','sprotdatfile') )
        self.cachedir = os.path.expanduser( config.get('ontologyplugin','cachedir') )
        excodes = config.get('ontologyplugin' , 'excluded_evidence_codes', fallback=[] ).split(',')
        excodes = [ x.strip() for x in excodes]
        self.excluded_evidence_codes = excodes
        self.sprotdf = None
        self.udf = None
        self.tdf = pd.read_csv(self.taxid_mapfile, index_col=0)
        
        # Create easy lookup mappings from taxon data frame...
        itdf = self.tdf.set_index('taxonid')
        self.taxiddict = itdf.to_dict(orient='index')
        
        isdf = self.tdf.set_index('species')
        self.specdict = isdf.to_dict(orient='index')         
        self.log.debug("UniProtGOlugin initialized.")
        

    def cafa_execute(self, dataframe, online=False):
        """
        Takes inbound dataframe of orthologs and adds in GO terms and evidence codes from 
        uniprot/swissprot.
        For a given ortholog protein, one row is added for each GO term.
        Returns new dataframe with all info.   
        
        """
        #
        # inbound:
        #            cafaid         evalue  score  bias  db proteinacc protein species cafaprot cafaspec
        # 0   T100900000001  1.100000e-156  523.6   8.5  sp    Q9CQV8   1433B   MOUSE    1433B    MOUSE
        # 1   T100900000001  4.100000e-155  518.4   7.7  sp    P35213   1433B     RAT    1433B    MOUSE    
        # 2   T100900000001  5.400000e-155  518.0   7.2  sp    A4K2U9   1433B   PONAB    1433B    MOUSE
        # 3   T100900000001  5.400000e-155  518.0   7.2  sp    P31946   1433B   HUMAN    1433B    MOUSE
        
        # Get all unique target accession numbers.        
        entries = dataframe['proteinacc'].unique().tolist()
        # Look up GOterms in uniprot...
        if online:
            self.uniprotapi = UniProt()
            self.log.debug("Querying uniprot API for %d unique entries" % len(entries))
            self.udf = self.uniprotapi.get_df(entries)
            self.log.debug(f"\n{self.udf}")
            self.udf.to_csv("%s/uniprot.csv" % self.outdir)
            udfslim = self.udf[['Entry','Gene ontology IDs']] 
            # df.tacc corresponds to udf.Entry  ...
            #  entry == proteinid
            #  gene ontology id = goterm
            #
            self.log.debug("Making new rows for each goterm.")
            newrowdict = {} 
            ix = 0 
            for row in udfslim.itertuples():
                (entry, golist) = row[1:] 
                for goterm in golist: 
                    #print("creating new row: %s : %s %s %s" % (ix, entry, gene, goterm)) 
                    newrow = [ entry, goterm ] 
                    newrowdict[ix] = newrow 
                    ix += 1       
            
            godf = pd.DataFrame.from_dict(newrowdict, orient='index', columns=['entry','goterm']) 

        else:
            self.log.debug("Using offline functionality...")
            godf = self.get_swissprot_df(usecache = True)
            self.log.debug(f"GO DataFrame:\n{godf}")
            #    proteinid   proteinacc    goterm      goaspect goevidence
            # 0  001R_FRG3G  Q6GZX4      GO:0046782    bp        IEA
            # 1  002L_FRG3G  Q6GZX3      GO:0033644    cc        IEA
            
        # For each go term add row...
        newdfdict= {}
        ix = 0
        for row in dataframe.itertuples():
            self.log.debug("inbound row = %s" % str(row))
            #(query, evalue, score, bias, db, tacc, protein, species) = row[1:]
            (cafaid, evalue, score, bias, db, proteinacc, protein, species, cafaprot, cafaspec) = row[1:]
            self.log.debug(f"Searching for match for '{proteinacc}'")
            gomatch = godf[ godf.proteinacc == proteinacc ]
            self.log.debug(f"gomatch is:\n {gomatch}")
            for gr in gomatch.itertuples():
                (entry, proteinacc, protein, species, goterm, goaspect, goevidence) = gr[1:]
                newrow = [cafaid, evalue, score, bias, db, 
                          proteinacc, protein, species, cafaprot, 
                          cafaspec, goterm, goaspect, goevidence ]
                newdfdict[ix] = newrow
                ix += 1
        
        newdf = pd.DataFrame.from_dict(newdfdict, orient='index', columns = ['cafaid', 
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
                                                                             'goevidence'                                                                             
                                                                             ])
        for xc in self.excluded_evidence_codes:
            self.log.debug(f"{len(newdf.index)} rows. Removing evidence code {xc}...")
            #newdf = newdf[newdf.goevidence != xc]
            newdf.drop(newdf.loc[newdf['goevidence'] == xc ].index, inplace=True)
            self.log.debug(f"{len(newdf.index)} rows after.")
            self.log.debug(f"\n{str(newdf)}")  
     
        return newdf
        # Output:
        #             cafaid         evalue  score  bias  db proteinacc protein species cafaprot cafaspec      goterm goaspect goevidence
        # 0    T100900000001  1.100000e-156  523.6   8.5  sp     Q9CQV8   1433B   MOUSE    1433B    MOUSE  GO:0005737       cc        ISO
        # 1    T100900000001  1.100000e-156  523.6   8.5  sp     Q9CQV8   1433B   MOUSE    1433B    MOUSE  GO:0005829       cc        ISO
        # 2    T100900000001  1.100000e-156  523.6   8.5  sp     Q9CQV8   1433B   MOUSE    1433B    MOUSE  GO:0042470       cc        IEA
        #
     

    def _dat2upr(self):
        self.log.debug("opening swissprot dat file %s" % self.sprotdatfile)
        rgen = SeqIO.parse(self.sprotdatfile,"swiss")
        i = 0
        uprlist = []
        self.log.debug("Completed SeqIO.parse(). Handling records...")
        for record in rgen:
            upr = UniProtRecord(record)
            uprlist.append(upr)
            #print(record)
            i += 1
            if i % 10000 == 0:
                self.log.debug("Handled %d records..." % i)
            #    break
        self.log.debug("parsed dat file of %d records" % len(uprlist))
        return uprlist


    def get_annotation_df(self):
        self.log.debug("opening swissprot dat file %s" % self.sprotdatfile)
        rgen = SeqIO.parse(self.sprotdatfile,"swiss")
        self.log.debug("rgen type is %s" % type(rgen))
        #self.log.debug("Created generator with %d records" % len(rgen))
        i = 0
        alltuples = []
        for record in rgen:
            #print(record)
            i += 1
            if i % 1000 == 0:
                self.log.debug("Handled %d records..." % i)
            goterms = []
            for xf in record.dbxrefs:
                if xf.startswith("GO:"):
                    gt = xf[3:]
                    goterms.append(gt)
            if len(goterms) > 0:
                proteinid = record.id
                protein = record.name
                taxonid = record.annotations['ncbi_taxid'][0]                
                for gt in goterms:
                    t = (taxonid, proteinid, protein, gt)
                    alltuples.append(t)
                # fan out over goterms
            else:
                # ignore un-annotated entries. 
                pass
                        
            if i >= 1000:
                break
        #self.log.debug("generated %d tuples" % len(alltuples)) 
        self.log.debug( f"Generated { len(alltuples) } tuples") 
        df = pd.DataFrame(alltuples, columns=['taxonid','proteinid','protein','goterm'])
        
        return df


##########################################
#
#   Non-cafalib usage (NOT using API)
#
##########################################

    def get_swissprot_df(self, usecache = True):
        """
        Get swissprot info as dataframe from files, without API, one row per GOterm.
       
        Fields:
           proteinid protein taxonid goterm goaspect goevidence 
      
        self.proteinid = record.id
        self.proteinacc = record. ?
        self.protein = record.name
        self.goterms = []
        for xf in record.dbxrefs:
            if xf.startswith("GO:"):
                gt = xf[3:]
                self.goterms.append(gt)
        self.accessions = record.annotations['accessions']
        self.taxonid = record.annotations['ncbi_taxid'][0]
        
        """

        cachepath = f"{self.cachedir}/sprotgolist.csv"
        if usecache:
            if os.path.exists(cachepath):
                self.sprotdf = pd.read_csv(cachepath, index_col=0)
                self.log.debug(f"Loaded dataframe from cache: {cachepath}")
        if self.sprotdf is not None:
            self.log.debug("Cache hit. Using DataFrame from cache...")
        else:
            self.log.debug("Getting dictionary list...")
            dlist = self._handle_swissprot_file()
            self.log.debug(f"Got dict list of {len(dlist)} entries. Creating dataframe...")
            self.sprotdf = pd.DataFrame(dlist)
            #self.sprotdf.set_index('proteinacc', inplace = True)
            self.log.debug(f"Made dataframe:\n {str(self.sprotdf)}")
            self.log.info(f"Saving dataframe to cache file: {cachepath}")        
            self.sprotdf.to_csv(cachepath)
        return self.sprotdf
    
    
    def _handle_swissprot_file(self):
        '''
         Read uniprot_sprot.dat and return list of dicts of relevant fields.
    

        '''
        self.log.debug("Handling swissprot file...")
        filehandle = None
        try:
            self.log.info(f"Opening file {self.sprotdatfile}" )
            filehandle = open(self.sprotdatfile, 'r')
            self.log.debug("File opened. Parsing...")
            dlist = self._parsefile(filehandle)
            filehandle.close()
                
        except FileNotFoundError:
            self.log.error("No such file %s" % filename)        
        
        finally:
            if filehandle is not None:
                filehandle.close()
        self.log.debug("Parsed data file.")
        return dlist


    def _parsefile(self, filehandle):
        """
        Parses sprot DAT file and fans out goterms to list of dicts. 
    
        """
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
            #for line in filehandle:
                if line.startswith("ID "):
                    # ID   001R_FRG3G              Reviewed;         256 AA.
                    #      <prot_name>_<prot_spec>
                    proteinid = line[5:16].strip()
                    current = defaultdict(dict)
                    current['proteinid'] = proteinid
                    (protein, species) = proteinid.split('_')
                    current['protein'] = protein
                    current['species'] = species
                    self.log.debug("Handling ID. New entry.")                
                
                elif line.startswith("AC "):
                    # AC   Q6GZX4;
                    # AC   Q91896; O57469;
                    self.log.debug("Handling AC.")
                    accession = line[5:11].strip()
                    current['proteinacc'] = accession

                elif line.startswith("OX   "):
                    #OX   NCBI_TaxID=654924;
                    self.log.debug("Handling OX.")
                    taxonid = ""
                    val = line[5:]
                    fields = val.split('=')
                    if fields[0] == 'NCBI_TaxID':
                        taxonid = fields[1].strip().replace(';','')
                    current['taxonid'] = taxonid
                    
                elif line.startswith("DR   GO;"):
                    # DR   GO; GO:0046782; P:regulation of viral transcription; IEA:InterPro.
                    # P biological process, C cellular component, F molecular function.  
                    self.log.debug("Handling DR.")
                    fields = line.split(';')
                    goterm = fields[1].strip()
                    goinfo = fields[2]
                    aspcode = goinfo.split(':')[0].strip()
                    goaspect = UniProt.ASPECTMAP[aspcode]
                    goevsrc = fields[3]
                    (goevidence, evsrc) = goevsrc.split(':') 
                    goevidence = goevidence.strip()
                    current['goterms'][goterm] = [ goaspect, goevidence]

                elif line.startswith("SQ   SEQUENCE"):
                    self.log.debug("Handling SQ:  XXX")
                    # line = filehandle.readline()

                elif line.startswith("GN   "):
                    # Examples:
                    #  GN   ABL1 {ECO:0000303|PubMed:21546455},
                    #  GN   Name=BRCA1; Synonyms=RNF53;
                    #  GN   ORFNames=T13E15.24/T13E15.23, T14P1.25/T14P1.24;
                    #   
                    
                    self.log.debug("Handling GN.")
                    val = line[5:]

            
                elif line.startswith("//"):
                    self.log.debug("End of entry.")
                    clist = self._handle_current(current)
                    current = None
                    allentries.extend(clist)
                    self.log.debug(f"All entries list now {len(allentries)} items... ")
                    if len(allentries) >= repthresh:
                        self.log.info(f"Processed {len(allentries)} entries... ")
                        sumreport +=1
                        repthresh = sumreport * suminterval
                
        except Exception as e:
            traceback.print_exc(file=sys.stdout)                
        
        self.log.info(f"Parsed file with {len(allentries)} goterms" )
        return allentries

    
    def _handle_current(self, currentinfo):
        """
        takes dictionary:
        currentinfo = { 'proteinid' : 'x', 'protein' : 'xxx' , 'goterms' :  { 'GO:0005634' : [ 'C' , 'HDA' ],
                                                                              'GO:0005886' : [ 'C' ,'HDA'],
                                                                              }                                                                                              
                        } 
        
        returns list of dicts:
                     [  { 'proteinid' : 'x', 'protein' : 'xxx' , 'goterm' : 'GO:0005634',
                                                                           'goaspect':'cc',
                                                                           'goevidence': 'HDA' },
                       { 'proteinid' : 'x', 'protein' : 'xxx' , 'goterm' : 'GO:0005886',
                                                                           'goaspect':'cc',
                                                                           'goevidence': 'HDA' },                                                                           
                      ]
        """
        self.log.debug(f'handling {currentinfo} ')
        newlist = []
        gtdict = currentinfo['goterms']
        for gt in gtdict.keys():
            self.log.debug(f"Handling term {gt}")
            newdict = {}
            newdict['proteinid'] = currentinfo['proteinid']
            newdict['proteinacc'] = currentinfo['proteinacc']
            newdict['protein'] = currentinfo['protein']
            newdict['species'] = currentinfo['species']
            newdict['goterm'] = gt
            newdict['goaspect'] = currentinfo['goterms'][gt][0]             
            newdict['goevidence'] = currentinfo['goterms'][gt][1]        
            newlist.append(newdict)
        
        self.log.debug(f"Created fanout of length: {len(newlist)}")
        return newlist
              

    def _make_species_map(self):
        '''
        Parses uniprot speclist.txt    https://www.uniprot.org/docs/speclist.txt
        to local .CSV
        
        taxonid   species   lineanname       commonname
        72259      ABANI    Abaeis nicippe   Sleepy orange butterfly
                                             
        OXYMO E  475340: N=Oxytenis modestia
                         C=Costa Rica leaf moth
                         S=Dead-leaf moth
        
        '''
        listfile = self.speciesmap
        self.log.debug("Opening species map file %s" % listfile)
        try:
            fh = open(listfile, 'r')
        except FileNotFoundError:
            self.log.error("No such file %s" % filename)                
       
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
                #self.log.debug("handling line %s" % line)
                if 'N=' in line and not line.startswith('Code')  :
                    #self.log.debug("handling N= line. taxonid is %s" % taxonid)
                    if species is not None:
                        tup = (species, kingdom, taxonid, lineanname, commonname)
                        #self.log.debug("Adding tuple: %s" % str(tup))
                        datalist.append( tup )
                        # reset all varaiables
                        species = kingdom = taxonid = lineanname = commonname = None
                    species = line[:5]
                    kingdom = line[6]
                    taxonid = line[7:15].strip()
                    lineanname = line[19:].strip()
                    #self.log.debug("handling N= line. taxonid is %s" % taxonid)         
                elif 'C=' in line :
                    commonname = line[19:].strip()
                elif 'S=' in line :
                     pass
        except Exception as e:
            traceback.print_exc(file=sys.stdout)                
        finally:
            fh.close()
        
        self.log.debug("Parsed file with %d terms" % len(datalist) )
    
        df = pd.DataFrame(datalist, columns = columnnames)
        outfile = "%s/speclist.csv" % self.outdir
        self.log.debug("Writing dataframe to %s" % outfile )
        df.to_csv(outfile )
        print(str(df))
        return df

    @classmethod
    def get_default_df(cls, usecache=True):
        cp = ConfigParser()
        cp.read(os.path.expanduser('~/git/cafa4/etc/cafa4.conf'))
        upg = UniProt(cp)
        df = upg.get_swissprot_df(usecache = usecache)

        return df 

    @classmethod
    def calculate_prior(cls, dataframe, species=None, goaspect=None):
        """
        @arg 
           dataframe :  standard internal dataframe, 
           species  :  NCBI species code   e.g. MOUSE | HUMAN
           goaspect : internal aspect code   e.g. [cc | bp | mf ]
           
           proteinid proteinacc protein species      goterm goaspect goevidence
           11K_PAVHV     P0DJZ0     11K   PAVHV  GO:0030430       cc        IDA
           ...

        returns:
            dataframe w/ ranked list of goterms, within the specified species/aspect if supplied.
            otherwise globally 
            
            goterm      goaspect    count    prob
            GO:0045735  cc           3679    .142
            GO:0030433  bp           1256    .086

        """
        df = dataframe
        if species is not None:
            df = df[df.species == species] 
        if goaspect is not None:
            df = df[df.goaspect == goaspect]

        totalterms = df.goterm.count()
        newdf = pd.DataFrame( df.goterm.value_counts() ).reset_index()
        newdf.columns = ['goterm','counts']
        
        

def test_uniprot(config):
    upg = UniProt(config)
    entrylist = ['Q9CQV8', 'P35213', 'A4K2U9', 'P31946', 'Q4R572', 'P68250']
    out = upg._query_entries(entrylist)
    print(out)     

def test_datparse(config):
    upg = UniProt(config)
    df = upg.get_swissprot_df()
    return df

def test_speciesmap(config):
    upg = UniProt(config)
    upg._make_species_map()

def test_testset(config):
    upg = UniProt(config)
    df = upg.get_annotation_df()
    return df

def test_swissprot(config):
    logging.debug("Running test_swissprot")
    upg = UniProt(config)
    out = upg.get_swissprot_df()
    print(str(out))

    


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
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    cp = ConfigParser()
    cp.read(args.conffile)

    #test_uniprot(cp)
    #df = test_testset(cp)
    #df.to_csv('testset.csv')
    #print(str(df))
    test_speciesmap(cp)
    #test_swissprot(cp)

