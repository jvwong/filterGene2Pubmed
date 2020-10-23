
from Bio import Entrez
import csv
import random
import math
import json
import numpy as np
import pandas as pd
import pprint

pp = pprint.PrettyPrinter()
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('max_colwidth', None)

ALL_INPUT_FILE_NAME = 'all_pmids.txt'
HITS_INPUT_FILE_NAME = 'hits_pmids.txt'

TEST_OUTPUT_FILE_NAME = 'test_set_2.txt'
TRAIN_OUTPUT_FILE_NAME = 'train_set_2.txt'

TEST_HITS_FILE_NAME = 'test_hits_pmids.txt'
TRAIN_HITS_FILE_NAME = 'train_hits_pmids.txt'

Entrez.email = "support@biofactoid.org"

######################################################################################################
########################     RETRIEVE JOURNAL PMIDS  ##################################################
######################################################################################################

# journal = '1083-351X'
# volumes = range(295, 296)
# issues = range(31, 38, 1)

# def esearch( journal, volume, issue=None ):
#   term = '{journal}[ta] {volume}[vi] Journal Article[pt]'.format( journal=journal, volume=volume )
#   if issue:
#     term += ' {issue}[ip]'.format(issue=issue)
#   print (term)
#   handle = Entrez.esearch( db='pubmed', term=term, retmax=10000 )
#   record = Entrez.read( handle )
#   handle.close()
#   return record

# def getPmidsForJournal():
#   ids = []
#   for volume in volumes:
#     for issue in issues:
#       record = esearch( journal, volume, issue )
#       ids += record['IdList']

#   with open('{issn}.txt'.format(issn=journal), 'a') as f:
#     for item in ids:
#       f.write('{item}\n'.format(item=item))


######################################################################################################
########################     CREATE          ########################################################
######################################################################################################


def efetch( ids ):
  id = ','.join( ids )
  handle = Entrez.efetch( db='pubmed', id=id, retmode='xml' )
  record = Entrez.read( handle )
  handle.close()
  return record

def sample( seq, count ):
  return set( random.sample( list( seq ), count ) )

def split( ids_A, ids_H ):
  """
  Notation:
    - A: set of all article ids
    - H: subset of article ids designated as a 'hit' (i.e. Good for Biofactoid)
    - M: subset of article ids designated as not a 'hit' (H complement or 'miss')
    - T: subset of articles ids for training
    - E: subset of articles used for evaluation (aka testing)
  """
  A = set( ids_A )
  H = set( ids_H )
  num_hits = len( H )
  T_size = math.floor( num_hits / 2 )
  M = A - H
  T_hits = sample( H, T_size )
  T_miss = sample( M, T_size )
  E_hits = H - T_hits
  E_miss = M - T_miss
  return {
    'T': list( T_hits | T_miss ),
    'E': list( E_hits | E_miss )
  }

def efetch2Tsv( efetch_response, ids_H, filename ):
  articles = efetch_response['PubmedArticle']
  H = set( ids_H )
  pmids = []
  texts = []
  hits = []
  for PubmedArticle in articles:
    try:
      MedlineCitation = PubmedArticle['MedlineCitation']
      PMID = MedlineCitation['PMID']
      Article = MedlineCitation['Article']
      ArticleTitle = Article['ArticleTitle']
      Abstract = ''
      if 'Abstract' in Article:
        Abstract = ' '.join( Article['Abstract'][ 'AbstractText' ] )
    except KeyError as err:
      print(f'KeyError {err}')
      print(f'PubmedArticle {json.dumps(PubmedArticle, indent=2)}')
    else:
      text = ArticleTitle + Abstract
      isHit = 1 if PMID in H else 0
      pmids.append( PMID )
      texts.append( text )
      hits.append( isHit )
  df = pd.DataFrame({
    'pmid': pmids,
    'text': texts,
    'hit': hits
  })
  df.to_csv( filename, sep='\t', index=False, quoting=csv.QUOTE_NONNUMERIC )


def createTestTrainPubMedData():
  """
  Notation:
    - T: subset of articles ids for training
    - E: subset of articles used for evaluation (aka testing)
  """
  ids_A = pd.read_csv( ALL_INPUT_FILE_NAME, dtype={'pmids': str} )['pmids']
  ids_H = pd.read_csv( HITS_INPUT_FILE_NAME, dtype={'pmids': str} )['pmids']
  id_lists = split( ids_A, ids_H )
  T_efetch_response = efetch( id_lists['T'] )
  efetch2Tsv( T_efetch_response, ids_H, TRAIN_OUTPUT_FILE_NAME )
  E_efetch_response = efetch( id_lists['E'] )
  efetch2Tsv( E_efetch_response, ids_H, TEST_OUTPUT_FILE_NAME )
  return id_lists

# def readJsonFromFile( filename ):
#   data = None
#   with open(filename, 'r') as f:
#     data = json.load(f)
#   return data

# def getHits( data ):
#   hits = readDataCol( HITS_INPUT_FILE_NAME )
#   pmids =  [ PubmedArticle['MedlineCitation']['PMID'] for PubmedArticle in data['PubmedArticle'] ]
#   common = set.intersection( set( hits ), set( pmids ) )
#   return common

# def write2ColFile( filename, data ):
#   with open(filename, 'w') as f:
#     for item in data:
#       f.write("%s\n" % item)

# def getTestTrainHits():
#   train = readJsonFromFile( TRAIN_OUTPUT_FILE_NAME )
#   train_hits = getHits( train )
#   write2ColFile( TRAIN_HITS_FILE_NAME, list( train_hits ) )
#   test = readJsonFromFile( TEST_OUTPUT_FILE_NAME )
#   test_hits = getHits( test )
#   write2ColFile( TEST_HITS_FILE_NAME, list( test_hits ) )


if __name__ == "__main__":
  id_lists = createTestTrainPubMedData()
  # getTestTrainHits()