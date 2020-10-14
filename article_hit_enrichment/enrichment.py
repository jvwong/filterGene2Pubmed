import pandas as pd
from Bio import Entrez
import csv
import random
import math
import pprint
import json
pp = pprint.PrettyPrinter()
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('max_colwidth', None)


ALL_INPUT_FILE_NAME = './article_hit_enrichment/all_pmids.txt'
HITS_INPUT_FILE_NAME = './article_hit_enrichment/hits_pmids.txt'

TEST_OUTPUT_FILE_NAME = './article_hit_enrichment/test_pubmed_info.json'
TRAIN_OUTPUT_FILE_NAME = './article_hit_enrichment/train_pubmed_info.json'
TEST_HITS_FILE_NAME = './article_hit_enrichment/test_hits_pmids.txt'
TRAIN_HITS_FILE_NAME = './article_hit_enrichment/train_hits_pmids.txt'

Entrez.email = "support@biofactoid.org"

######################################################################################################
########################     CREATE          ########################################################
######################################################################################################


def efetch( ids ):
  id = ','.join( ids )
  handle = Entrez.efetch( db='pubmed', id=id, retmode='xml' )
  record = Entrez.read( handle )
  handle.close()
  return record

def readDataCol( filename ):
  ids = []
  with open(filename, newline='') as f:
    csvread = csv.reader( f )
    for line in csvread:
      ids.append( line[0] )
  return ids

def split( data, k=2 ):
  size = math.floor( len(data) / k )
  test = random.sample(data, size)
  train = set( data ) - set( test )
  return {
    'test': test,
    'train': train
  }

def writeJson2File( filename, data ):
  with open(filename, 'w') as f:
    json.dump(data, f, indent=2, sort_keys=True)


# def main():
#   ids = readDataCol( ALL_INPUT_FILE_NAME )
#   sets = split( ids, k=2 )
#   test_efetch_response = efetch( sets['test'] )
#   writeJson2File( TEST_OUTPUT_FILE_NAME, test_efetch_response )
#   train_efetch_response = efetch( sets['train'] )
#   writeJson2File( TRAIN_OUTPUT_FILE_NAME, train_efetch_response )

# main()

######################################################################################################
########################     EVALUATE         ########################################################
######################################################################################################

def readJsonFromFile( filename ):
  data = None
  with open(filename, 'r') as f:
    data = json.load(f)
  return data

def getHits( data ):
  hits = readDataCol( HITS_INPUT_FILE_NAME )
  pmids =  [ PubmedArticle['MedlineCitation']['PMID'] for PubmedArticle in data['PubmedArticle'] ]
  common = set.intersection( set( hits ), set( pmids ) )
  return common

def write2ColFile( filename, data ):
  with open(filename, 'w') as f:
    for item in data:
      f.write("%s\n" % item)

def main():
  train = readJsonFromFile( TRAIN_OUTPUT_FILE_NAME )
  train_hits = getHits( train )
  write2ColFile( TRAIN_HITS_FILE_NAME, list( train_hits ) )
  test = readJsonFromFile( TEST_OUTPUT_FILE_NAME )
  test_hits = getHits( test )
  write2ColFile( TEST_HITS_FILE_NAME, list( test_hits ) )


main()