
from Bio import Entrez
import csv
import random
import math
import json

ALL_INPUT_FILE_NAME = 'all_pmids.txt'
HITS_INPUT_FILE_NAME = 'hits_pmids.txt'

TEST_OUTPUT_FILE_NAME = 'test_pubmed_info.json'
TRAIN_OUTPUT_FILE_NAME = 'train_pubmed_info.json'
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


def createTestTrainPubMedData():
  ids = readDataCol( ALL_INPUT_FILE_NAME )
  sets = split( ids, k=2 )
  test_efetch_response = efetch( sets['test'] )
  writeJson2File( TEST_OUTPUT_FILE_NAME, test_efetch_response )
  train_efetch_response = efetch( sets['train'] )
  writeJson2File( TRAIN_OUTPUT_FILE_NAME, train_efetch_response )


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

def getTestTrainHits():
  train = readJsonFromFile( TRAIN_OUTPUT_FILE_NAME )
  train_hits = getHits( train )
  write2ColFile( TRAIN_HITS_FILE_NAME, list( train_hits ) )
  test = readJsonFromFile( TEST_OUTPUT_FILE_NAME )
  test_hits = getHits( test )
  write2ColFile( TEST_HITS_FILE_NAME, list( test_hits ) )


if __name__ == "__main__":
  pass
  # createTestTrainPubMedData()
  # getTestTrainHits()