import sys
import csv
import json
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pprint
pp = pprint.PrettyPrinter()
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('max_colwidth', None)

TEST_HITS_FILE_NAME = 'test_hits_pmids.txt'
TEST_OUTPUT_FILE_NAME = 'test_pubmed_info.json'

######################################################################################################
########################     EVALUATE         ########################################################
######################################################################################################

def readDataCol( filename ):
  ids = []
  with open(filename, newline='') as f:
    csvread = csv.reader( f )
    for line in csvread:
      ids.append( line[0] )
  return ids


def readJsonFromFile( filename ):
  data = None
  with open(filename, 'r') as f:
    data = json.load(f)
  return data


def getHits( filename ):
  proposed = set( readDataCol( filename ) )
  actual = set( readDataCol( TEST_HITS_FILE_NAME ) )
  identified = set.intersection( actual, proposed )
  return {
    'actual': actual,
    'proposed': proposed,
    'identified': identified
  }

def getTestArticleCount( filename ):
  PubMedResponse = readJsonFromFile( filename )
  return len( PubMedResponse['PubmedArticle'] )

def getPmfDistr( M, n, N, x ):
  rv = hypergeom( M, n, N )
  xVals = np.arange( 0, n + 1 )
  pmf_hits = rv.pmf( xVals )
  return pmf_hits

def getCdfDistr( M, n, N, x ):
  rv = hypergeom( M, n, N )
  xVals = np.arange( 0, n + 1 )
  cdf_hits = rv.cdf( xVals )
  return cdf_hits

def calCdfDistr( M, n, N, x ):
  prb = hypergeom.cdf(x, M, n, N)
  print(f'x: {x}')
  print(f'F(x): {prb}')
  print(f'1 - F(X): {1 - prb}')

def plotResults(  M, n, N, x ):
  cdf_hits = getCdfDistr( M, n, N, x )
  pmf_hits = getPmfDistr( M, n, N, x )
  xVals = np.arange(0, n + 1)
  fig = plt.figure()
  ax = fig.add_subplot(121)
  ax.plot(xVals, pmf_hits)
  ax.set_xlabel('x (Number hits)')
  ax.set_ylabel('f(x)')
  ax2 = fig.add_subplot(122)
  ax2.plot(xVals, cdf_hits, 'r--')
  ax2.set_xlabel('x (Number hits)')
  ax2.set_ylabel('F(x)')
  plt.show()

def createDist( filename ):
  hits = getHits( filename )
  numTotalArticles = getTestArticleCount( TEST_OUTPUT_FILE_NAME )
  numTotalHits = len( hits['actual'] )
  numProposedArticles = len( hits['proposed'] )
  [ M, n, N ] = [ numTotalArticles, numTotalHits, numProposedArticles ]
  x = len( hits['identified'] )
  print(f'M:{M}; N:{N}; n:{n}; x:{x}')
  calCdfDistr( M, n, N, x )
  plotResults( M, n, N, x )

def evaluate( filename ):
  hits = getHits( filename )
  numActual = len( hits['actual'] )
  numIdentified = len( hits['identified'] )
  print(f'Number of hits available: {numActual}')
  print(f'Number of hits identified: {numIdentified}')
  print(f'Proportion of hits identified: {numIdentified / numActual}')


## Provide the filename with a list of proposed PubMed uids, newline-separated
## useage: python eval.py evaluate.txt
if __name__ == "__main__":
  filename = sys.argv[1]
  if not filename:
    print(f"Useage: python eval.py <path to file>")
  createDist( filename )

