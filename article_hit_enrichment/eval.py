import sys
import csv

TEST_HITS_FILE_NAME = './article_hit_enrichment/test_hits_pmids.txt'
TRAIN_HITS_FILE_NAME = './article_hit_enrichment/train_hits_pmids.txt'


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


def getHits( filename ):
  proposed = readDataCol( filename )
  actual = readDataCol( TEST_HITS_FILE_NAME )
  identified = set.intersection( set( actual ), set( proposed ) )
  return {
    'actual': actual,
    'identified': identified
  }

def evaluate( filename ):
  hits = getHits( filename )
  numActual = len( hits['actual'] )
  numIdentified = len( hits['identified'] )
  print(f'Number of hits available: {numActual}')
  print(f'Number of hits identified: {numIdentified}')
  print(f'Proportion of hits identified: {numIdentified / numActual}')


## Provide the filename with a list of proposed PubMed uids, newline-separated
## useage: python article_hit_enrichment/eval.py article_hit_enrichment/evaluate.txt
if __name__ == "__main__":
  filename = sys.argv[1]
  if not filename:
    print(f"Useage: python eval.py <path to file>")
  evaluate( filename )

