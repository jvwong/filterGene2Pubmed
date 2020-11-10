import requests
import csv
import json
import numpy as np
import pandas as pd
import pprint

pp = pprint.PrettyPrinter()
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('max_colwidth', None)

QUERY = '9887103'

pos_examples = [
'''9813169''',
'''9887103''',
'''10320478''',
'''10357889''',
'''15473904''',
'''10831835''',
'''12581521''',
'''16098524''',
'''16437159''',
'''16951053''',
'''17119021''',
'''17119022''',
'''17434471''',
'''17651718''',
'''18054857''',
'''18077144''',
'''18171686''',
'''18588885''',
'''18820452''',
'''19004854''',
'''19766717''',
'''20346940''',
'''21613324''',
'''22563507''',
'''22874913''',
'''23293296''',
'''23345225''',
'''23416177''',
'''23485686''',
'''23664613''',
'''24244197''',
'''24462039''',
'''24680699''',
'''24706779''',
'''25255438''',
'''25284780''',
'''25633355''',
'''28108836''',
'''28123053''',
'''28178568''',
'''28739899''',
'''28748922''',
'''28853393''',
'''30049242'''
]
neg_examples = [
'''10712925''',
'''10731132''',
'''10733523''',
'''12537568''',
'''12537572''',
'''16110336''',
'''12537573''',
'''12537574''',
'''14605208''',
'''15018941''',
'''15046717''',
'''17569856''',
'''17569867''',
'''20220848''',
'''20800474''',
'''22936248''',
'''23071443''',
'''23072462''',
'''25312911''',
'''26109356''',
'''26109357''',
'''24793180''',
'''24912777''',
'''28130362''',
'''26667894''',
'''28421183''',
'''29615416'''
]

def doRanking( query, documents ):
  SEMANTIC_SEAERCH_URL = 'http://localhost:8000/'
  data = { 'query': query, 'documents': documents }
  r = requests.post( SEMANTIC_SEAERCH_URL, json = data )
  return r.json()


def createTestData():
  df_pos = pd.DataFrame({ 'pmid': pos_examples, 'indicator': [1] * len( pos_examples )  })
  df_neg = pd.DataFrame({ 'pmid': neg_examples, 'indicator': [0] * len( neg_examples )  })
  df_examples = pd.concat( [ df_pos, df_neg ], axis = 0, ignore_index = True )
  # df_examples.to_csv( 'examples.csv', sep='\t', index=False, quoting=csv.QUOTE_NONNUMERIC )
  return df_examples


def evaluateTestData( df_test ):
  documents = list( df_test['pmid'] )
  ranks = doRanking( QUERY, documents )
  df_results = pd.DataFrame.from_records( ranks )
  df_results = df_results.rename(columns={'uid': 'pmid'})
  return df_results


if __name__ == "__main__":
  OUTPUT_FILE = 'declutr-sci-base'
  df_test = createTestData()
  df_results = evaluateTestData( df_test )
  df_merged = df_results.merge( df_test, how = 'outer' )
  df_merged.to_csv( OUTPUT_FILE, sep='\t', index=False, quoting=csv.QUOTE_NONNUMERIC )
