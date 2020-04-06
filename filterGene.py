import pandas as pd

GI1 = 831748
GI2 = 831889

PARTICIPANTS_CLASSIC = [
  [7046, 4087],
  [35900,31378],
  [851680,851076],
  [948544,947170],
  [181697,172981],
  [831748,831889],
  [30590,114401]
]

PARTICIPANTS_BIOFACTOIDS = [
  [6573,75947],
  [29843,23410],
  [8520,8367],
  [17869,13555],
  [1612,23586],
  [10419,8607],
  [5230,8678],
  [830909,838236],
  [2475,2932],
  [5324,854557],
  [102800311,3091],
  [85451,8536],
  [283677,256006],
  [84844,55818]
]

#import the data
df = pd.read_csv( './data/gene2pubmed', sep='\t', header=0 , names=['tax_id', 'GeneID', 'PubMed_ID'] )

for ppts in PARTICIPANTS_BIOFACTOIDS:
  # Get rows with Gene ID 1/2
  df_GI1 = df[ df.GeneID == ppts[0] ]
  df_GI2 = df[ df.GeneID == ppts[1] ]

  df_merged = pd.merge( df_GI1, df_GI2, on='PubMed_ID', how='inner')

  print( f'GENEIDs {ppts[0]},{ppts[1]}: Cocitations: {df_merged.shape[0]}' )