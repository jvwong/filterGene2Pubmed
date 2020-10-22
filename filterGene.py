import pandas as pd
from Bio import Entrez
import pprint
pp = pprint.PrettyPrinter()
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('max_colwidth', None)

######################################################################################################
########################     ESEARCH          ########################################################
######################################################################################################

# Entrez.email = "jvwong@outlook.com"
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

# def main():
#   ids = []
#   for volume in volumes:
#     for issue in issues:
#       record = esearch( journal, volume, issue )
#       ids += record['IdList']

#   with open('{issn}.txt'.format(issn=journal), 'a') as f:
#     for item in ids:
#       f.write('{item}\n'.format(item=item))


# main()

######################################################################################################
########################     gene_info    ############################################################
######################################################################################################

df = pd.read_csv( './data/gene_info', sep='\t', header=0 , names=[
  'tax_id',
  'GeneID',
  'Symbol',
  'LocusTag',
  'Synonyms',
  'dbXrefs',
  'chromosome',
  'map_location',
  'description',
  'type_of_gene',
  'Symbol_from_nomenclature_authority',
  'Full_name_from_nomenclature_authority',
  'Nomenclature_status',
  'Other_designations',
  'Modification_date',
  'Feature_type'
  ])

# PMID = 9518
taxid = 9606
GeneID = 9518

df_GI = df[ df.GeneID == GeneID ]
# df_GI_index = list( df[ df.tax_id == taxid ].index )
# print( df.loc[df_GI_index, ['tax_id', 'GeneID', 'Symbol', 'description', 'dbXrefs', 'type_of_gene']] )



# ######################################################################################################
# ########################     gene2pubmed    ##########################################################
# ######################################################################################################

# Entrez.email = "jvwong@outlook.com"
# PMID = 9887103

# # PARTICIPANTS_CLASSIC = [
# #     [177393, 175071],
# #   [7046, 4087],
# #   [7046, 4087],
# #   [35900,31738],
# #   [851680,851076],
# #   [948544,947170],
# #   [181697,172981],
# #   [831748,831889],
# #   [30590,114401]
# # ]

# PARTICIPANTS_BIOFACTOIDS = [
# #   [6573,75947],
#   [29843,23410],
# #   [8520,8367],
# #   [17869,13555],
#   # [1612,23586],
# #   [10419,8607],
# #   [5230,8678],
# #   [830909,838236],
# #   [2475,2932],
# #   [5324,854557],
# #   [102800311,3091],
#   # [85451,8536],
# #   [283677,256006],
# #   [84844,55818]
# ]

# ### import the data
# df = pd.read_csv( './data/gene2pubmed', sep='\t', header=0 , names=['tax_id', 'GeneID', 'PubMed_ID'] )

# ## Get all GeneID for given PubMed_ID
# # df_PMID = df[ df.PubMed_ID == PMID ]
# # print( df_PMID )

# # def pubmed_text( pmid ):

# #   handle = Entrez.efetch( db='pubmed', id=pmid, retmode='xml' )
# #   record = Entrez.read( handle )
# #   handle.close()
# #   MedlineCitation = record['PubmedArticle'][0]['MedlineCitation']
# #   PMID = MedlineCitation['PMID']
# #   Article = MedlineCitation['Article']
# #   ArticleTitle = Article['ArticleTitle']
# #   abstr = ''
# #   if 'Abstract' in Article:
# #     Abstract = Article['Abstract']
# #     AbstractText = Abstract['AbstractText']
# #     abstr = ''.join( AbstractText )
# #   return ''.join( [ PMID, ': ', ArticleTitle, abstr ] )

# # pmtext = pubmed_text( '32259912' )
# # print( f'{pmtext}' )

# for ppts in PARTICIPANTS_BIOFACTOIDS:
#   # Get rows with Gene ID 1/2
#   df_GI1 = df[ df.GeneID == ppts[0] ]
#   df_GI2 = df[ df.GeneID == ppts[1] ]
#   df_merged = pd.merge( df_GI1, df_GI2, on='PubMed_ID', how='inner')
#   # raw_pmids = df_merged['PubMed_ID'].values

#   # abstracts = []
#   # pmids = [ str(i) for i in raw_pmids ]
#   # pmids = pmids[60:]
#   # for pmid in pmids:
#   #   abstract = pubmed_text( pmid )
#   #   print( f'\'\'\'{abstract}\'\'\',' )

#   ### Intersection of GENEIDs
#   print( f'GENEIDs {ppts[0]},{ppts[1]}: Cocitations: {df_merged.shape[0]}' )
#   print( list(df_merged.PubMed_ID) )

#   ### Individual participants
#   print( f'GENEID {ppts[0]}: PMIDS: {df_GI1.shape[0]}')
#   print( list(df_GI1.PubMed_ID) )

#   print( f'GENEID {ppts[1]}: PMIDS: {df_GI2.shape[0]}')
#   print( list(df_GI2.PubMed_ID) )

