# Biofactoid: Article Hit Enrichment

## Background

A valuable goal is to identify article authors who could be invited to contribute their article to Biofactoid via the homepage. Currently, articles are manually pre-screened by looking at a given journal, issue by issue, article by article. Typically, the article title, abstract and figures are skimmed to determine suitability. This can take as little as a few seconds to reject whereas it may require tens of minutes to flag an article as a hit - containing at least one, low-throughput, experimentally-verified interaction that is compatible with the Biofactoid model.

## Goal

Ideally, it would be desirable to have a way to directly identify articles suitable for Biofactoid - 'hits'. However, a more realistic goal would be to provide a set of articles for pre-screening enriched for hits. In more formal terms, the goal is to increase the precision of article detection relative to manual prescreening.

As described below, the current hit rate is ~10% for a selected subset of journals. This folder contains files intended to support development of means for enriching article hits.


# Journal analysis

|  **ISSN** | **Journal** | **Coverage [vol(iss)]** | **all** | **hits*** | **rate** |
| :--- | :--- | :--- | :--- | :--- | :--- |
|  2211-1247 | Cell Reports | 30(1) - 32(11) | 953 | 109 | 10.3% |
|  1097-4164 | Molecular Cell | 73(1) - 79(6) | 725 | 85 | 10.5% |
|  1549-5477 | Genes & Development | 34(1-2) - 34(17-18) | 93 | 15 | 13.9% |
|  1476-4679 | Nature Cell Biology | 22(4) - 22(9) | 84 | 10 | 10.6% |
|  1083-351X | Journal of Biological Chemistry | 295(31) - 295(37) | 210 | 21 | 9.1% |
|   |  |  |  |  |  |
|  **Total** |  |  | 2065 | 240 | 10.4% |


# Data

## PubMed

- `all_pmids.txt`: A newline-delimited list of PubMed uids including all the articles considered (N=2065)
- `hits_pmids.txt`: A newline-delimited list of PubMed uids for those articles deemed appropriate for inclusion in Biofactoid (N=240)
- `test_pubmed_info.json`: The EUTILS EFETCH response for each PubMed uid in a sample of `all_pmids.txt` (N=1032)
- `train_pubmed_info.json`: The EUTILS EFETCH response for each PubMed uid that are not part of `test_pubmed_info.json` from (N=1033)
- `test_hits_pmids.txt`:  A newline-delimited list of PubMed uids in `test_pubmed_info.json` that are hits (N=117)
- `train_hits_pmids.txt`:  A newline-delimited list of PubMed uids in `train_pubmed_info.json` that are hits (N=123)


## Scripts

- `main.py`: The script used to generate the test and train files
  - createTestTrainPubMedData: Generate `test_pubmed_info.json` and `train_pubmed_info.json` files
  - getTestTrainHits: Generate `test_hits_pmids.txt` and `train_hits_pmids.txt` files


# Evaulation

Model as a hypergeometric distribution. Try to beat 10% to save the screener time.

//todo

---

Last updated: October 14, 2020