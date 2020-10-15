# Biofactoid: Article Hit Enrichment

This folder contains a set of data files intended to facilitate automated detection of articles that are suitable for Biofactoid. The goal is to identify article authors who could be invited to contribute their article to Biofactoid via the homepage.

# Journal overview

|  **issn** | **journal** | **coverage [vol(iss)]** | **all** | **hit** | **rate** |
| :--- | :--- | :--- | :--- | :--- | :--- |
|  2211-1247 | Cell Reports | 30(1) - 32(11) | 953 | 109 | 10.3% |
|  1097-4164 | Molecular Cell | 73(1) - 79(6) | 725 | 85 | 10.5% |
|  1549-5477 | Genes & Development | 34(1-2) - 34(17-18) | 93 | 15 | 13.9% |
|  1476-4679 | Nature Cell Biology | 22(4) - 22(9) | 84 | 10 | 10.6% |
|  1083-351X | Journal of Biological Chemistry | 295(31) - 295(37) | 210 | 21 | 9.1% |
|   |  |  |  |  |  |
|  **Total** |  |  | 2065 | 240 | 10.4% |

* 'hit' defined as an article containing at least one experimentally-verified interaction that is supported by the Biofactoid model

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

---

Last updated: October 14, 2020