# Biofactoid: Article Hit Enrichment

## Background

A valuable goal is to identify article authors who could be invited to contribute their article to Biofactoid via the homepage. Currently, articles are manually pre-screened by looking at a given journal, issue by issue, article by article. Typically, the article title, abstract and figures are skimmed to determine suitability. This can take as little as a few seconds to reject whereas it may require tens of minutes to flag an article as a hit - containing at least one, low-throughput, experimentally-verified interaction that is compatible with the Biofactoid model.

## Goal

Ideally, it would be desirable to have a way to directly identify articles suitable for Biofactoid - 'hits'. However, a more realistic goal would be to provide a set of articles for pre-screening enriched for hits. In more formal terms, the goal is to increase the precision of article detection relative to manual prescreening.

As described below, the current hit rate is ~10% for a selected subset of journals. This folder contains files intended to support development of means for enriching article hits.

<br/>

# Journal analysis

|  **ISSN** | **Journal** | **Coverage [vol(iss)]** | **all** | **hit** | **Rate** |
| :--- | :--- | :--- | :--- | :--- | :--- |
|  2211-1247 | Cell Reports | 30(1) - 32(11) | 953 | 109 | 10.3% |
|  1097-4164 | Molecular Cell | 73(1) - 79(6) | 725 | 85 | 10.5% |
|  1549-5477 | Genes & Development | 34(1-2) - 34(17-18) | 93 | 15 | 13.9% |
|  1476-4679 | Nature Cell Biology | 22(4) - 22(9) | 84 | 10 | 10.6% |
|  1083-351X | Journal of Biological Chemistry | 295(31) - 295(37) | 210 | 21 | 9.1% |
|   |  |  |  |  |  |
|  **Total** |  |  | 2065 | 240 | - |
|  **Average** |  |  |  |  | 11.6% |
|  **Weighted** |  |  |  |  | 10.4% |

<br/>

## Data

- `all_pmids.txt`: A newline-delimited list of PubMed uids including all the articles considered (N=2065).
  - Gathered using NCBI EUTILS ESEARCH, filtering for journal ([ta]), volume ([vi]), issue ([ip]) and publication type = Journal Artice([pt])
- `hits_pmids.txt`: A newline-delimited list of PubMed uids for those articles deemed appropriate for inclusion in Biofactoid (N=240)
- `train_set.txt` and `test_set.txt`: Tab-delimited files with the following columns
  - `pmid` The PubMed uid
  - `text` The title and abstract concatenated
  - `hit` Indicator variable for 'hit' or 'miss' with respect to suitability for Biofactoid

## Scripts

- `main.py`: The script used to generate the test and train files

<br/>

# Evaluate

- The null model is a hypergeometric distribution `F(x; N, n, M)` where:
  - M: The total number of articles in the test set
  - n: The total number of hits in the test set
  - N: The number of articles proposed (i.e. by the candidate algorithm)
  - x: The actual number of hits among the set of proposed articles

## Scripts

- `evaluate.py`: This script calculates the cumulative distribution function (`F(x)`) and outputs the value `1 - F(x)`. It also plots the CDF and probability distribution (`f(x)`).
  - Useage: Simply provide a file path containing a newline-delimited list of proposed PubMed uids
    - `python eval.py <path to file>`
  - Example: `python eval.py evalute.txt`
    - Output:
      ```
        M:1032; N:100; n:117; x:23
        x: 23
        F(x): 0.9998821316509936
        1 - F(X): 0.00011786834900640031
      ```
      - ![sample distributions](Figure_1.png)


---

Last updated: October 22, 2020