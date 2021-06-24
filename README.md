# Chirona

![Relative mutation score](Relative%20significance%20score.png)

**NOTE: As of 24 Jun 2021, Chriona data has identified the mutation K417N mutation present in the Delta Plus variant (Pango AY.1) with ~24 points. Since this data was taken on 23 May 2021, newer mutations can only be screened with newer data from GISAID/NCBI.**

## Overview

Current SARS-Cov-2 variants and lineages are identified by PANGOLIN.

Emerging variants show resistance against vaccines previously produced, with neutralising effects.

Chirona attempts to find hotspots of mutations in SARS-CoV-2 Pango lineages and generate permutations of sequences that can adapt to said mutations.

## Current mutations of concern as identified by Chirona
| Mutation	| Alignment position	| Score	| Adjusted incidence	| BLOSUM62 delta |
| --- | --- | --- | --- | --- |
| W152C	| 153	| 127.3998575	| 9.79998904	| 13 |
| L452R	| 456	| 121.4137761	| 20.23562935	| 6 |
| N501Y	| 505	| 85.02822072	| 10.62852759	| 8 |
| Y145*	| 145	| 75.00904159	| 6.819003781	| 11 |
| E484K	| 488	| 74.37119842	| 18.59279961	| 4 |
| S13I	| 13	| 59.21420352	| 9.86903392	| 6 |
| P681R	| 685	| 56.39980273	| 6.266644748	| 9 |
| T95I	| 95	| 54.82820976	| 9.138034961	| 6 |

## Credits

Original nucleotide sequences are directly from NCBI's Betacoronavirus database and are last updated here on **23 May 2021 UTC**.

Moderna and Pfizer-BioNTech vaccine mRNA sequences by Dae-Eun Jeong et al., 2021.

Multiple sequence alignments produced solely using MAFFT.

BLOSUM62 JSON by @camillescott.
