## MAAMOUL: A method for detecting microbiome-metabolome alterations in disease using metabolic networks

**Table of contents:**
 - [Method overview](#ch1)
 - [Installation](#ch2)
 - [Instructions - Running MAAMOUL on your own data](#ch3)
 - [Usage example](#ch4)

<a id="ch1"></a>
## Method overview

MAAMOUL is a knowledge-based computational method that integrates metagenomic and metabolomic data to identify custom data-driven microbial metabolic modules associated with disease states. Unlike traditional statistical approaches, MAAMOUL leverages prior biological knowledge about bacterial metabolism to link genes to metabolites through a global, microbiome-wide metabolic network, and then projects genes' and metabolites' disease-association scores onto this network. The identified 'modules' are sub-networks in this graph that are significantly enriched with disease-associated features, both metagenomic and metabolomic.

For further details see: Muller E, Baum S, and Borenstein E. __"Detecting Microbiome-Metabolome Alterations in Disease Using Metabolic Networks."__ _In preparation_ 

<img src="docs/wiki_figure.png" width="700">

***

<a id="ch2"></a>
## Installation

MAAMOUL can be installed directly from GitHub, by running the following:

```
install.packages("devtools")  
library(devtools)   
install_github("efratmuller/MAAMOUL")   
library(MAAMOUL)
```

***
   
<a id="ch3"></a>
## Instructions - Running MAAMOUL on your own data

_Coming soon..._

***

<a id="ch4"></a>
## Usage example

```
library(MAAMOUL)
write_test_files()
maamoul(global_network_edges = 'test_input/enzyme_compound_edges_kegg.csv',
  ec_pvals = 'test_input/ec_pvals.tsv',
  metabolite_pvals = 'test_input/mtb_pvals.tsv',
  out_dir = 'test_outputs',
  N_REPEATS = 100,
  N_VAL_PERM = 9,
  N_THREADS = 4
)
```

*** 

For questions about the pipeline, please open an issue (https://github.com/efratmuller/MAAMOUL/issues) or contact Prof. Elhanan Borenstein at elbo@tauex.tau.ac.il.

***
