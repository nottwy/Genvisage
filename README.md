# *Genvisage* - Rapid Identification of Discriminative Feature Pairs for Genomic Analysis
Silu Huang [shuang86@illinois.edu] Charles Blatti [blatti@illinois.edu], Saurabh Sinha, and Aditya Parameswaran
KnowEnG BD2K Center of Excellence  
University of Illinois Urbana-Champaign  

## Table of Contents
1. [Motivation](#motivation)
2. [Installation](#installation)
3. [Tutorial](#tutorial)
    1. [Creating Gene Sets](#creating-gene-sets)
    2. [Building Heterogeneous Networks](#building-heterogeneous-networks)
    3. [Example DRaWR Runs](#example-drawr-runs)
4. [DRaWR Resources](#drawr-resources)
    1. [Function Parameters](#function-parameters)
    2. [Input File Formats](#input-file-formats)
    3. [Output File Formats](#output-file-formats)
## Motivation

A fundamental, but challenging problem in the analysis of genomics datasets is the separability task: discovering features that most strongly capture the difference between two different classes (sets) of biological objects with high dimensional feature representations, such as gene signatures described by their transcriptomic profiles, or genes described by their functional annotations. This problem is typically addressed by finding discriminative, significant single features using univariate statistical methods, or by identifying difficult-to-interpret multi-feature combinations using time-intensive machine learning algorithms. We sought a middle ground that is simultaneously rapid, interpretable, salient, and significant.

Our tool, Genvisage, enables researchers to interactively identify visually interpretable, significant feature pairs that strongly separate two classes. This webserver embodies one possible use case of Genvisage, where the user wants to find pairs of Gene Ontology functional annotation terms that when considered together are especially good at distinguishing the genes of two separate user-submitted gene sets.

![Method Overview](images/DRaWR_method.small.png)

[Return to TOC](#table-of-contents)

## Installation

### Local copy of DRaWR repository

If you wish to use the sample files necessary to complete this tutorial or the datasets from the paper, first clone this repository from github:
```
git clone https://github.com/cblatti3/DRaWR.git
```

To use the 5 species networks from the paper, you must first unzip their contents:
```
gunzip networks/5ins_cdhmw.names.edge.gz
gunzip networks/5sp_adhiw.names.edge.gz


[Return to TOC](#table-of-contents)

## Tutorial

This section of the README is meant to walk a user through a process of using Genvisage to rapidly identify pairs of interacting functional annotations that are able to discriminate between two distinct gene sets.

### Feature-gene Matrix
The user can use their own feature-gene matrix for analysis, where the value in each cell denotes the feature value for each corresponding gene. Our program can take in different [forms of feature matrix](#forms-feature-matrix). 

User can also use [our provided feature-object matrix](data/feature_gene_scale.selected.txt) with 3,632 Gene Ontology annotation terms as the features and with 22,210 gene objects. Rather than being a 0/1 membership indicator matrix, the features of this matrix represent the diffusion of the gene across a heterogeneous network of prior knowledge about the annotation of and the relationships between genes (see [DRaWR](https://www.ncbi.nlm.nih.gov/pubmed/27153592) method for more details). The prior knowledge in the heterogeneous network used here included annotations from [Gene Ontology](http://www.geneontology.org/), [KEGG](https://www.genome.jp/kegg/), [Reactome](https://reactome.org/), and [Pfam](https://pfam.xfam.org/) as well as gene-gene relationships from protein similarity defined by [BLASTP](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins). After diffusion was performed, 3,632 Gene Ontology terms were extracted creating a term-gene matrix where each value represents the scaled probability that a random walk started at the gene would be at the Gene Ontology term. 

### Creating Gene Sets

The first step is to create two gene sets for discriminative analysis, i.e., positive gene set vs. negative gene set, each stored in one file. The file format should list one gene name on each row. We also allow users to only input one positive gene set. In this instance, Genvisage will use all of the remaining genes in your feature-gene matrix as the negative gene set.   

As an example, uses can also [our provided positive gene set](data/DELYS_THYROID_CANCER), which contains 675 differentially expressed genes between papillary thyroid carcinoma (PTC) compared to normal tissue as the positive gene set. These genes were downloaded from the Molecular Signature Database (MSigDB) ([up regulated](http://software.broadinstitute.org/gsea/msigdb/cards/DELYS_THYROID_CANCER_UP.html), [down regulated](http://software.broadinstitute.org/gsea/msigdb/cards/DELYS_THYROID_CANCER_DN.html)) and the original publication can be found [here](https://www.ncbi.nlm.nih.gov/pubmed/17621275). Since only a positive gene set is provided, the remaining genes in the feature-object matrix will be used as the negative set for the example.

### Selecting Genvisage Algorithm

The user can select different modes of the Genvisage algorithm to mine the underlying feature-object matrix for pairs of Gene Ontology terms whose gene scores separate the positive and negative classes of genes. All of these methods make use of the Rocchio-based separability metric and a transformation of the underlying feature-object matrix for greater computational efficiency. The user can select one of the three algorithmic modes (briefly described below) depending on their desire for greater accuracy or speed in finding the top feature pairs.

#### EarlyOrdering Mode
This method is the most accurate method for identifying top ranking pairs of separating features. It evaluates all possible feature pair candidates using the most accurate scoring that considers all possible genes. It has optimization over baseline methods by maintaining an upper bound for the separability error of the top-k feature pairs and terminating early when a feature pair's error exceeds this upper bound. We also further enhance this mode by scoring the "problematic" genes first.

#### SampOpt Mode
This method is faster that the EarlyOrdering mode because rather than fully evaluating each feature pair, we first calculate the confidence interval of its separability accuracy by performing a faster evaluation on a subsampled set of genes rather than all genes. We can extract the likely top feature pair candidates by comparing the confidence intervals and only perform the full accuracy evaluation for those feature pairs. We further enhanced this mode by ordering the candidate feature pairs which enables the early termination possibility of not having to evaluate all candidates before finding the top-k. This mode is less accurate than EarlyOrdering because true top feature pairs may not be selected ascandidates from their evaluation on the sampled genes.

#### HorizSampOpt Mode
This mode is the fastest method because it rather than sampling and evaluating all possible feature pairs, it greedily only examines 500,000 feature pairs for possible candidates. These feature pairs will be selected by having at least one of their two features being a good single feature separator of the positive and negative gene sets in its own right. This mode is the least accurate because not only might true top feature pairs not be identified as candidates, they might not be considered at all.



[Return to TOC](#table-of-contents)