atSNP
======

OVERVIEW
--------

atSNP ( *Affinity Test for regulatory SNP detection* ) package is a bioinformatics tool for computing and testing large-scale motif-SNP interactions. It provides three main functions:

- Compute the binding affinity scores for both the reference and the SNP alleles based on position weight matrices;

- Compute the p-values of the affinity scores for each allele;

- Compute the p-values of the affinity score changes between the reference and the SNP alleles.

atSNP implements the importance sampling algorithm to compute the p-values. Existing tools, such as FIMO and is-rSNP, compute the p-values analytically. This is computationally intensive because the probability sample space is a exponential order of the motif length. By implementing the importance sampling algorithm, atSNP is able to evaluate the p-value without exhausting the sample space, thereby significantly reduces running time.

In one of our research projects, we have used atSNP to evaluate interactions between 26K SNPs and 2K motifs within 5 hours. To our knowledge, no other existing tool can finish the analysis of such a scale.

INSTALLATION
------------

atSNP will be available at Bioconductor. Currently you can download the development version here and install in R by:

    library(devtools)
    install_github("chandlerzuo/atSNP")


REFERENCES
----------

Chandler Zuo, Sunyoung Shin and Sunduz Keles (2014). "atSNP: affinity test for regulatory SNP detection". *To appear*.
