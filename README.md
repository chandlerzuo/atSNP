atSNP
======

OVERVIEW
--------

atSNP ( *Affinity Test for regulatory SNP detection* ) package is a bioinformatics tool for computing and testing motif-SNP interaction. It provides three main functions:

- Compute the binding affinity scores for both reference and SNP allele based on position weight matrices;

- Compute the p-value of the affinity score for each allele;

- Compute the p-value of the affinity score changes between the reference and SNP allele.

atSNP implements the importance sampling algorithm to compute the p-values. Existing tools, such as FIMO and is-rSNP, computes the p-values analytically. This is computationally intensive because the probability sample space is exponential of the motif length. By implementing the importance sampling algorithm, atSNP is able to compute the p-values at incredible speeds while maintaining accuracy.

In one of our research projects, we have used atSNP to evaluate interactions between 20K SNPs and 2K motifs within 5 hours. To our knowledge, no other existing tool can finish the analysis of such a scale.

INSTALLATION
------------

atSNP will be available at Bioconductor. Currently you can download the development version here and install in R by:

    library( devtools )
    install_github( "chandlerzuo/atSNP" )


REFERENCES
----------

Chandler Zuo, Sunyoung Shin and Sunduz Keles (2014). "atSNP: affinity test for regulatory SNP detection". *To appear*.
