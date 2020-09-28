---
title: 'ShinyGAStool: A user-friendly tool for candidate gene association studies '
tags:
  - genetic association
  - regression analysis
  - single nucleotide polymorphism
  - R shiny
authors:
  - name: Thomas J. Hoffmann^[Custom footnotes for e.g. denoting who the corresponding author is can be included like this.]
    orcid: XXXX
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
  - name: Christine Miaskowski
    orcid: XXXX
    affiliation: "3,4" # (Multiple affiliations must be quoted)
  - name: Kord M. Kober
    orcid: 0000-0001-9732-3321
    affiliation: "3,4,5" # (Multiple affiliations must be quoted)
affiliations:
 - name: Department of Epidemiology and Biostatistics, University of California, San Francisco, CA
   index: 1
 - name: Institute for Human Genetics, University of California, San Francisco, CA
   index: 2
 - name: School of Nursing, University of California, San Francisco, CA
   index: 3
 - name: Helen Diller Family Comprehensive Cancer Center, University of California, San Francisco, CA
   index: 4
 - name: Bakar Computational Health Sciences Institute, University of California, San Francisco, CA 
   index: 5
date: 28 September 2020
bibliography: paper.bib

---

# Summary

A rapid increase in the development of genotyping and sequencing technologies along with a decrease in the cost of implementation and data collection has led to a large increase in the availability and consequently the number of genetic association studies.`[@Buniello:2019]`

The major advancements in genotyping and sequencing technologies obtained over the past decade along with a decrease in the cost of implementation and data collection has resulted in the collection of enormous amounts of sequencing data from a variety of biological sources. Along with these new datasets, a large variety of research questions can now be addressed which previously were unattainable. These datasets are becoming more accessible and more fields are utilizing these data. For example, the era of personalized health has seen a dramatic increase in genetic association studies.`[@Gwas:2018]` Despite this increase in usage, many of these new research directions require a deep understanding of the new technologies and their applications, the implementation of novel adaptations to existing tools, or the development of entirely new tools. 

# Statement of need 

A major barrier to the broader adaptation of genomics analyses are the relative complexity of the tools needed to access the data and perform the analysis. In reality, not all research groups have the skills or resources to learn or utilize yet another computational interface (e.g., command line to Unix/Linux systems and tools). For example, although the field of nursing research is rapidly growing in its utilization of systems biology,`[@Founds:2018]` data science,`[@Dreisbach:2020]` and genomics `[@Singh:2018; @Conley:2009]` the informatics support systems and training are still lacking.`[@Calzone:2013]`

The need for tools to support current research as well as teach the next generation of scientists is high. This gap in technical skills limits the ability of researchers to perform the analyses. To address this gap, we introduce shinyGAStool, which was developed with ease of access (i.e., open source) and ease of use (i.e., simple interface) in mind, while also providing a powerful analytic approach. This tool is particularly useful for exploratory analysis or candidate gene analyses. We envision the shinyGAStool will be used for undergraduate and graduate training and exploratory candidate gene analyses (e.g., `[@Illi:2012; @Eshragh:2017; @Kober:2016]`).

We present shinyGAStool, a simple and powerful tool that enables the user to access large genome-wide genotype datasets, genomic annotations, and perform a genetic association analysis on their desktop.

# Implementation and available data
ShinyGAStool is written in R (https://www.R-project.org/) and implemented in shiny (https://CRAN.R-project.org/package=shiny), and runs in all major web browsers tested, making the interface comfortable for a wide audience. The package requires the other R packages DT (https://CRAN.R-project.org/package=DT) and snpStats (https://www.bioconductor.org/packages/snpStats) User sample phenotype/metadata are provided as a spreadsheet in a comma-separated value format (CSV) and sample genotypes are provided as binary PLINK files (Figure 1). A user selects the version the human genome reference to use from the annotations provided for the GRCh38/hg19 and GRCh39/hg38 assemblies. The genetic association analysis is performed using either a linear or logistic regression. A user identifies the dependent variable and independent variables from the sample data provided. Loci for evaluation are identified by gene symbol, dbSNP id, and/or chromosomal location. Three gene annotation tables from the UCSC Genome Browser '[@Karolchik:2011; @Rosenbloom:2015]' are provided for symbol selection (i.e., CCDS,`[@Pujar:2018]` RefSeq, `[@OLeary:2016]` and GENCODE `[@Frankish:2019]`). Loci are modeled as doses (i.e., zero, one, or two copies) of the rare allele. Three genetic models are available for testing (i.e., additive, dominant, and recessive).

A release is also provided as a stand-alone binary with installer for Windows using R-portable.

# Acknowledgements

This study was supported by the National Cancer Institute (NCI, CA233774). It's contents are solely the reponsibility of the authors and do not represent the official views of the NIH.

# References
