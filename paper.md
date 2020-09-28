---
title: 'ShinyGAStool: A user-friendly tool for candidate gene association studies '
tags:
  - genetic association
  - regression analysis
  - single nucleotide polymorphism
  - R shiny
authors:
  - name: Thomas J. Hoffmann
    orcid: 0000-0001-6893-4449
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
  - name: Christine Miaskowski
    orcid: 0000-0001-5170-2027
    affiliation: "3,4" # (Multiple affiliations must be quoted)
  - name: Kord M. Kober^[Corresponding author.]
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

Here we present shinyGAStool, a simple and powerful tool that enables the user to access large genome-wide genotype datasets, genomic annotations, and perform a genetic association analysis on their desktop. ShinyGAStool is freely available under a GPL license from https://github.com/kordk/shinyGAStool.

# Statement of need 

A rapid increase in the development of genotyping and sequencing technologies along with a decrease in the cost of implementation and data collection has led to a large increase in the availability and consequently the number of genetic association studies.[@Buniello:2019] The major advancements in genotyping and sequencing technologies obtained over the past decade along with a decrease in the cost of implementation and data collection has resulted in the collection of enormous amounts of sequencing data from a variety of biological sources. Along with these new datasets, a large variety of research questions can now be addressed which previously were unattainable. These datasets are becoming more accessible and more fields are utilizing these data. For example, the era of personalized health has seen a dramatic increase in genetic association studies.[@Gwas:2018] Despite this increase in usage, many of these new research directions require a deep understanding of the new technologies and their applications, the implementation of novel adaptations to existing tools, or the development of entirely new tools. 

A major barrier to the broader adaptation of genomics analyses are the relative complexity of the tools needed to access the data and perform the analysis. In reality, not all research groups have the skills or resources to learn or utilize yet another computational interface (e.g., command line to Unix/Linux systems and tools). For example, although the field of nursing research is rapidly growing in its utilization of systems biology,[@Founds:2018] data science,[@Dreisbach:2020] and genomics [@Singh:2018; @Conley:2009] the informatics support systems and training are still lacking.[@Calzone:2013]

The need for tools to support current research as well as teach the next generation of scientists is high. This gap in technical skills limits the ability of researchers to perform the analyses. To address this gap, we introduce shinyGAStool, which was developed with ease of access (i.e., open source) and ease of use (i.e., simple interface) in mind, while also providing a powerful analytic approach. This tool is particularly useful for exploratory analysis or candidate gene analyses. We envision shinyGAStool will be used for undergraduate and graduate training and exploratory similary candidate gene analyses in the future. 

# Implementation

Following previous studies that have implemented these methods (e.g., [@Illi:2012; @Eshragh:2017; @Kober:2016]), shinyGAStool implements a regression approach to evaluate for genetic associations. ShinyGAStool is written in R (https://www.R-project.org/) and implemented in shiny (https://CRAN.R-project.org/package=shiny), and runs in all major web browsers tested, making the interface comfortable for a wide audience. The package requires the other R packages DT (https://CRAN.R-project.org/package=DT) and snpStats (https://www.bioconductor.org/packages/snpStats) User sample phenotype/metadata are provided as a spreadsheet in a comma-separated value format (CSV) and sample genotypes are provided as binary PLINK files (\autoref{fig:fig1}). A user selects the version the human genome reference to use from the annotations provided for the GRCh38/hg19 and GRCh39/hg38 assemblies. The genetic association analysis is performed using either a linear or logistic regression. A user identifies the dependent variable and independent variables from the sample data provided. Loci for evaluation are identified by gene symbol, dbSNP id, and/or chromosomal location. Three gene annotation tables from the UCSC Genome Browser [@Karolchik:2011; @Rosenbloom:2015] are provided for symbol selection (i.e., CCDS,[@Pujar:2018] RefSeq, [@OLeary:2016] and GENCODE [@Frankish:2019]). Loci are modeled as doses (i.e., zero, one, or two copies) of the rare allele. Three genetic models are available for testing (i.e., additive, dominant, and recessive). A release is also provided as a stand-alone binary with installer for Windows using R-portable. 

A three-step workflow to perform a genetic association analysis. Data files for demonstration and testing are included in the repository (https://github.com/kordk/shinyGAStool/demo).

_**Variable selection.**_ A user will load phenotype (meta) data and identify dependent and independent variables for evaluation (\autoref{fig:fig2}A, simulated LDL phenotype data). The sample information is provided in a CSV formatted file including the sample identifier, the outcome of interest, and any additional co-variates to include in the model (e.g., eigenvalues of ancestry informative markers). The user then identifies the variable used to connect with the genotype data (provided in the next step) and the outcome of interest. Finally, the co-variates that are identified are plotted for reference. Summary statistics are provided.

_**Genotype loci selection.**_ A user selects the appropriate annotation data from either the GRCh37/hg19 or the GRCh38/hg38 human reference genome assemblies provided with shinyGAStool (\autoref{fig:fig2}B). The study genotype data should be in the common PLINK file format (i.e., bed/bim/fam).[@Purcell:2007] The loci for evaluation are then selected by gene, loci name (e.g., dbSNP ‘rs’ number, selectable only from the bim file, to help avoid any typing mistakes), and/or chromosomal region (chromosome and base pair position). When the gene name is provided, the user can select between gene regions (i.e., Transcription start – end, exons, or coding only) and the type of transcript (i.e., all or specific). The option to include loci upstream and downstream of the gene’s named are also provided. Loci can also be selected based on symbol from the CCDS, RefSeq, and/or GENCODE sources. Once the selection is complete, genotypes for these loci are extracted from the user genotype file and summarized in a table (e.g., frequencies). 

_**Genetic association analysis.**_ Finally, the association test is selected and performed (\autoref{fig:fig2}C) (e.g., [@Illi:2012; @Eshragh:2017; @Kober:2016]). A user will select either a linear or logistic regression appropriate to the trait characteristics of their outcome and select the genetic model(s) to evaluate. The genetics are modeled as the dosage of the rare allele (as discussed above). The results of the analysis can then be saved as a CSV file.

# Conclusion
ShinyGAStool is a simple and powerful tool that enables the user to access large genome-wide genotype datasets, genomic annotations, and perform a genetic association analysis on their desktop.

# Acknowledgements

This study was supported by the National Cancer Institute (NCI, CA233774). It's contents are solely the reponsibility of the authors and do not represent the official views of the National Institue of Health (NIH).

![An overview of the data flow in the ShinyGAStool.\label{fig:fig1}](fig1.tiff)

![Screenshots of the ShinyGAStool workflow. (A) First, the sample characteristics data are loaded and the sample ID, outcome variable, and covariates are identified. The selected variables are presented in a summary table and figures for exploration. (B) The second step is to select the genome build and identify the loci to extract and evaluate. (C) The final step is to select the regression model and genetic model(s) to perform.\label{fig:fig2}](fig2.tiff)

# References
