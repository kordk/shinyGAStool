# shinyGAStool
shinyGAStool: an open-source user-friendly tool for candidate gene association studies implemented in R-shiny

## Citation
Please cite the tool using this website repository and the manuscript:

Hoffmann TJ, Miaskowski C, Kober KM. ShinyGAStool: A user-friendly tool for candidate gene association studies. SoftwareX. 2023;21:101274. doi: 10.1016/j.softx.2022.101274.

## Installation instructions

### From a Release

Periodically, we provide a release of the tool as a stand-alone Windows executable with an installation tool. Please see the <a href="https://github.com/kordk/shinyGAStool/releases">releases</a> to download the latest version.

### From Source Code

#### Download the app

Make sure you download the following files & folders, making sure that the files in the data/ folder are put in a data/ subfolder where the main application R source file is (shinyGeneticsApp.R, located in the *src* directory on the repository, but download the files to the following directory structure):

- src/shinyGeneticsApp.R
- data/anno_ccds_hg19.txt.gz  
- data/anno_gencode_attrs_hg38.txt.gz
- data/anno_gencode_attrs_hg19.txt.gz
- data/anno_gencode_basic_hg19.txt.gz
- data/anno_refFlat_hg19.txt.gz
- data/anno_ccds_hg38.txt.gz
- data/anno_gencode_basic_hg38.txt.gz
- data/anno_refFlat_hg38.txt.gz

You can optionally (recommended) download the demo files from (and recommended to put in a demo/ subfolder where the main application source file is):

- demo/kgp-eur.bed
- demo/kgp-eur.bim
- demo/kgp-eur.fam
- demo/kgp-eur-ldl-pheno.csv


#### Installing R and/or RStudio Desktop

You can install R from:  
https://cran.r-project.org/

Or install RStudio for a slightly more polished user interface to R:  
https://rstudio.com/products/rstudio/


#### Installing R packages

After launching the newly installed R, from the command prompt, type the following to install necessary dependencies:

    install.packages("shiny")
    install.packages("shinyFiles")
    install.packages("DT")
    install.packages("compiler")
    install.packages("heatmaply")
    install.packages("rio")

    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("snpStats")

#### Launching the app

Navigate to the folder that the shiny app (shinyGeneticsApp.R) is located, and type:

    source("shinyGeneticsApp.R")
    
This should launch the interface. A demo dataset is included in the demo directory.

Note: If you are running the shiny app from within RStudio, make sure to set the working directory:
Session -> Set Working Directory -> To Source File Location

## Authors

Thomas Hoffmann: Thomas.Hoffmann@ucsf.edu

Kord M. Kober: Kord.Kober@ucsf.edu

## Acknowledgements
Support for this project was provided by the National Cancer Institute (CA233774). Its contents are solely the responsibility of the authors and do not necessarily represent the official views of the National Institute of Health. 

