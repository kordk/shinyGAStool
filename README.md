# shinyGAStool
shinyGAStool: an open-source user-friendly tool for candidate gene association studies implemented in R-shiny

## Required libraries
The following libraries must be installed prior to running the code:
<pre>
library(shiny)
library(shinyFiles)
library(DT)
library(snpStats)
library(compiler)
library(batch)
install.packages("BiocManager")
BiocManager::install("snpStats")
</pre>

## Running the interface
Run the code
 <pre>
cd shinyGAStool
R --vanilla
source('shinyGeneticsApp.R')
shinyApp(ui, server)
</pre>

## Authors
Thomas Hoffmann: Thomas.Hoffmann@ucsf.edu

Kord M. Kober: Kord.Kober@ucsf.edu

## Acknowledgements
Support for this project was provided by the National Cancer Institute (CA233774). Its contents are solely the responsibility of the authors and do not necessarily represent the official views of the National Institute of Health. 
