# shinyGeneticsApp.R
# Version: 2022-05-12
# Authors: Thomas Hoffmann and Kord Kober
# License: GPL v3

# Try LDLR with 10K boundaries

# Required libraries
library(shiny)
library(shinyFiles)
library(DT)
library(snpStats)
library(compiler)
library(heatmaply)
# library(plotly) #?
library(rio)

#library(batch)
#USE_INTERNAL_BROWSER = 0
#parseCommandArgs()

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("snpStats")


# Remaining issues:
# (- Hardening for bad input - done somewhat. Probably could have error messages still be more informative.)

# Maybe future things:
# - Determination of genotype boundaries code is better now, could probably make even better [send starts of all exons at once, not just the current transcript]
# - Print summary of phenotype regressed on covariates (?) [would have to be at the base for the pheno-geno summary, where the model is chosen...???]
# - Pathways
# - Allow transformations of phenotypes, covariates
# - Including gene name in output. [Depends how clear this is with multiple transcripts...]
# - maybe HWE (issue is has to be in controls)
# - Allowing less strict input: e.g., case/control as character?
# - Only autosomes are currently allowed...
# - Autodetect if phenotype is dichotomous (logistic regression) or continuous (linear regression)


# Constants 
DURATION = 10 # How long error messages should stay displayed
#W1 = 4 # Width 1 in display - not used anymore (these have to add up to 12)
#W2 = 8 # Width 2 in display - not used anymore

GA_ALLOWABLE_CHROMOSOMES = 1:22


# 2020-05-06 -- NEW, getting the path of this so that you can souce from *anywhere*

library(base)

thisFilePath = function(){
  fname = ''
  
  cArgs = commandArgs(trailingOnly=FALSE)
  m = grep('--file=', cArgs)
  if(length(m) > 0){
    fname = normalizePath(gsub('--file=', '', cArgs[m]))
  }else{
    fname = normalizePath(sys.frames()[[1]]$ofile)
  }
  
  paste0(dirname(fname), '/')
}

SHINY_GAS_TOOL_THIS_FILE_PATH <<- thisFilePath()

print(SHINY_GAS_TOOL_THIS_FILE_PATH)

SHINY_GAS_TOOL_THIS_FILE_PATH <<- gsub("src/$", "", SHINY_GAS_TOOL_THIS_FILE_PATH)
SHINY_GAS_TOOL_THIS_FILE_PATH <<- gsub("/src$", "", SHINY_GAS_TOOL_THIS_FILE_PATH)

print(SHINY_GAS_TOOL_THIS_FILE_PATH)

# Extra functions

# Read in a plink bim file (has the chr-pos information) into a data.frame
read.bim = function(bimfile){
  data = read.table(bimfile, header=FALSE, stringsAsFactors=FALSE, colClasses=c("integer", "character", "numeric", "integer", "character", "character"))
  names(data) = c("chr", "rs", "morgan", "pos", "a1", "a2")
  return(data)
}
# Write a plink bim file
write.bim = function(bim, bimfile)
  write.table(bim, bimfile, row.names=FALSE, col.names=FALSE, quote=FALSE)
read.map = function(mapfile){
  data = read.table(mapfile, header=FALSE, stringsAsFactors=FALSE)
  names(data)[1:4] = c("chr", "rs", "morgan", "pos")
  return(data)
}
# This is used to handle issues with X&Y chromosomes
fix.bim = function(bim){
  # Use this if we were to test them...
  # bim$chr[bim$chr == 23] = 'X'
  # bim$chr[bim$chr == 24] = 'Y'
  # bim$chr[bim$chr == 25] = 'X' # XY
  
  # ... but for now just eliminate them...
  bim = bim[is.element(bim$chr, GA_ALLOWABLE_CHROMOSOMES), ]
  
  bim
}

# Load in the database of genes. It's not really refFlat anymore, but the name persists in the code
# build = 'hg19' or 'hg38'
loadRefFlat = function(build){
  PATH = SHINY_GAS_TOOL_THIS_FILE_PATH
  
  # refFlat = read.table(paste0('data/', build, '/refFlat.txt.gz'), stringsAsFactors=FALSE)
  # names(refFlat) = c('geneName','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds')
  # refFlat$chrom = gsub('chr', '', refFlat$chrom)
  # refFlat$chrom = gsub('_.*', '', refFlat$chrom)
  # refFlat$chrom[refFlat$chrom == 'X'] = 23
  # refFlat$chrom[refFlat$chrom == 'Y'] = 24
  # 
  # refFlat = refFlat[!duplicated(refFlat$geneName), ] ## THIS IS A REAL PROBLEM I DON'T KNOW HOW TO DEAL WITH YET...
  # 
  # return(refFlat)
  
  refFlat = read.csv(paste0(PATH, 'data/', build, '/ccds.csv.gz'), stringsAsFactors=FALSE) ## id, chr, posStart, posEnd, plusStrand
  refFlat = refFlat[is.element(refFlat$chr, GA_ALLOWABLE_CHROMOSOMES), ]
  return(refFlat)
}

# This is **ONLY** for the developer to run to download necessary database files...
downloadGeneDB = function(){
  library(utils)

  # CCDS (Kord prefers this for the coding regions)  
  download.file('ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/15/CCDS.current.txt' , 'anno_ccds_hg19.txt')
  download.file('ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/23/BuildInfo.current.txt' , 'anno_ccds_hg38.txt')
  
  # RefSeq
  download.file('http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz' , 'anno_refFlat_hg19.txt.gz')
  download.file('http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz' , 'anno_refFlat_hg38.txt.gz')
  
  # Gencode seems to require two files
  download.file('http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeAttrsV33lift37.txt.gz', 'anno_gencode_attrs_hg19.txt.gz')
  download.file('http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV33lift37.txt.gz', 'anno_gencode_basic_hg19.txt.gz')
  download.file('http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeAttrsV33.txt.gz', 'anno_gencode_attrs_hg38.txt.gz')
  download.file('http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV33.txt.gz', 'anno_gencode_basic_hg38.txt.gz')

  # Compress if we can (on a linux machine)
  try(system('gzip anno_ccds_hg19.txt'))
  try(system('gzip anno_ccds_hg38.txt'))
}

strrep = function(v, from_eq_to){
  m = match(v, names(from_eq_to))
  mna = which(!is.na(m))
  v[mna] = from_eq_to[m[mna]]
  v
}
#strrep(1:10, c('1'='foobar', '3'='barfoo'))

# Loads in gene dbs for current build
loadGeneDB = function(build){
  PATH = SHINY_GAS_TOOL_THIS_FILE_PATH
  
  ## ccds$ id, gene, plusStrand, cdsStart, cdsEnd, exonStarts, exonEnds
  ccds = read.table(paste0(PATH, 'data/anno_ccds_', build, '.txt.gz'), header=TRUE, sep='\t', comment.char='', stringsAsFactors=FALSE)
  # table(ccds$ccds_status)
  ccds = ccds[ccds$ccds_status=='Public', ]
  ccds$id = paste0(ccds$gene, ':', ccds$ccds_id)
  ccds$plusStrand = as.numeric(ccds$cds_strand=='+')
  ccds$cds_locations = gsub('[[]|[]]', '', ccds$cds_locations)
  # ccds$exonStarts = ''
  # ccds$exonEnds = ''
  # for(i in 1:nrow(ccds)){ ## ***THIS IS REALLY SLOW...***
  #   toks = strsplit(ccds$cds_locations[i], ', ')[[1]]
  #   ttoks = strsplit(toks, '-')
  #   m = matrix(unlist(ttoks), ncol=2, byrow=TRUE)
  #   ccds$exonStarts[i] = paste(m[,1], collapse=',')
  #   ccds$exonEnds[i] = paste(m[,2], collapse=',')
  # }
  # - Regular expressions to the rescue! Very fast!
  ccds$exonStarts = gsub('[-][0-9]*[,][ ]', ',', paste0(ccds$cds_locations,', '))
  ccds$exonStarts = gsub(',$', '', ccds$exonStarts)
  ccds$exonEnds = gsub('[0-9]*[-]', '', ccds$cds_locations)
  ccds$exonEnds = gsub(' ', '', ccds$exonEnds, fixed=TRUE)
  ccds = ccds[, c('X.chromosome', 'id', 'gene', 'plusStrand', 'cds_from', 'cds_to', 'exonStarts', 'exonEnds')]
  names(ccds) = strrep(names(ccds), c('X.chromosome'='chr', cds_from='cdsStart', cds_to='cdsEnd'))
  ccds$cdsStart = as.numeric(ccds$cdsStart)
  ccds$cdsEnd = as.numeric(ccds$cdsEnd)
  ccds$txStart = ccds$cdsStart  ## prevent any crashes
  ccds$txEnd = ccds$cdsEnd      ## prevent any crashes
  #print("Loaded CCDS")
  
  ## refseq: id, gene, plusStrand, cdsStart, cdsEnd, txStart, txEnd, exonStarts, exonEnds
  #refseq = read.table(paste0('anno_refFlat_', build, '.txt.gz'), header=FALSE, sep='\t', comment.char='', stringsAsFactors=FALSE)
  refseq = read.table(paste0(PATH, 'data/anno_refFlat_', build, '.txt.gz'), header=FALSE, sep='\t', comment.char='', stringsAsFactors=FALSE, colClasses=c('character', 'character', 'character', 'character', 'integer', 'integer', 'integer', 'integer', 'integer', 'character', 'character'))
  names(refseq) = c('gene','id','chr','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds')
  refseq$plusStrand = as.numeric(refseq$strand=='+')
  refseq = refseq[, setdiff(names(refseq), 'exonCount')]
  refseq$exonStarts = gsub(',$', '', refseq$exonStarts)
  refseq$exonEnds = gsub(',$', '', refseq$exonEnds)
  refseq$id = paste0(refseq$gene, ':', refseq$id)
  refseq$chr = gsub('chr', '', refseq$chr)
  #print("Loaded Refseq")
  
  ## gencode:
  #g1 = read.table(paste0('anno_gencode_attrs_', build, '.txt.gz'), header=FALSE, sep='\t', comment.char='', stringsAsFactors=FALSE)
  g1 = read.table(paste0(PATH, 'data/anno_gencode_attrs_', build, '.txt.gz'), header=FALSE, sep='\t', comment.char='', colClasses=c('character', 'character', 'character', 'logical', 'character', 'character', 'character', 'logical', 'character', 'character', 'character', 'integer', 'character', 'character'))
  names(g1) = c(
    'geneid', 'gene', 'gtype', 'gstatus', 'transcriptId',
    'id', #transcriptName
    'ttype', 'tstatus', 'havanaGene', 'havanaTranscript', 'ccdsId', 'level', 'transcriptClass')
  g1 = g1[, c('gene', 'id', 'transcriptId')]
  g2 = read.table(paste0(PATH, 'data/anno_gencode_basic_', build, '.txt.gz'), header=FALSE, sep='\t', comment.char='')
  names(g2) = c(
    'unk', 'transcriptId', 'chr', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'unk1',
    'id',
    'unk2', 'unk3', 'unk4'
  )
  g2$plusStrand = as.numeric(g2$strand == '+')
  #encode = cbind(g1, g2[match(g1$id, g2$id), ])
  #encode = cbind(g1, g2[match(g1$transcriptId, g2$transcriptId), ])
  encode = cbind(g2[, c('chr', 'plusStrand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonStarts', 'exonEnds')], g1[match(g2$transcriptId, g1$transcriptId), c('gene', 'id')])
  encode$exonStarts = gsub(',$', '', encode$exonStarts)
  encode$exonEnds = gsub(',$', '', encode$exonEnds)
  encode$chr = gsub('chr', '', encode$chr)
  #print("Loaded gencode")
  
  # New, allowable chromosomes
  ccds   = ccds[  is.element(ccds$chr,   GA_ALLOWABLE_CHROMOSOMES), ]
  refseq = refseq[is.element(refseq$chr, GA_ALLOWABLE_CHROMOSOMES), ]
  encode = encode[is.element(encode$chr, GA_ALLOWABLE_CHROMOSOMES), ]
  
  genodef = list(ccds=ccds, refseq=refseq, gencode=encode)
  assign("GA_genodef", genodef, envir=.GlobalEnv)

  return(invisible())
}
#loadGeneDB('hg19')

## No benefit to compiling this function
#loadGeneDBC = cmpfun(loadGeneDB)
#loadGeneDBC('hg19')

#loadPlinkGeno = function(genoPre, subset=''){
#  
#}
#tmp = loadGeno('kgmp-eur-19.')





# Assumes that bim is sorted!
calculateChrIndex = function(bim){
  w = which(!duplicated(bim$chr))
  w = c(w, nrow(bim)+1) # explicitly for markOverlap(...) routine
  list(w=w, wchr=bim$chr[w], chr=bim$chr, pos=bim$pos, n=nrow(bim))
}

markOverlap = function(chrIndex, chr, starts, ends){
  o = order(starts)
  starts = starts[o]
  ends = ends[o]
  J = length(starts)
  
  keep = rep(FALSE, chrIndex$n)
  w = which(chr == chrIndex$wchr)
  #cat('w', w, '\n')
  indexStart = chrIndex$w[w]
  indexEnd = chrIndex$w[w+1] - 1
  #cat('indexStart', indexStart, 'indexEnd', indexEnd, '\n')
  i = indexStart # index into the bim
  j = 1  # index into the starts
  while(i<=indexEnd & j<=J){
    #cat('i:', i, '\n')
    #print(keep)
    if(chrIndex$pos[i] >= starts[j]){
      if(chrIndex$pos[i] <= ends[j]){
        keep[i] = TRUE
        i = i + 1
      }else{
        j = j + 1
      }
    }else{
      i = i + 1
    }
  }
  keep
}

markOverlapC = cmpfun(markOverlap)
  
# bim = data.frame(
#   chr=rep(1,10),
#   pos=1:10,
#   stringsAsFactors=FALSE)
# cci = calculateChrIndex(bim)
# cci
# which(markOverlap(cci, 1, c(2,7), c(5,7))) # should give 2, 3, 4, 5, 7
# which(markOverlapC(cci, 1, c(2,7), c(5,7))) # should give 2, 3, 4, 5, 7

fixplinkgeno = function(geno){
  z = NULL
  if(is.matrix(geno)){
    z = matrix(as.numeric(geno), nrow=nrow(geno))
  }else{
    z = as.numeric(geno)
  }
  z[z==0] = NA
  z = z - 1
  z
}
#fixplinkgeno(0:3)
#fixplinkgeno(cbind(0:3, 0:3))


# Shiny functions

# This controls groups of where files are found
getVolumes2 = function(){
  volumes = getVolumes()()
  volumes = c(wd='.', home='~', volumes)
  volumes
}
#getVolumes2()
#getVolumes2 = getVolumes


# Shiny: The layout
ui = fluidPage(
  h1("shinyGAStool: Shiny Genetic Association Tool"),
  #p("Thomas Hoffmann and Kord Kober\n"),
  p(a("Thomas Hoffmann", href="https://profiles.ucsf.edu/thomas.hoffmann"), " and ", a("Kord Kober", href="https://profiles.ucsf.edu/kord.kober")),
  a("Github page", href="https://github.com/kordk/shinyGAStool.git"),

  hr(),
  
  p(strong("Instructions"), ": Load in phenotype data, genotype data, and then run a statistical analysis of the two."),
  p(strong("Please cite:"), em("(In progress)")),
  
  hr(),
  
  h2("Phenotype"),
  
  fluidRow(
    # Interactive input
    #column(
    #  width=W1,
    sidebarPanel(
      shinyFilesButton("phefile", label="Click to select phenotype file...", title="Select phenotype file...", multiple=FALSE),
      textOutput("phefile_selected"),
      selectizeInput("list_id", "ID variable", choices=c(), multiple=!TRUE),
      selectizeInput("list_pheno", "Phenotype (or time for survival)", choices=c(), multiple=!TRUE),
      selectizeInput("list_covar", "Covariate", choices=c(), multiple=TRUE),
      selectizeInput("list_event", "Event (survival only)", choices=c(), multiple=!TRUE)
    ),
    
    # Interactive output
    # Row 1: summary(phefile) [textOutput] [as a data.table, turning summary(...) on it's side?]
    # Row 2: Phenotype distribution, or residual distribution
    # Row 3: Current covariate distribution (or them all?)
    #column(
    #  width=W2, 
    #  offset=W1,
    mainPanel(
      "Start by selecting a phenotype file (csv format, you can export to this from excel) by pressing the first button. Then choose values from that file: the 'ID variable' should be used to match the individual ID in the genotype data (the family ID variable is not used); the 'Phenotype' should be the outcome of interest; the 'Covariate'(s) should be any covariates of interest.",
      #dataTableOutput("phe_summary"), #,
      #plotOutput("plot_phe"),
      #plotOutput("plot_covar")
      uiOutput("phebar") # For drawing distributions [OR, could just be a text box with summary, and the latest distribution...]
      #"4 offset 2"
    )
  ),
  
  hr(),
  
  h2("Genotype"),
  
  fluidRow(
    # Interactive input
    #column(
    #  width=W1,
    sidebarPanel(
      radioButtons("rb_build", "Genotype build", choices=c("GRCh37/hg19", "GRCh38/hg38")),
      shinyFilesButton("genofile", label="Click to select genotype file...", title="Select genotype file...", multiple=FALSE),
      textOutput("genofile_selected"),
      
      
      h4("Select at least one of the following:"),
      
      fluidRow(
        #style = "border: 2px dashed blue;",
        style = "border: 2px solid #bbbbbb;",
        tags$b("GENES"),
        radioButtons("rb_generegion", "Gene region", choices=c("Transcription start - end", "Exons", "Coding only")),
        radioButtons("rb_transcript", "All/Specific transcripts", choices=c("All transcripts", "Specific transcript")),
        fluidRow(
          column(width=6, numericInput("num_gene_upstream", "Upstream:", value=0, min=0)),
          column(width=6, numericInput("num_gene_downstream", "Downstream:", value=0, min=0))),
        tags$b("Three different gene databases can be selected from:"),
        selectizeInput("list_ccds", "CCDS (**Coding only**)", choices=c(), multiple=TRUE),
        selectizeInput("list_refseq", "AND/OR RefSeq", choices=c(), multiple=TRUE),
        selectizeInput("list_gencode", "AND/OR GENCODE", choices=c(), multiple=TRUE)
      ),
      
      # ## DELETE BEGIN ######
      # h4("Select at least one of the following:"),
      # #selectInput("list_pathways", "Pathway", choices=c("KEGG"), multiple=TRUE),
      # #selectizeInput("list_genes", "Genes", choices=c("APOE", "LDLR"), multiple=TRUE),
      # selectizeInput("list_genes", "Genes", choices=c(), multiple=TRUE),
      # fluidRow(
      #   column(width=6, numericInput("num_gene_upstream", "Upstream:", value=0, min=0)),
      #   column(width=6, numericInput("num_gene_downstream", "Downstream:", value=0, min=0))),
      # ## DELETE END ######
      
      
      #selectizeInput("list_snps", "AND/OR SNP names (e.g., rs #)", choices=c("rs111", "rs222"), multiple=TRUE),
      
      fluidRow(
        #style = "border: 2px dashed red;",
        style = "border: 2px solid #bbbbbb;",
        tags$b("AND/OR SNP NAMES:"),
        selectizeInput("list_snps", "SNP names (e.g., rs #)", choices=c(), multiple=TRUE),
      ),
      
      fluidRow(
        #style = "border: 2px dashed green;",
        style = "border: 2px solid #bbbbbb;",
        tags$b("AND/OR BY CHOROMOSOME & POSITION:"),
        textInput("text_chrpos", "e.g., <chr>:<pos1>-<pos2>, <chr>:<pos1>-<pos2>")
        ),
      
      actionButton("button_extract", "Extract genotypes...")
    ),
    
    # Interactive output -> SNP, (Gene), (Pathway), A1, A2, Freq(A1), HWE? [tableOutput, renderTable(...) or renderDataTable(...)]
    # --> Text on how to interpret
    #column(
    #  width=W2,
    #  offset=W1,
    mainPanel(
      "First make sure that you have the correct genotype build. To spot check this, open up, e.g., your .bim file in a text editor, find an rs number, and search for that on dbSNP, and match it to your build.",
      "Second load in your .bim file (plink format), by pressing the 'Click to select genotype file...' button. It is assumed that their is a corresponding .fam and .bed file. This will populate the SNP name selection list.",
      "Third, choose either genes (populated from refSeq.txt from the UCSC genome browser), or SNP names (from your dataset), and then press the 'Extract genotype...' button.",
      #uiOutput("genobar")
      dataTableOutput("geno_summary"),
      #uiOutput("ld") 
      #plotOutput("plot_ld")
      plotlyOutput("plot_ld"),
      a("For LD in other populations, you might consider LD Link", href="https://ldlink.nci.nih.gov/?tab=ldmatrix"),
      renderUI("<br>"),
      a("For a view of gene annotation, you might consider the UCSC genome browser", href="https://genome.ucsc.edu/cgi-bin/hgGateway")
    )
  ),

  hr(),
  
  h2("Phenotype genotype analysis"),
  
  fluidRow(
    # Model selection parameters
    #column(
    #  width=W1,
    sidebarPanel(
      radioButtons("rb_model", "Model", choices=c("Linear regression", "Logistic regression", "Cox PH", "Custom")),
      textInput("text_custom_model", "Custom model, e.g., lm(outcome ~ g + covariate, data=phe), where g represents the SNP.", ""),
      checkboxGroupInput("cb_gmodel", "Genetic model", choices=c("Additive", "Dominant", "Recessive"), selected=c("Additive")),
      actionButton("button_process", "Run analysis..."),
      #actionButton("button_save", "Save analysis results...")
      shinySaveButton("button_save", "Save analysis results...", "Save analysis results as...", filetype=list(csv="csv"))
    ),
    
    # Output data frame [tableOutput]
    #column(
    #  width=W2,
    #  offset=W1,
    mainPanel(
      #uiOutput("phenogenobar")
      "Choose linear regression for a quantitative trait, logistic regression for a dichotomous/binary trait, cox for survival analysis, or enter a custom R code to analyze any model R can by entering the snp as 'g' (this allows for any other regression command, e.g., lmer, glmer, quantile regression, as well as gene-environment interaction). Choose one or more genetic models. Press the 'Run analysis...' button to start the analysis, and 'Save analysis results...' to save those results.\n",
      textOutput("text_model"),
      dataTableOutput("analysis_summary"),
      "Results of the additive, dominant, and recessive models have suffixes of _a, _d, and _r, respectively.",
      plotlyOutput("plot_man"),
      a("For more advanced manhattan (and local manhattan) plots, and a comprehensive look at gwas hits, eqtls, etc., you might consider FUMA.", href="https://fuma.ctglab.nl/"),
      a("For more advanced plots of a locus, you might consider locuszoom", href="http://locuszoom.org/")
    )
  )

)


# Update to rio, now we can handle lots of data formats!
# https://www.rdocumentation.org/packages/rio/versions/0.5.29

# Shiny: Handling interaction
# (Main place to harden for errors...)
server = shinyServer(function(input, output, session){
  #output$phetext = renderPrint({input$phefile})
  
  observe({
    shinyFileChoose(input, 'phefile', root=getVolumes2(), session=session, filetypes=c('csv', 'psv', 'tsv', 'csvy','sas7bdat', 'sav', 'zsav', 'dta', 'xpt', 'por', 'xls', 'xlsx', 'rds', 'RData', 'rda', 'rec', 'mtp', 'syd', 'dbf', 'orff', 'dif', 'fwf', 'csv.gz', 'parquet', 'wft', 'feather', 'fst', 'json', 'mat', 'ods', 'html', 'xml', 'yml', 'pzfx'))
    #shinyFileChoose(input, 'phefile', root=getVolumes(), defaultRoot='wd', session=session)
    #shinyFileChoose(input, phefile, root=normalizePath(getwd()), session=session)
    if(!is.null(input$phefile)){
      fname = parseFilePaths(getVolumes2(), input$phefile)$datapath
      output$phefile_selected = renderText(paste("Phenotype file:", fname))
      
      if(!is.null(fname) && length(fname)>0){
        #print(fname) ## temp
        
        assign("GA_phename", fname, envir=.GlobalEnv)
        
        ## Load in the bim, and then set the possible SNP options
        #assign("GA_phe", NULL)
        tryCatch({
          #phe = read.csv(fname, stringsAsFactors=FALSE)
          phe = import(fname)
          assign("GA_phe", phe, envir=.GlobalEnv)
          
          # Temporary, fill in with actual names of data.frame once load in data
          #phenames = c('LDL', 'HDL', 'Age', 'Sex')
          #phenames = setdiff(names(phe), c('FID', 'IID'))
          phenames = names(phe)
          
          #updateSelectizeInput(session=session, "list_id", choices=phenames)
          #updateSelectizeInput(session=session, "list_pheno", choices=phenames)
          #updateSelectizeInput(session=session, "list_covar", choices=phenames)
          updateSelectizeInput(session=session, "list_id", choices=phenames, server=TRUE)
          updateSelectizeInput(session=session, "list_pheno", choices=phenames, server=TRUE)
          updateSelectizeInput(session=session, "list_covar", choices=phenames, server=TRUE)
          updateSelectizeInput(session=session, "list_event", choices=phenames, server=TRUE)
          
          # Output a summary of the variables...
          pheSum = data.frame(
            Variable = names(phe),
            Class = sapply(phe, 'class'),
            Min=NA, First_Qu=NA, Median=NA, Mean=NA, Third_Qu=NA, Max=NA)
          for(i in 1:ncol(phe)){
            if(pheSum$Class[i] != 'character'){
              pheSum$Min[i]      = min(phe[[i]], na.rm=TRUE)
              pheSum$First_Qu[i] = quantile(phe[[i]], .25, na.rm=TRUE)
              pheSum$Median[i]   = median(phe[[i]], na.rm=TRUE)
              pheSum$Mean[i]     = mean(phe[[i]], na.rm=TRUE)
              pheSum$Third_Qu[i] = quantile(phe[[i]], 0.75, na.rm=TRUE)
              pheSum$Max[i]      = max(phe[[i]], na.rm=TRUE)
            }
          }
          
          # OUTPUT FORMAT
          for(k in c('Min', 'First_Qu', 'Median', 'Mean', 'Third_Qu', 'Max')){
            pheSum[[k]] = sprintf('%0.4g', pheSum[[k]])
          }
        
          output$phe_summary = DT::renderDataTable(DT::datatable(pheSum))
          
          #for(i in phenames){
          #  output$phetext = renderText(paste(phenames, collapse="\n"))
          #}
          ##output$phebar = renderUI({
          ##  tagList(
          ##    sliderInput("n", "N", 1, 1000, 500),
          ##    textInput("label", "Label")
          ##  )
          ##})
          
          #tl = list()
          #for(i in phenames)
          #  tl[[i]] = textInput(phenames[[i]], phenames[[i]])
          #class(tl) = c("shiny.tag.list", "list")
          #output$phebar = renderUI({tagList(tl)})
          
          
          ## New, now expand this here, instead of initially?
          output$phebar = renderUI({
            mainPanel(
              dataTableOutput("phe_summary"),
              plotOutput("plot_phe"),
              plotOutput("plot_covar"),
              textOutput("text_event")
            )
          })
        }, error=function(e){
          showNotification(paste("Error loading csv phenotype file, please check it's format. Error:", as.character(e)), duration=DURATION, session=session, type='error')
        })
      }
    }
  })
  
  observe({
    shinyFileChoose(input, 'genofile', root=getVolumes2(), session=session, filetypes=c('bim'))
    if(!is.null(input$genofile)){
      fname = parseFilePaths(getVolumes2(), input$genofile)$datapath
      
      if(!is.null(fname) && length(fname)>0){
        #print(fname) ## temp
        
        # Remove the file extension (keep the '.')
        fname = gsub('bed$|bim$|fam$', '', fname)
        # Assign it to a global variable
        assign("GA_genofile", fname, envir=.GlobalEnv)

        # Also output it
        output$genofile_selected = renderText(paste("Genotype file:", fname, "[bed/bim/fam]"))

        tryCatch({
          withProgress(message="Loading list of possible SNPs...", value=0, {        
            ## Load in the bim, and then set the possible SNP options
            bim = read.bim(paste0(fname, 'bim'))
            bim = fix.bim(bim) ## A little questionable...
            assign("GA_bim", bim, envir=.GlobalEnv)
            
            incProgress(0.50, detail="Populating list of SNPs...")
    
            #updateSelectizeInput(session=session, "list_snps", choices=bim$rs)
            updateSelectizeInput(session=session, "list_snps", choices=bim$rs, server=TRUE)
          })
        },
        error = function(e){
          showNotification(paste("Error loading in the list of possible SNPs from the bimfile. Error:", as.character(e)), type='error', duration=DURATION, session=session)
        })
      }
    }
  })
  
  observe({
    if(!is.null(input$rb_build) & nchar(input$rb_build)>0){
      #print(input$rb_build) # looks good!!!
      
      tryCatch({
        withProgress(message="Loading list of genes...", value=0, {
          build = tail(strsplit(input$rb_build, '/', fixed=TRUE)[[1]], n=1)
          
          loadGeneDB(build)
          
          #incProgress(0.50, detail="Populating list of genes...") # Done elsewhere, and so fast itdoesn't actually matter...
        })
      },
      error = function(e){
        showNotification(paste("Error loading in the list of genes... Did you download supplementary data files? Error:", as.character(e)), type='error', duration=DURATION, session=session)
      })
    }
  })
  
  #radioButtons("rb_transcript", "All/Specific transcripts", choices=c("All transcripts", "Specific transcript")),
  observe({
    if(!is.null(input$rb_transcript) & nchar(input$rb_transcript)>0){
      tryCatch({
        genodef = get("GA_genodef")
        if(input$rb_transcript == "All transcripts"){
          updateSelectizeInput(session=session, 'list_ccds', choices=unique(genodef$ccds$gene), server=TRUE)
          updateSelectizeInput(session=session, 'list_refseq', choices=unique(genodef$refseq$gene), server=TRUE)
          updateSelectizeInput(session=session, 'list_gencode', choices=unique(genodef$gencode$gene), server=TRUE)
        }else{
          updateSelectizeInput(session=session, 'list_ccds', choices=genodef$ccds$id, server=TRUE)
          updateSelectizeInput(session=session, 'list_refseq', choices=genodef$refseq$id, server=TRUE)
          updateSelectizeInput(session=session, 'list_gencode', choices=genodef$gencode$id, server=TRUE)
        }
      })
    }
  })
  
  
  observeEvent(input$list_pheno, {
    tryCatch({
      #print(input$list_pheno)
      if(!is.null(input$list_pheno) && nchar(input$list_pheno) > 0){
        if(!is.element(class(get("GA_phe")[[input$list_pheno]]), c('numeric', 'integer')))
          stop("Phenotype must be numeric.")
        output$plot_phe = renderPlot({
          hist(get("GA_phe")[[input$list_pheno]], xlab=input$list_pheno, main='', col='pink') # what about # breaks?
        })
      }
    },
    error = function(e){
      showNotification(paste("Error in handling the plotting of the phenotype. If this was a character, it might need to be dichotomous... Error:", as.character(e)), type='error', duration=DURATION, session=session)
    })
  }) # list_covar
  
  observeEvent(input$list_covar, {
    if(!is.null(input$list_covar)){
      tryCatch({
        covars = input$list_covar
        #par()
        #print(covars)
        C = length(covars)
        ncol = ceiling(sqrt(C))
        nrow = ceiling(C/ncol)
        output$plot_covar = renderPlot({
          par(mfrow=c(nrow, ncol))
          for(i in 1:C)
            hist(get("GA_phe")[[covars[i]]], xlab=covars[i], main='', col='light blue')
        }) #, height="1200px")
      },
      error = function(e){
        showNotification(paste("Error in trying to plot the covariates. If any of those were of character formats, it might need to be numeric... Error:", as.character(e)), type='error', duration=DURATION, session=session)
      })
    }
  })

  observeEvent(input$list_event, {
    if(!is.null(input$list_event)){
      tryCatch({
        event = input$list_event
        eventv = get("GA_phe")[[event]]
        outstr = sprintf("Event (vs. censor) variable selected (%s): 0(%i), 1(%i), NA(%i)", event, sum(event==0,na.rm=TRUE), sum(event==1,na.rm=TRUE), sum(is.na(event)))
        output$text_event = renderText(outstr)
      },
      error = function(e){
        showNotification(paste("Error in trying to summarize the event variable, which is expected to have values 0, 1, and NA if missing. Error:", as.character(e)), type='error', duration=DURATION, session=session)
      })
    }
  })
  
  observeEvent(input$button_extract, {
    #print("Load genotypes button pressed...")

    tryCatch({
      # Resolve the list of SNPs that should be loaded
      
      # Start with the list of SNPs a user specificed
      snps = input$list_snps 
      
      # Add in genes [ccds, refseq, gencode]
      upstream = as.numeric(input$num_gene_upstream)
      downstream = as.numeric(input$num_gene_downstream)
      #print(paste("NOT YET HANDLED UPSTREAM/DOWNSTREAM:", upstream, "/", downstream))

      # Get variables that determine how we should extract
      genodef = get("GA_genodef")
      allTranscripts = (input$rb_transcript == "All transcripts")
      #radioButtons("rb_generegion", "Gene region", choices=c("Transcription start - end", "Exons", "Coding only")),
      generegion = input$rb_generegion

      # Get the set of available SNPs in the genotype file
      bim = get("GA_bim")
      bimKeep = rep(FALSE, nrow(bim))
      
      # so plusStrand[i] + 1 indexes properly into this
      listUpstream = c(downstream, upstream)
      listDownstream = c(upstream, downstream)
      
      #print("listUpstream")
      #print(listUpstream)

      withProgress(message="Determining genotype boundaries...", value=0, {
        cci = calculateChrIndex(bim)
        
        # Go through the different gene definitions the user might have kept
        for(defn in c('ccds', 'refseq', 'gencode')){
          # List of genes in the current definition
          genes = input[[paste0('list_', defn)]]
          
          if(length(genes) > 0){
            # Subset out to the definition being used
            sgenodef = NULL
            if(allTranscripts){
              sgenodef = genodef[[defn]][is.element(genodef[[defn]]$gene, genes), ]
            }else{
              sgenodef = genodef[[defn]][is.element(genodef[[defn]]$id, genes), ]
            }
            #print("dim(sgenodef)")
            #print(dim(sgenodef))
            #print(sgenodef)
            
            for(i in 1:nrow(sgenodef)){
              curPlusStrandP1 = sgenodef$plusStrand[i] + 1
              #print("curPlusStrandP1")
              #print(curPlusStrandP1)
              curUpstream = listUpstream[curPlusStrandP1]
              curDownstream = listDownstream[curPlusStrandP1]
              #print("curUpstream")
              #print(curUpstream)
              #print("curDownstream")
              #print(curDownstream)
              
              # Go through each in level of restrictiveness
              ## THIS CODE IS SUPER SLOW... (Somewhat expected...)
              if(generegion == "Transcription start - end"){
                bimKeep[
                  bim$chr == sgenodef$chr[i] &
                    (sgenodef$txStart[i] - curUpstream) <= bim$pos &
                    bim$pos <= (sgenodef$txEnd[i] + curDownstream)
                  ] = TRUE
              }else if(generegion == "Exons"){ # ignores up/downstreams
                starts = as.numeric(unlist(strsplit(sgenodef$exonStarts[i], ',')))
                ends = as.numeric(unlist(strsplit(sgenodef$exonEnds[i], ',')))
                # for(j in 1:length(starts))
                #  bimKeep[
                #    bim$chr == sgenodef$chr[i] &
                #      starts[j] <= bim$pos &
                #      bim$pos <= ends[j]
                #    ] = TRUE
                
                # More optimized, but we really could do even better
                ## cci calculated above
                bimKeep = bimKeep | markOverlapC(cci, sgenodef$chr[i], starts, ends)
              }else if(generegion == "Coding only"){
                starts = as.numeric(unlist(strsplit(sgenodef$exonStarts[i], ',')))
                ends = as.numeric(unlist(strsplit(sgenodef$exonEnds[i], ',')))
                # for(j in 1:length(starts))
                #   bimKeep[
                #     bim$chr == sgenodef$chr[i] &
                #       #starts[j] <= bim$pos &
                #       #bim$pos <= ends[j] &
                #       #sgenodef$cdsStart[i] <= bim$pos &
                #       #bim$pos <= sgenodef$cdsEnd[i]
                #       pmax(starts[j], sgenodef$cdsStart[i]) <= bim$pos &
                #       bim$pos <= pmin(ends[j], sgenodef$cdsEnd[i])
                #     ] = TRUE
                
                # Again, more optimized, but we really could do even better
                ## cci calculated above
                starts = pmax(starts, sgenodef$cdsStart[i])
                ends = pmin(ends, sgenodef$cdsEnd[i])
                bimKeep = bimKeep | markOverlapC(cci, sgenodef$chr[i], starts, ends)
              }else{
                stop("Coding error, generegion has an impossible value!")
              }
            }
          }
        
          snps = unique(c(snps, bim$rs[bimKeep]))
          
          incProgress(0.30, paste0(defn))      
        }
      })
      
      # Add in ranges
      if(nchar(input$text_chrpos) > 0){
        bim = get("GA_bim")
        chrpos = trimws(unlist(strsplit(input$text_chrpos, ',')))
        #print(chrpos) 
        for(cpos in chrpos){
          toks = unlist(strsplit(cpos, ':|-'))
          curChr = toks[1]
          curPos1 = as.numeric(toks[2])
          curPos2 = as.numeric(toks[3])
          
          snps = unique(c(
            snps,
            bim$rs[bim$chr==curChr & curPos1<=bim$pos & bim$pos<=curPos2]
          ))
        }
      }
      
      if(length(snps) < 1)
        stop("No genotypes meet the restrictions in the list...")
      
      # Then load in and summarize the genotype data...
      withProgress(message="Extracting genotypes...", value=0, {
        # Load in the SNPs
        incProgress(0.01, detail="Loading in the genotypes...")
        genofile = get("GA_genofile")
        plinkdata = read.plink(bed=paste0(genofile,'bed'), bim=paste0(genofile,'bim'), fam=paste0(genofile,'fam'), select.subjects=NULL, select.snps=snps) # $genotypes (snpMatrix), $fam, $map
        
        # Calculate the frequency
        incProgress(0.30, detail="Determining if alleles need to be flipped...")
        snpsum = plinkdata$map[, c('snp.name', 'chromosome', 'position', 'allele.1', 'allele.2')]
        snpsum$freq1 = NA
        for(i in 1:nrow(snpsum)){
          #snpsum$freq1[i] = mean(as(plinkdata$genotypes[,i], 'numeric'), na.rm=TRUE) / 2
          snpsum$freq1[i] = mean(fixplinkgeno(plinkdata$genotypes[,i])) / 2  ## 2020-04-30 -- unexpected coding fix
        }
        
        # Flip to coding the minor allele
        if(any(snpsum$freq1 > 0.5))
          plinkdata$genotypes = switch.alleles(plinkdata$genotypes, which(snpsum$freq1>0.5))
    
        # Store it as a global variable
        assign("GA_plinkdata", plinkdata, envir=.GlobalEnv)
  
        incProgress(0.30, "Calculating final summaries...")      
        # Now run the internal summaries
        ## TODO: better to run HWE only on the controls...
        snpsum = plinkdata$map[, c('snp.name', 'chromosome', 'position', 'allele.1', 'allele.2')]
        cs = col.summary(plinkdata$genotypes)
        cs$p.HWE = sprintf("%0.2g", 2 * pnorm(abs(cs$z.HWE), lower.tail=FALSE))
        cs = cs[, setdiff(names(cs), c('MAF', 'z.HWE'))] # Confirmed RAF=MAF
        
        ## REFORMAT COLUMNS
        for(k in c('RAF', 'P.AA', 'P.AB', 'P.BB'))
          cs[[k]] = sprintf('%0.4f', cs[[k]])
        
        snpsum = cbind(snpsum, cs)
        
        incProgress(0.30, "Rendering final output...")
        output$geno_summary = DT::renderDataTable(DT::datatable(snpsum))

        # NEW 2022-05-11
        ld.mydata = ld(plinkdata$genotypes, stats=c("D.prime", "R.squared"), depth=ncol(plinkdata$genotypes)-1)
        #ld.mydata = ld(plinkdata$genotypes, stats=c("R.squared"), depth=nrow(cs)-1)
        spectrum = rainbow(10, start=0, end=1/6)[10:1]
        #output$plot_ld = renderPlot(plot(0:1, 0:1, xlab="ld")) # testing

        m.r2 = as.matrix(ld.mydata$R.squared)
        ltri = lower.tri(m.r2)
        m.r2[ltri] = t(m.r2)[ltri]
        diag(m.r2) = 1
        m.dp = as.matrix(ld.mydata$D.prime)
        m.dp[ltri] = t(m.dp)[ltri]
        diag(m.dp) = 1
        nsnp = nrow(m.r2)
        m.cell = matrix("", nrow=nsnp, ncol=nsnp)
        for(i in 1:nsnp)
          for(j in 1:nsnp)
            m.cell[i,j] = sprintf(
              "Variant 1: %s\nVariant 2: %s\nR^2=%0.3f\nD'=%0.3f",
              rownames(m.r2)[i],
              colnames(m.r2)[j],
              m.r2[i,j],
              m.dp[i,j]
            )
        m = m.r2
        m[lower.tri(m)] = -t(m.dp)[lower.tri(t(m.dp))]
        
        #m.cell = as.matrix(as.vector(sprintf("R^2=%0.2f\nD'=%0.2f", m.r2, m.dp)), nrow=nrow(m.r2))
        #m.cell = as.matrix(paste0("R^2=", sprintfm.r2, '\n', "D'=", m.dp), nrow=nrow(m.r2), ncol=nrow(m.r2))
        #m[lower.tri(m)] = NA
        #m = m[ncol(m):1, ]
        #ncut=100

        rrgb = function(r, g, b)
          rgb(r/255, g/255, b/255)

        #output$plot_ld = renderPlot(
        output$plot_ld = renderPlotly(
          #{
            heatmaply(
              m,
              colors=rrgb(
                  c(0:255,  rep(255, 255)),
                  c(0:255,  254:0),
                  c(rep(255,255),  255:0)
              ),
              custom_hovertext = m.cell,
              dendrogram="none",
              xlab="", ylab="",
              main="D' (lower triangle) \\ R^2 (upper triangle)",
              grid_color="white",
              grid.width=0.00001,
              hide_colorbar=!TRUE,
              labCol = colnames(m),
              labRow = rownames(m),
              plot_method = "plotly"   # controls the custom_hovertext!
            )

            #heatmap(m, Rowv=NA, Colv=NA, col=)

            #image(ld.mydata$R.squared, lwd=0, cuts=9, col.regions=spectrum, colorkey=TRUE, main=expression(R^2), axis=FALSE)


            # library(LDheatmap)
            # LDheatmap(
            #   plinkdata$genotypes,
            #   genetic.distances=plinkdata$map$position,
            #   distances="physical",
            #   LDmeasure="r",
            #   title="Pairwise LD with R^2",
            #   add.map=TRUE, add.key=TRUE,
            #   geneMapLocation=0.15,
            #   SNP.name=rownames(plinkdata$map),
            #   color=NULL
            # )
          #}
        )
      })
    },
    error = function(e){
      showNotification(paste("Error in fully loading in the genotypes. Error:", as.character(e)), type='error', duration=DURATION, session=session)
    })
  })
  
  observeEvent(input$button_process, {
    #print("Button process...")
    tryCatch({
      r_model = input$rb_model
      g_model = input$cb_gmodel
      #print(r_model)  # "Linear regression", "Logistic regression"
      #print(g_model)  # Additive, Dominant, Recessive
  
      #ga = is.element("Additive", g_model)
      #gd = is.element("Dominant", g_model)
      #gr = is.element("Recessive", g_model)
      ggmods = c(
      'a'[is.element("Additive", g_model)],
      'd'[is.element("Dominant", g_model)],
      'r'[is.element("Recessive", g_model)]
      )

      
      rl = (r_model == "Linear regression")
      rlogistic = (r_model == "Logistic regression")
      rcox = (r_model == "Cox PH")
      rcustom = (r_model == "Custom")
      
      pheno = input$list_pheno
      covar = input$list_covar
      event = input$list_event
      
      id_var = input$list_id
      
      phe = get("GA_phe")
      
      if(rl || rlogistic){
        phe = phe[complete.cases(phe[, c(pheno, covar)]), ]
      }else if(rcox){
        phe = phe[complete.cases(phe[, c(pheno, covar, event)]), ]
      }else if(rcustom){
        ## Nothing done!!!
      }
      
      plinkdata = get("GA_plinkdata")
      
      #snpres = plinkdata$map[, c('snp.name', 'chromosome', 'position', 'allele.1', 'allele.2')]
      #snpres$model = NA
      #snpres$gene = NA
      #snpres$beta = NA
      #snpres$se = NA
      #snpres$p = NA
      
      G = nrow(plinkdata$map)
      freq = rep(NA, G)
      
      fitl = list()

      cmd = ''
      pstr = ''
      coefstr = 'Estimate'
      sestr = 'Std. Error'
      if(rl){
        cmd = paste0('lm(', pheno, ' ~ ', paste(c(covar, 'g'), collapse=' + '), ', data=phe)')
        pstr = 'Pr(>|t|)'
      }else if(rlogistic){
        cmd = paste0('glm(', pheno, ' ~ ', paste(c(covar, 'g'), collapse=' + '), ', data=phe, family="binomial")')
        pstr = 'Pr(>|z|)'
      }else if(rcox){
        cmd = paste0('coxph(Surv(', pheno, ',', event, ') ~ ', paste(c(covar, 'g'), collapse=' + '), ', data=phe)')
        pstr = 'Pr(>|z|)'
        coefstr = 'coef'
        sestr = 'se(coef)'
      }else if(rcustom){
        cmd = input$text_custom_model
      }
      output$text_model = renderText(paste0("\nMODEL:\n", paste0("", cmd), "\n")) ## The paste0 looks stupid, but otherwise it does it by reference, and the value changes...
      if(!rcustom)
        updateTextInput(session, "text_custom_model", value=paste0("", cmd))

      # But then actually want this
      cmd = paste0("cf = coef(summary(", cmd, "))")

      ##m = match(phe$IID, plinkdata$fam$member)
      m = match(phe[[id_var]], plinkdata$fam$member)
      withProgress(message="Analyzing each genotype...", value=0, {
        for(i in 1:G){
          incProgress(1/G)

          curfit = c(snp_index=i) # index so that we can merge in despite crashed models. we don't want to keep this in the final output though
          
          # Pull in the next genotype
          #phe$g = as(plinkdata$genotypes[,i], 'numeric')[m]
          gcurrent = fixplinkgeno(plinkdata$genotypes[,i])[m]
          phe$g = gcurrent 
          
          freq[i] = mean(phe$g, na.rm=TRUE) / 2

          try({
            for(ggmod in ggmods){ #c('a', 'd', 'r')){
              # Handle the coding to the genotype
              if(ggmod == 'd'){
                phe$g = as.numeric(gcurrent != 0)
              }else if(ggmod == 'r'){
                phe$g = as.numeric(gcurrent == 1)
              }
              # else it's additive, but nothing to be done

              # Fit the model
              cf = NULL
              #print(cmd) # debug only
              eval(parse(text=cmd)) # debug only
              #print(cf) # debug only

              if(!rcustom){
                ccurfit = c(
                  beta = cf['g', coefstr],
                  se = cf['g', sestr],
                  p = cf['g', pstr]
                )
                names(ccurfit) = paste0(names(ccurfit), '_', ggmod)

                curfit = c(curfit, ccurfit)
              }else{
                # Stick in the entirety of the results...
                for(jj in 1:nrow(cf)){
                  ccurfit = cf[jj,]
                  names(ccurfit) = paste0(rownames(cf)[jj], "_", colnames(cf), "_", ggmod)
                  curfit = c(curfit, ccurfit)
                }
              }# if(!rcustom)
            } # for ggmod

            # done a little specially for crashed fits
            if(length(curfit) > 1)
              fitl[[length(fitl) + 1]] = curfit
          })
        } # for i in 1:G
      })
      
      # Put together the results
      snpres = plinkdata$map[, c('snp.name', 'chromosome', 'position', 'allele.1', 'allele.2')]
      snpres = cbind(snpres, freq.1 = freq)

      # TAKE OUT STUFF BELOW, rbind it, and go
      fitl = as.data.frame(do.call('rbind', fitl), stringsAsFactors=FALSE)

      #print(fitl)
      #snpres = cbind(snpres, fitl[match(1:nrow(snpres), fitl$snp_index), ])  ## debug, leave in snp_index to make sure merged in properly - check!
      snpres = cbind(snpres, fitl[match(1:nrow(snpres), fitl$snp_index), setdiff(names(fitl), "snp_index")])

      assign('GA_snpres', snpres, envir=.GlobalEnv)
      
      ## REFORMAT COLUMNS, just for the display...
      for(k in intersect(names(snpres), c('freq.1', 'beta_a', 'se_a', 'p_a', 'beta_d', 'se_d', 'p_d', 'beta_r', 'se_r', 'p_r')))
        snpres[[k]] = sprintf('%0.4g', snpres[[k]])
      
      output$analysis_summary = DT::renderDataTable(DT::datatable(snpres))

        # snpres = get("GA_snpres")
      pval = "p_a"
      if(!is.element(pval, names(snpres))){
        pval = "p_d"
        if(!is.element(pval, names(snpres))){
          pval = "p_r"
          if(!is.element(pval, names(snpres))){
            pval = names(snpres)[grep("^g[_]Pr[(].*[_]a$", names(snpres))]
          }
        }
      }

      # snpres = get("GA_snpres") ## for debugging the plot only
      ssnpres = snpres[!is.na(as.numeric(snpres[[pval]])), ]

      chr = ssnpres$chromosome
      pos = ssnpres$position
      pvalue = as.numeric(ssnpres[[pval]])

      getGroup = function(chr, pos, dist_break=1e6){
        group = rep(1, length(chr))
        for(ij in 2:length(pos)){
          #print(ij)
          group[ij] = group[ij-1]
          distprev = pos[ij] - pos[ij-1]
          samechr = (chr[ij] == chr[ij-1])
          #cat("distprev", distprev, "\n")
          #cat("samechr", samechr, "\n")
          #cat("distprev", (distprev >= dist_break), "\n")
          if(!samechr || (distprev >= dist_break)){
            group[ij] = group[ij-1] + 1
          }
        }
        #print(group)
        group
      }

      group = getGroup(chr, pos)

      posadj = pos
      ngroup = length(unique(group))
      if(ngroup > 1){
        max_within_diff = 0
        for(g in 1:ngroup){
          wc = which(group == g)
          lwc = length(wc)
          ccurpos = pos[wc]
          cur_diff = max(ccurpos[2:lwc] - ccurpos[1:(lwc-1)])
          max_within_diff = max(c(max_within_diff, cur_diff))
        }

        posadj = pos - min(pos)
        for(g in 2:ngroup){
          wc = which(group == g)
          wp = which(group == (g-1))
          posadj[wc] = posadj[wc] - min(posadj[wc]) + max(posadj[wp]) + 2*max_within_diff #1000
        }
      }

      #neglogp = -log(pvalue)
      #keep = which(!is.na(neglogp) & !is.infinite(neglogp))
      #print("pvalue")
      #print(pvalue)
      #print("-log(pvalue)")
      #print(-log(pvalue))

      #print(snpres)
      output$plot_man = renderPlotly(
        plot_ly(
          type="scatter",
          mode="markers",
          x=posadj, # xlab="Chromosome position",
          y=-log(pvalue, 10), # ylab=expression(-log[10](p)),
          text=sprintf("SNP: %s\nChr: %i\nPos: %i\nP: %s", ssnpres$snp.name, ssnpres$chromosome, ssnpres$position, pvalue),
          marker = list(
            size=10,
            #color = 'rgba(255, 182, 193, 0.9)',
            color=rainbow(ngroup)[group],
            #line = list(color = 'rgba(152, 0, 0, 0.8)', width=2)
            line = list(color="black", width=2)
          )
        ) %>%
          layout(
            xaxis=list(title="Chromosome position (relative, by group)", zerolinecolor='#ffff', zerolinewidth=2, gridcolor='fff'),
            yaxis=list(title=expression(-log[10](p)), zerolinecolor='#ffff', zerolinewidth=2, gridcolor='fff')
          )
      )
    },
    error = function(e){
      showNotification(paste("Error in analyzing the genotypes. Error:", as.character(e)), type='error', duration=DURATION, session=session)
    })
  })
  
  ## fs::path_home()
  #observeEvent(input$button_save, {
  #print("Button save...")
#  shinyFileChoose(input, 'phefile', root=getVolumes2(), session=session, filetypes=c('csv'))
  ## This could be an issues with R or with Shiny... unclear...
  observe({
    volumes = getVolumes2() # c("UserFolder"="c:/data")
    shinyFileSave(input, "button_save", roots=volumes, session=session, filetypes=c('csv'))
    fileinfo = parseSavePath(volumes, input$button_save)
    filename = fileinfo$datapath
    #print(paste("Filename:", filename))
    if(length(filename)>0 && !is.null(filename) && nchar(filename) > 0){
      if(file.exists(filename)){
        showNotification("Error in writing results - output file already exists, and will not be overwritten.", type='error', duration=DURATION, session=session)
      }else{
        write.csv(get("GA_snpres"), filename)
      }
    }
  })
  
  
  ## 2020-04-12: Added so that it actually terminates...  
  session$onSessionEnded(function() {
    stopApp()
  })
})


# Start the application
app = shinyApp(ui, server) # this won't run if sourced...

# 2020-04-12
#BROWSER = paste0(getwd(), '/ChromiumPortable/ChromiumPortable.exe')
BROWSER = ""
if(file.exists(BROWSER))
  options(browser=BROWSER)

runApp(app, launch.browser=TRUE)

# source('shinyGeneticsApp.R')
