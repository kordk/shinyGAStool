# Introduction

To compile this document:  

    pandoc --standalone --self-contained --toc github-create-demo.md -o github-create-demo.html
     
Global variables.
You can download plink and install it from
https://www.cog-genomics.org/plink2

    PLINK=plink_20200219


# Download summary statistics from the GWAS catalog, calculate a PRS, and simulate some data
    
    wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/HoffmannTJ_29507422_GCST007141/GERA-LDL.tsv.gz
    
```r
R --vanilla
source('scripts/plink-helper-v2.R')
sumstat = read.table(
  'GERA-LDL.tsv.gz', stringsAsFactors=FALSE, sep='\t', header=TRUE)
kgp = read.bim('kgp-eur.bim')
w = which(is.element(sumstat$SNP_ID, kgp$rs))
length(w)
sumstat = sumstat[w, ]
dim(sumstat)

ssumstat = sumstat[which(sumstat$P<=5e-8), ]
dim(ssumstat)

DIST = 1e6
keep = c()
while(nrow(ssumstat) > 0){
  w = which(ssumstat$P.value == min(ssumstat$P.value))[1]
  keep = c(keep, ssumstat$SNP_ID[w])
  ssumstat = ssumstat[
    which(ssumstat$chromosome!=ssumstat$chromosome[w] | 
          abs(ssumstat$position-ssumstat$position[w])>DIST), ]
}
length(keep)

# The input file should have one line per scored variant. By default, 
#  the variant ID is read from column 1, 
#  an allele code is read from the following column, 
#  and the score associated with the named allele is 
#   read from the column after the allele column;
osumstat = sumstat[
  is.element(sumstat$SNP_ID, keep), 
             c('SNP_ID', 'Allele.1', 'Estimate.Effect')]
dim(osumstat)
write.table(osumstat, 'kgp-eur-prs-score.txt', row.names=FALSE, 
            col.names=FALSE, quote=FALSE)
q()
```

    ${PLINK} \
      --bfile kgp-eur \
      --score kgp-eur-prs-score.txt \
      --out kgp-eur-prs-score-calc
      
    ${PLINK} \
      --bfile kgp-eur \
      --pca \
      --out kgp-eur-pca

```r
R --vanilla
score = read.table('kgp-eur-prs-score-calc.profile', 
                   header=TRUE, stringsAsFactors=FALSE)
pcs = read.table('kgp-eur-pca.eigenvec', header=!TRUE, stringsAsFactors=FALSE)
names(pcs) = c('FID', 'IID', paste0("PC", 1:20))
score$SCOREN = (score$SCORE - mean(score$SCORE))/sd(score$SCORE)

set.seed(13)
phe = score[, c('FID', 'IID')]
phe$age = rnorm(n=nrow(phe), mean=40, sd=10)
summary(phe)

phe$female = rbinom(n=nrow(phe), size=1, prob=0.55)
table(phe$female)

phe$ldl = rnorm(
  n = nrow(phe),
  mean = 70 + phe$age * 1 + phe$female * 1 + score$SCOREN * 5,
  sd = 1) # not sure about this last parameter
summary(phe)

cphe = cbind(phe, pcs[match(phe$IID, pcs$IID), paste0('PC', 1:6)])
write.csv(cphe, 'kgp-eur-ldl-pheno.csv', row.names=FALSE, quote=FALSE)
write.table(cphe, 'kgp-eur-ldl-pheno.tsv', 
            row.names=FALSE, col.names=TRUE, quote=FALSE)
q()
```

# Confirm our package matches plink results in this specific type of analysis.

    mkdir results
    ${PLINK} --bfile kgp-eur \
      --allow-no-sex \
      --linear \
      --pheno kgp-eur-ldl-pheno.tsv \
      --pheno-name ldl \
      --covar kgp-eur-ldl-pheno.tsv \
      --covar-name age female PC1 PC2 PC3 PC4 PC5 PC6 \
      --hide-covar \
      --out results/kgp-eur-plink-res
      
    gzip -f results/kgp-eur-plink-res.assoc.linear

```r
R --vanilla
score = read.table('kgp-eur-prs-score.txt', 
                   header=FALSE, stringsAsFactors=FALSE)
names(score) = c('rs', 'allele', 'beta')
res = read.table('results/kgp-eur-plink-res.assoc.linear.gz', 
                 stringsAsFactors=FALSE, header=TRUE)
sres = res[is.element(res$SNP, score$rs), ]
sres
q()
```


          CHR        SNP        BP A1 TEST NMISS     BETA    STAT         P
    3888     1 rs61775174  25787313  A  ADD   503 -0.55380 -1.7030 8.917e-02
    8358     1  rs6663252  55630151  C  ADD   503 -1.62900 -3.7390 2.063e-04
    10604    1  rs6678483  63074442  A  ADD   503 -0.90860 -2.6130 9.258e-03
    18039    1  rs1886683  92976023  T  ADD   503 -1.39200 -3.5240 4.639e-04
    21812    1   rs599839 109822166  G  ADD   503 -3.05100 -8.6080 1.003e-16
    41414    1  rs2642438 220970028  A  ADD   503 -0.09154 -0.2536 7.999e-01
    44963    1   rs508293 234848712  G  ADD   503 -0.44220 -1.3830 1.673e-01
    55276    2   rs541041  21294975  G  ADD   503 -3.14000 -8.4190 4.175e-16
    56366    2  rs1260326  27730940  T  ADD   503  0.79660  2.3330 2.003e-02
    61025    2  rs6544713  44073881  T  ADD   503  1.66300  4.6250 4.795e-06
    66507    2  rs2710642  63149557  G  ADD   503 -0.53910 -1.6360 1.025e-01
    76365    2 rs17512204 118732831  A  ADD   503 -0.46300 -0.7718 4.406e-01
    87174    2  rs2241342 169827918  T  ADD   503  0.66870  1.9390 5.308e-02
    119667   3 rs13317400  58408792  G  ADD   503 -0.78390 -1.3840 1.671e-01
    134006   3  rs2877601 122117566  C  ADD   503  0.10760  0.3200 7.491e-01
    135861   3 rs17412738 132257419  G  ADD   503 -0.54510 -1.1190 2.635e-01
    213076   5    rs12916  74656539  C  ADD   503  1.71200  5.4310 8.797e-08
    226286   5  rs2522062 131805416  G  ADD   503 -0.56590 -1.4130 1.584e-01
    231400   5  rs4704727 156380067  T  ADD   503 -0.83670 -2.3550 1.892e-02
    242669   6  rs6459450  16124560  C  ADD   503 -0.35360 -0.9789 3.281e-01
    246132   6  rs1800562  26093141  A  ADD   503 -1.48200 -1.8170 6.982e-02
    251476   6  rs6920323  31352060  C  ADD   503  1.11200  2.8190 5.011e-03
    253250   6 rs61670843  32574453  T  ADD   503  1.22100  2.2610 2.422e-02
    257790   6  rs2395943  42940673  A  ADD   503  0.03694  0.1099 9.126e-01
    277845   6  rs7776054 135418916  G  ADD   503 -0.51380 -1.3640 1.731e-01
    284938   6 rs10455872 161010118  G  ADD   503  2.42600  4.0360 6.299e-05
    296586   7 rs11764322  25952710  A  ADD   503  1.00700  2.9800 3.024e-03
    329228   8  rs4841133   9183664  A  ADD   503 -1.01900 -1.6560 9.836e-02
    340618   8 rs56204645  55421769  C  ADD   503  0.48420  1.1520 2.500e-01
    341612   8  rs6985620  59370159  T  ADD   503  0.41700  1.2320 2.185e-01
    353576   8  rs2737252 116663898  A  ADD   503 -0.45330 -1.1870 2.359e-01
    356164   8 rs17321515 126486409  G  ADD   503 -1.63700 -5.1480 3.801e-07
    384014   9 rs12686004 107653426  A  ADD   503 -0.36900 -0.7123 4.766e-01
    390574   9   rs651007 136153875  T  ADD   503  1.44600  3.7370 2.079e-04
    417430  10  rs2792751 113940329  T  ADD   503  1.15800  3.2950 1.055e-03
    427842  11 rs10832962  18656271  C  ADD   503 -0.92660 -2.5980 9.648e-03
    436001  11   rs174546  61569830  T  ADD   503 -1.97300 -6.0210 3.395e-09
    448187  11   rs964184 116648917  G  ADD   503  0.14800  0.3447 7.304e-01
    450648  11 rs77603146 126200090  A  ADD   503  1.46100  2.5840 1.006e-02
    480094  12  rs3184504 111884608  T  ADD   503 -0.89960 -2.7520 6.134e-03
    482517  12  rs7953249 121403724  G  ADD   503  0.40400  1.1800 2.386e-01
    512170  14 rs11620783  24871530  T  ADD   503  0.38260  1.1800 2.386e-01
    523830  14 rs11158869  71057558  C  ADD   503 -0.54120 -1.1740 2.410e-01
    529742  14    rs17580  94847262  A  ADD   503  1.06900  1.5070 1.326e-01
    562136  16 rs12149545  56993161  A  ADD   503 -0.51150 -1.4140 1.579e-01
    564918  16  rs2000999  72108093  A  ADD   503  0.75980  1.8690 6.224e-02
    570143  16 rs67890964  83979317  C  ADD   503 -0.49900 -1.4730 1.415e-01
    572709  17 rs55714927   7080316  T  ADD   503 -0.61440 -1.4220 1.558e-01
    578908  17  rs7210738  45757598  G  ADD   503  0.18090  0.5454 5.858e-01
    581300  17  rs2645492  57875554  A  ADD   503 -0.80740 -1.9310 5.402e-02
    582061  17  rs1801689  64210580  C  ADD   503  3.41300  4.5080 8.178e-06
    582806  17 rs77542162  67081278  G  ADD   503  2.15200  1.3530 1.767e-01
    607305  19 rs55997232  11188117  T  ADD   503 -3.80400 -8.1530 2.961e-15
    608356  19 rs10500212  19723215  T  ADD   503 -1.13700 -1.9070 5.708e-02
    611855  19  rs4420638  45422946  G  ADD   503  2.98100  7.6450 1.098e-13
    612387  19   rs492602  49206417  G  ADD   503  0.76780  2.3730 1.802e-02
    614833  19 rs34503352  58651296  A  ADD   503 -0.15150 -0.3458 7.296e-01
    620897  20  rs2252247  17842006  A  ADD   503  0.59910  1.2360 2.170e-01
    624042  20  rs2207132  39142516  A  ADD   503  2.03100  2.4380 1.512e-02

# Now run the shiny interface, and compare...

    R --vanilla
    source('shinyGeneticsApp.R')
    shinyApp(ui, server)


