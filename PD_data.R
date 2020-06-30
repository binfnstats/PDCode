# Required Packages

library('GEOquery')
library('Biobase')
library('limma')
library('tidyverse')

##########

# Parametres for analysis

## Top k d.e. probes
k = 100

##########

# Load Data

## load series and platform data from GEO
gset = GEOquery::getGEO("GSE99039", GSEMatrix =TRUE, AnnotGPL=TRUE)

if (length(gset) > 1) idx = grep("GPL570", attr(gset, "names")) else idx = 1
gset = gset[[idx]]

## Make coloumn names match toptable
fvarLabels(gset) = make.names(fvarLabels(gset))

## group names for all samples
gsms = base::paste0(
  "011011111X1010010XX00001X000X00001X111XXX001101111",
  "XX0X01X1010XX1111111110XXXX1X1110110XXXXXXXXX01010",
  "XXX01000010010000X00000000001X0X1111X1001011111000",
  "00010011001XX10000X000010XXXXX11111X11X11111XXXX11",
  "0110XX10111110110111110111011X11011011X110X0111011",
  "1100011000X101000000001000X11000000000001011110000",
  "0010XX00110X0100000010000000XXX110XX1111X100001011",
  "10000X1101011000100111000010000X0X001X1X1011001111",
  "0010001010111X00X000100011101111101X00000X1XXXXXXX",
  "XXXXXXXXXXXX10XX11X0XX000X0000000000111101XX011X11",
  "111111111111X10000000000X00XX0X011X1X0110X00XX000X",
  "XX0X0XXX"
)
sml = c()
for (i in 1:nchar(gsms)) {
  sml[i] = base::substr(gsms, i, i)
}

## Eliminate samples marked as "X"
sel = which(sml != "X")
sml = sml[sel]
gset = gset[, sel]

## log2 transform
ex = Biobase::exprs(gset)
qx = as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
LogC <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
}

## set up the data and proceed with analysis
sml = paste("G", sml, sep = "")
fl = as.factor(sml)
gset$description = fl
design = stats::model.matrix( ~ description + 0, gset)
colnames(design) = levels(fl)
fit = limma::lmFit(gset, design)
cont.matrix = limma::makeContrasts(G1 - G0, levels = design)
fit2 = limma::contrasts.fit(fit, cont.matrix)
fit2 = limma::eBayes(fit2, 0.01)
tT = limma::topTable(fit2,
               adjust = "fdr",
               sort.by = "B",
               number = k)

tT =
  subset(
    tT,
    select = c(
      "ID",
      "adj.P.Val",
      "P.Value",
      "t",
      "B",
      "logFC",
      "Gene.symbol",
      "Gene.title"
    )
  )

names = Biobase::featureNames(gset)
tT$ID = factor(tT$ID, levels = names)
gene_new = c(tT$ID)
gset_new = gset[gene_new]

data = gset_new %>% 
  Biobase::exprs(.) %>% 
  t(.)

y = design[,1] %>% 
  as.numeric(.)

