## load all libraries needed for the project

list.of.packages <- c(
  "ape", #for reag.gff
  "cAIC4",
  "ComplexUpset", # for prettier upset plots
  "cowplot",
  "dendextend", # The package dendextend contains many functions for changing the appearance of a dendrogram and for comparing dendrograms.
  "devtools",
  "dplyr",
  "emmeans", ## for post-hoc Tukey tests
  "factoextra", # color PCA plots
  "FactoMineR", # for PCA
  "forcats", # keeps characters in previous order axis ggplot (for bubble plot)
  "ggeffects", # to plot random effects predictions
  "ggplot2",
  "ggpubr", ## to merge ggplot2 plots
  "ggrepel",
  "goeveg", # find the best number of dimensions for NMDS
  "ggsignif", ## for significance bars on ggplot
  "ggVennDiagram",## Venn diagram in ggplot
  "grid",
  "gridExtra",
  "lme4", ## for mixed models
  "lmtest", # for lrtests
  "lmerTest", # for stepwise analysis of lmer
  "magrittr",      # provides the %>% operator
  "MuMIn", # participation of variables to the variance
  "missMDA",# PCA for incomplete data
  "nlme", ## for mixed models
  "pheatmap", # for heatmaps
  "plyr", # for join (keep row order",
  "png",
  "purrr",
  "qualpalr",# extra palettes
  "RColorBrewer", # for colors in Venn diagrams
  "rentrez", # to extract info from NCBI Entrez
  "reshape2",
  "sjPlot", # plot interaction effects
  "slider", # for slidding windows
  "splitstackshape", # to spread the V9 column of gff into columns by key
  "stringr", # to modify characters
  "tidyverse",  # tidyverse will pull in ggplot2, readr, other useful libraries
  "UpSetR", # for upset plots
  "VCA",
  "vegan", ## for Adonis
  "VennDiagram")

###################################################################
## install from CRAN and require all libraries from CRAN and github
ipak <- function(pkg){
#  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
# if (length(new.pkg))
#    install.packages(new.pkg, dependencies = TRUE,repos = "http://cran.us.r-project.org")
  sapply(pkg, require, character.only = TRUE)
}

ipak(list.of.packages)

##########################################
## install packages from github if not yet

#if (length("rentrez"[!("rentrez" %in% installed.packages()[, "Package"])])){
#    install_github("ropensci/rentrez")
#}

#if (length("goEnrichment"[!("goEnrichment" %in% installed.packages()[, "Package"])])){
#    install_github("asishallab/goEnrichment")   
#}

#if (length("/pairwiseAdonis"[!("pairwiseAdonis" %in% installed.packages()[, "Package"])])){
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#}

#if (length("ggVennDiagram"[!("ggVennDiagram" %in% installed.packages()[, "Package"])])){
#install_github("gaospecial/ggVennDiagram")
#}


#####################################################
## install from biocmanager and require all libraries
## Biocmanager packages 
list.bioc <- c("AnnotationDbi", # gene annotation from online databases
               "biomaRt", # to retrieve genes descriptions
               "Category", # for hypergeometric GO test
               "genomation", ## for annotation
               "GenomicFeatures",## for annotation
               "goEnrichment",
               "GOstats", # for GO analysis
               "GSEABase",  # for GO term GeneSetCollection
               "methylKit",
               "org.Hs.eg.db" # gene annotation from online databases
               ) 
ipak2 <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#  if (length(new.pkg))
#   BiocManager::install(new.pkg)
  sapply(pkg, require, character.only = TRUE)
}

ipak2(list.bioc)

## offspring colors for all kind of plots
colOffs <- c("#ffe67f", "#ff6300","#a8caff","#a800d4")

theme_set(theme_pubr())
