## MethylKit object preparation
## A. Balard
## 3rd of March 2023

# Each script sources the previous script of the pipeline if needed
# NB: comment when on Apocrita
source("R02_prepBSBOLTForMethylkit_runInCLUSTER.R")

## Load unitecov objects
base::load("../../gitignore/uniteCovALL_woSexAndUnknowChr_March2023.RData")
base::load("../../gitignore/uniteCovALL_G2_woSexAndUnknowChr_March2023.RData")
base::load("../../gitignore/uniteCovHALFperSex_G2_woSexAndUnknowChr_March2023.RData")

##################
## To recalculate:
rerun = FALSE
if (rerun == TRUE){
  ##### Load prepared dataset (in APOCRITA) #####
  dataPath="/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/04BSBolt_methCall/BSBolt/MethylationCalling/Methylation_calling_splitted/formatCG4methylKit"
  
  temp = list.files(path=dataPath,
                    pattern = ".CG4methylkit.txt",
                    full.names = T)
  
  ## Add metadata on sex (the grouping factor)
    metadata <- readxl::read_xlsx("../../dataIn/raw data Joshka Kostas/Kostas_G2_info.xlsx")
    
    ## Males will be coded as 0, females as 1    
    metadata$sex_NUM <- 0
    metadata$sex_NUM[metadata$Sex %in% "F"] = 1
    
  ### Make methylkit object
  myobj=methylKit::methRead(as.list(temp),
                            mincov=10,
                            sample.id=as.list(metadata$ID),
                            assembly="Gynogen_pchrom_assembly_all",
                            treatment=metadata$sex_NUM,
                            context="CpG")
  
  #################################################
  ## We remove several samples from the raw dataset
  fullMetadata <- read.csv("../../dataIn/fullMetadata127_Alice.csv") 
  
  ## create a new methylRawList object
  print("Remove unwanted samples")
  
  myobj=reorganize(
    myobj,
    sample.ids=metadata$ID[metadata$ID %in% fullMetadata$ID],
    treatment=metadata$sex_NUM[metadata$ID %in% fullMetadata$ID])
  
  ###############################
  ## Filtering and normalising ##
  ###############################
  
  ## Filtering based on coverage:
  # It might be useful to filter samples based on coverage. Particularly, if our samples are suffering from PCR bias it would be useful to discard bases with very high read coverage. Furthermore, we would also like to discard bases that have low read coverage, a high enough read coverage will increase the power of the statistical tests. The code below filters a methylRawList and discards bases that have coverage below 10X and also discards the bases that have more than 99.9th percentile of coverage in each sample.
  print("Filter")
  filtered.myobj=filterByCoverage(myobj, lo.count=10, lo.perc=NULL,
                                  hi.count=NULL, hi.perc=99.9)
  
  ## normalise the coverage
  print("Normalise")
  normFil.myobj=normalizeCoverage(filtered.myobj)
  
  #####################
  ## MERGING SAMPLES ##
  #####################
  ##In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. The unite() function will return a methylBase object which will be our main object for all comparative analysis. The methylBase object contains methylation information for regions/bases that are covered in all samples.
  
  print("Add CpG present in ALL individuals")
  uniteCovALL= methylKit::unite(normFil.myobj, mc.cores=8)
  uniteCovALL=as(uniteCovALL,"methylBase")
    
  ## In all OFFSPRING:
  uniteCovALL_G2 = reorganize(
    normFil.myobj,
    sample.ids=metadata$ID[metadata$ID %in% fullMetadata$ID & metadata$Generat %in% "O"],
    treatment=metadata$sex_NUM[metadata$ID %in% fullMetadata$ID & metadata$Generat %in% "O"])
  
  uniteCovALL_G2 = methylKit::unite(uniteCovALL_G2, mc.cores=8)# try with 8 cores
  uniteCovALL_G2 = as(uniteCovALL_G2,"methylBase")
  
## Keep methylated CpG sites observed in at least 50% individual fish per sex (N=27) after filtering and normalising:
table(fullMetadata[fullMetadata$Generat %in% "O","Sex"])

## OFFSPRING
  uniteCovHALFperSex_G2 = reorganize(
    normFil.myobj,
    sample.ids=metadata$ID[metadata$ID %in% fullMetadata$ID & metadata$Generat %in% "O"],
    treatment=metadata$sex_NUM[metadata$ID %in% fullMetadata$ID & metadata$Generat %in% "O"])
  
  uniteCovHALFperSex_G2 = methylKit::unite(uniteCovHALFperSex_G2, min.per.group=27L, mc.cores=8)# try with 8 cores
  uniteCovHALFperSex_G2 = as(uniteCovHALFperSex_G2,"methylBase")
  
  ########################################################################################
  ## Remove reads from sex chromosome X ("Gy_chrXIX") and unmapped contigs ("Gy_chrUn") ##
  ########################################################################################
  print("nbr CpG on sex chromosome of unmapped:")
  nrow(uniteCovALL[uniteCovALL$chr %in% c("Gy_chrXIX", "Gy_chrUn"),])
  
  print("Keep CpG apart from sex chromosome XIX and unmapped (comprise Y chr)")
  uniteCovALL_woSexAndUnknowChr=uniteCovALL[!uniteCovALL$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]
  uniteCovALL_G2_woSexAndUnknowChr=uniteCovALL_G2[!uniteCovALL_G2$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]
  uniteCovHALFperSex_G2_woSexAndUnknowChr=uniteCovHALFperSex_G2[!uniteCovHALFperSex_G2$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]
    
  print("nbr CpG shared by all 127 samples:")
  length(uniteCovALL_woSexAndUnknowChr$chr)
    print("nbr CpG shared by all offsprings:")
  length(uniteCovALL_G2_woSexAndUnknowChr$chr)
  print("nbr CpG shared by 50% (>27) of offsprings per sex:")
  length(uniteCovHALFperSex_G2_woSexAndUnknowChr$chr)
  
  ######################################################################################################
  ## Save outcomes
  save(uniteCovALL_woSexAndUnknowChr,
       file = "/data/SBCS-EizaguirreLab/Alice/SexImmunity/article_sex_immunity_sticklebacks/gitignore/uniteCovALL_woSexAndUnknowChr_March2023.RData")
  save(uniteCovALL_G2_woSexAndUnknowChr,
       file = "/data/SBCS-EizaguirreLab/Alice/SexImmunity/article_sex_immunity_sticklebacks/gitignore/uniteCovALL_G2_woSexAndUnknowChr_March2023.RData")
  save(uniteCovHALFperSex_G2_woSexAndUnknowChr,
       file = "/data/SBCS-EizaguirreLab/Alice/SexImmunity/article_sex_immunity_sticklebacks/gitignore/uniteCovHALFperSex_G2_woSexAndUnknowChr_March2023.RData")
}
