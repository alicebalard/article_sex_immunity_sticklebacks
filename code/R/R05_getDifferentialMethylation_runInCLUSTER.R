## Differential methylation analyses
## A. Balard
## March 2023

## Each script sources the previous script of the pipeline if needed
source("R04_GlobalMethylation.R")

#######################################
### Differential Methylation Sites ####
#######################################

####################################################################
## Calculate DMS between sexes, in each treatment, per brother pair:
getDMSsexperBP <- function(BP, uniteCov=uniteCovHALFperSex_G2_woSexAndUnknowChr,metadata = fullMetadata_OFFS){

    ## subset for the BP considered
    metadata = metadata[metadata$brotherPairID %in% BP,]

    ## Update male/female codes 
    ## Males will be coded as 0, females as 1    
    metadata$sex_NUM <- 0
    metadata$sex_NUM[metadata$Sex %in% "F"] = 1

    ## remove individuals in treatment/BP where only one sex is present
    tab = table(metadata$Sex, metadata$trtG1G2)  %>% data.frame %>% dplyr::filter(Freq ==0)
    
    ## rm treatment without both males & females in this BP
    metadata = metadata[!metadata$trtG1G2 %in% tab$Var2,]

    ## Calculate DMS between sexes, in each treatment, per brother pair:
    
    ## Unite object for one treatment:
    metadataBP_CC = metadata[metadata$trtG1G2 %in% c("NE_control"), ]
    metadataBP_CT = metadata[metadata$trtG1G2 %in% c("NE_exposed"), ]
    metadataBP_TC = metadata[metadata$trtG1G2 %in% c("E_control"), ]
    metadataBP_TT = metadata[metadata$trtG1G2 %in% c("E_exposed"), ]
    
    ## Make 4 separate uniteCov:
    myuniteCovBP_CC = methylKit::reorganize(methylObj = uniteCov,
                                            treatment = metadataBP_CC$sex_NUM,
                                            sample.ids = metadataBP_CC$ID)
    myuniteCovBP_CT = methylKit::reorganize(methylObj = uniteCov,
                                            treatment = metadataBP_CT$sex_NUM,
                                            sample.ids = metadataBP_CT$ID)
    myuniteCovBP_TC = methylKit::reorganize(methylObj = uniteCov,
                                            treatment = metadataBP_TC$sex_NUM,
                                            sample.ids = metadataBP_TC$ID)
    myuniteCovBP_TT = methylKit::reorganize(methylObj = uniteCov,
                                            treatment = metadataBP_TT$sex_NUM,
                                            sample.ids = metadataBP_TT$ID)

    ## remove bases where NO fish in this BP has a coverage
    myuniteCovBP_CC = methylKit::select(myuniteCovBP_CC, which(!is.na(rowSums(percMethylation(myuniteCovBP_CC)))))
    myuniteCovBP_CT = methylKit::select(myuniteCovBP_CT, which(!is.na(rowSums(percMethylation(myuniteCovBP_CT)))))
    myuniteCovBP_TC = methylKit::select(myuniteCovBP_TC, which(!is.na(rowSums(percMethylation(myuniteCovBP_TC)))))
    myuniteCovBP_TT = methylKit::select(myuniteCovBP_TT, which(!is.na(rowSums(percMethylation(myuniteCovBP_TT)))))

    ## Calculate differential methylation:
    ## We select the bases that have q-value<0.01 and percent methylation difference larger than 15%, sex as covariate

    ## Prepare in case there are no individuals present to calculate DMS:
    DMS_15pc_BP_CC="both sex are not represented in this treatment/family group"
    DMS_15pc_BP_CT="both sex are not represented in this treatment/family group"
    DMS_15pc_BP_TC="both sex are not represented in this treatment/family group"
    DMS_15pc_BP_TT="both sex are not represented in this treatment/family group"

    if (length(myuniteCovBP_CC@sample.ids)>0){
        DMS_15pc_BP_CC = getDiffMethSex(myuniteCov = myuniteCovBP_CC, myMetadata = metadataBP_CC, mccores = 1, mydif = 15)
    }

    if (length(myuniteCovBP_CT@sample.ids)>0){
        DMS_15pc_BP_CT = getDiffMethSex(myuniteCov = myuniteCovBP_CT, myMetadata = metadataBP_CT, mccores = 10, mydif = 15)
    }

    if (length(myuniteCovBP_TC@sample.ids)>0){
        DMS_15pc_BP_TC = getDiffMethSex(myuniteCov = myuniteCovBP_TC, myMetadata = metadataBP_TC, mccores = 10, mydif = 15)
    }

    if (length(myuniteCovBP_TT@sample.ids)>0){
        DMS_15pc_BP_TT = getDiffMethSex(myuniteCov = myuniteCovBP_TT, myMetadata = metadataBP_TT, mccores = 10, mydif = 15)
    }

    return(list(DMSlist = list(DMS_15pc_BP_CC = DMS_15pc_BP_CC, DMS_15pc_BP_CT = DMS_15pc_BP_CT, DMS_15pc_BP_TC = DMS_15pc_BP_TC, DMS_15pc_BP_TT = DMS_15pc_BP_TT)))
}

run = TRUE

## Calculate DMS by brother pairs:

if (run == TRUE){

    ## We will apply the following function to all BP:
    vecBP <- unique(fullMetadata_OFFS$brotherPairID)

    ## Loop over all BP
    DMlist <- list() # empty plot list
    for (i in 1:length(vecBP)){
        DMlist[[i]] <- getDMSsexperBP(BP = vecBP[[i]])
    } 
    names(DMlist) <- vecBP

    saveRDS(DMlist, "../../dataOut/DiffMeth/DMperBP_SexwithinTrt_Mar2022.RDS")
}
