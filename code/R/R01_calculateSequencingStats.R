## Alice Balard
## November 2021 (updated Jan23)

# Each script sources the previous script of the pipeline if needed
source("R00_rawDataCleaning.R")

#########################
## Data loading & prep ##
#########################
# Raw data including Kostas previous results to compare and have the data
rawData <- read.csv("../../dataIn/cleanedRawData132fishG1G2.csv")
rawData$SampleID <- rawData$ID

## Raw reads quality check
file <- read.csv("../../dataIn/multiqc_general_stats_rawReads.txt") 
file$SampleID = gsub("_L00.*","", gsub(".*-L1_","", file$Sample.Name))
names(file) <- c("Sample.Name", "percent_duplicates_rawReads", "percent_GC_rawReads", "M.Seqs_rawReads", "SampleID")

## Trimmed reads quality check
fileT <- read.table("../../dataIn/multiqc_general_stats_trimmedReadsCutadapt.txt", header = T) 
names(fileT) <- gsub("FastQC_mqc.generalstats.fastqc.", "", names(fileT))
fileT$SampleID = gsub("_L00.*","", gsub(".*-L1_","", fileT$Sample))
names(fileT) <- c("Sample", "percent_duplicates_trimmedReads", "percent_gc_trimmedReads", "avg_sequence_length_trimmedReads", "percent_fails_trimmedReads", "total_sequences_trimmedReads", "SampleID")   

## Mapping efficiency after BSBolt
mappability <- read.table(file = "../../dataIn/Report_mapping_efficiency.txt", header = T)
mappability2 <- read.table(file = "../../dataIn/Report_mapping_efficiencySexUn.txt", header = T)

mappability <- merge(mappability, mappability2) %>% 
  mutate(SampleID = gsub("_L00*.","", gsub(".*-L1_","", gsub("_R1_001_trimmed_cutadapt.fastq.gz.bam.sorted.bam", "", filename))),
         mappabilityTotal = Mapped_read/Total_reads*100, 
         mappabilityNoSexUn = Mapped_reads_minusSexAndUnknownChr/Total_reads_minusSexAndUnknownChr*100,
         mappabilitySexChr = Mapped_read_sexchr/Total_reads*100,
         mappabilityUnChr = Mapped_reads_unknownchr/Total_reads*100)

## Nbr of methylated sites:
methylBSdf <- read.delim("../../dataIn/BSBolt_Methcall_Methstats_nbrMethylatedCpGperSample.txt", header = FALSE)
names(methylBSdf) <- c("SampleID", "NbrMethylatedCpG")
methylBSdf$SampleID <- gsub("_L00*.","", gsub(".*-L1_","", gsub("_R1_001", "", methylBSdf$SampleID)))

## Merge metadata:
fullMetadata <- merge(merge(merge(merge(methylBSdf, mappability[-1]), file[-1]), fileT[-1]), rawData)

#######################
# Mapping stats:
print("Number of raw reads:")
print(mean(fullMetadata$M.Seqs_rawReads)) # 11.02M reads
print("+/-")
print(qnorm(0.975)*sd(fullMetadata$M.Seqs_rawReads)/sqrt(nrow(fullMetadata))) # +/-0.38

print("Mapping efficiency with BSBolt in %:")
print(mean(fullMetadata$mappabilityTotal))  # 85.6
print("+/-")
print(qnorm(0.975)*sd(fullMetadata$mappabilityTotal)/sqrt(nrow(fullMetadata))) # +/- 0.51

## After exploration of raw data (multiQC files), we decide to remove 4 samples
# with bad quality and/or less than 6M reads after trimming ("S12", "S22", "S110", "S118", "S142")
fullMetadata <- fullMetadata[!fullMetadata$ID %in% c("S12", "S22", "S110", "S118", "S142"),]

print("Number of individuals considered:")
print(sum(table(fullMetadata$Generat)))
print("including parents and offspring:")
print(table(fullMetadata$Generat))
# N=127: 111 offspring + 16 parents

##############################
## Prepare sub metadatasets ##
##############################

## Change one term: NbrMethylatedCpG is the one calculated by BSBolt IN TOTAL
names(fullMetadata)[names(fullMetadata) %in% "NbrMethylatedCpG"] <- "NbrMethylatedCpG_global_BSBolt"

# give a numerical value to treatment, for Methylkit
fullMetadata$trtG1G2_NUM <- as.numeric(as.factor(fullMetadata$trtG1G2))

## relevel treatments for graphs
fullMetadata$trtG1G2 <- factor(as.factor(fullMetadata$trtG1G2), levels = c("Control", "Exposed","NE_control", "NE_exposed", "E_control", "E_exposed"  ))

## family as factor for models
fullMetadata$Family <- as.factor(fullMetadata$Family)

# paternal exposure
fullMetadata$PAT="Exposed father group"
fullMetadata$PAT[fullMetadata$trtG1G2 %in% c("Control", "NE_control", "NE_exposed")]="Control father group"

## Add brother pairs
fullMetadata$brotherPairID <- sapply(str_split(fullMetadata$clutch.ID, "_"), "[", 2 )
## avoid confusion with numeric
fullMetadata$brotherPairID <- paste0("BP",fullMetadata$brotherPairID)

########################
## Parents only metadata
fullMetadata_PAR <- fullMetadata[fullMetadata$Generat %in% "P",]

##########################
## Offspring only metadata
fullMetadata_OFFS <- fullMetadata[fullMetadata$Generat %in% "O",]
fullMetadata_OFFS$trtG1G2 <- droplevels(fullMetadata_OFFS$trtG1G2)

## Create variable for offsping and parents separated
fullMetadata_OFFS$offsTrt <- "controlO"
fullMetadata_OFFS$offsTrt[fullMetadata_OFFS$Tr %in% c("TT", "CT")] <- "infectedO"
fullMetadata_OFFS$patTrt <- "controlP"
fullMetadata_OFFS$patTrt[fullMetadata_OFFS$Tr %in% c("TC", "TT")] <- "infectedP"

## Sanity check
table(fullMetadata_OFFS$offsTrt, fullMetadata_OFFS$trtG1G2)
table(fullMetadata_OFFS$patTrt, fullMetadata_OFFS$trtG1G2)

## REORDER metadata by sample ID
fullMetadata = fullMetadata[order(as.numeric(gsub("S", "", fullMetadata$SampleID))),]
fullMetadata_PAR = fullMetadata_PAR[order(as.numeric(gsub("S", "", fullMetadata_PAR$SampleID))),]
fullMetadata_OFFS = fullMetadata_OFFS[order(as.numeric(gsub("S", "", fullMetadata_OFFS$SampleID))),]

## Export summary table:
write.csv(fullMetadata, "../../dataIn/fullMetadata127_Alice.csv", row.names=FALSE, quote=FALSE)

## Does mappability in offspring differs by sex (with all chr/without sex & unknown ones)? No.
ggplot(fullMetadata_OFFS, aes(x=Sex, y=mappabilitySexChr)) + geom_violin()
anova(lm(mappabilitySexChr~Sex, data = fullMetadata_OFFS))
# Sex         1 23.1742 23.1742  424.54 < 2.2e-16 ***

fullMetadata_OFFS %>% group_by(Sex) %>% dplyr::summarise(mean=mean(mappabilitySexChr),
                                                         n=n(),
                                                         margin_error = qt(p= 0.05/2, df=n - 1, lower.tail=F) *
                                                           sd(mappabilitySexChr) / sqrt(n()),
                                                         lowCI=mean-margin_error, highCI=mean+margin_error)
# Females: mean=3.74%, 95% CI[3.69-3.79], Males: mean=2.82%, 95% CI[2.75-2.90]

anova(lm(mappabilityUnChr~Sex, data = fullMetadata_OFFS))
anova(lm(mappabilityTotal~Sex, data = fullMetadata_OFFS))
anova(lm(mappabilityNoSexUn~Sex, data = fullMetadata_OFFS))

fullMetadata_OFFS %>% group_by(Sex) %>% dplyr::summarise(mean=mean(mappabilityNoSexUn),
                                                         n=n(),
                                                         margin_error = qt(p= 0.05/2, df=n - 1, lower.tail=F) *
                                                           sd(mappabilityNoSexUn) / sqrt(n()),
                                                         lowCI=mean-margin_error, highCI=mean+margin_error)
# Females: mean=84.4%, 95% CI[83.7-85.1], Males: mean=83.5%, 95% CI[82.7-84.4]

## clean workspace
rm(file, fileT, mappability, methylBSdf, rawData)
