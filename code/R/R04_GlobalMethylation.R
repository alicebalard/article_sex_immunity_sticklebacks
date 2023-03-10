## MethylKit object preparation
## A. Balard
## January 2023

# Each script sources the previous script of the pipeline if needed
# NB: comment when on Apocrita
source("R03_prepObjectMethylkit_runInCLUSTER.R")

## Calculate BCI
# Body condition of the G2 fish estimated by residuals from the regression of body mass on body length (Chellappaet al.1995).
fullMetadata_OFFS$BCI2 <- residuals(lmer(Wnettofin ~ Slfin + (1|brotherPairID), data=fullMetadata_OFFS))

# Does sex affect BCI? NOPE
mod = lmerTest::lmer(BCI2 ~ Sex * outcome * PAT + (1|brotherPairID), data=fullMetadata_OFFS)
step(mod)
plot(ggpredict(mod, terms = c("outcome", "Sex", "PAT")), add.data=T)

# Does sex affect tolerance? NOPE
mod = lmerTest::lmer(BCI2 ~ No.Worms * Sex * PAT + (1|brotherPairID), data=fullMetadata_OFFS)
step(mod)
# No.Worms:PAT              0  1   21085.1 394415 913.49  6.0432 0.01556 *
P1 <- plot(ggpredict(mod, terms = c("No.Worms", "Sex", "PAT")), add.data=T)+
  ylab("Body Condition Index")+xlab("Number of worms")+scale_x_continuous(breaks = 0:10)+
  ggtitle("Predicted values of Body Condition Index")

# Does sex affect resistance? NOPE
mod = lmerTest::lmer(No.Worms ~ Sex * PAT + (1|brotherPairID), 
                     data=fullMetadata_OFFS[fullMetadata_OFFS$outcome %in% "infected",])
step(mod)
P2 <- plot(ggpredict(mod, terms = c("PAT", "Sex")), add.data=T)+
  ylab("Number of worms")+xlab(NULL)+scale_y_continuous(breaks = 0:10)+
  ggtitle("Predicted values of number of worms")

pdf("../../dataOut/figures/figA_norestolSexeffect.pdf", width = 10, height = 4)
ggarrange(P1, P2, ncol=2, common.legend = TRUE, legend="right")
dev.off()

## Sagonas 2020: "or findings showed significant interactions between treatment and 
# 1) respiratory burst activity (measure of innate immune response)
# 2) liver weight
# 3) head kidney weight 
# on RMS. 
# However, body condition showed no significant association with RMS

# Correction for fish length
fullMetadata_OFFS$integration.RLU.resSlfin[!is.na(fullMetadata_OFFS$integration.RLU.)] = 
  residuals(lm(integration.RLU.~Slfin, data=fullMetadata_OFFS))

fullMetadata_OFFS$peak.RLU.sec.resSlfin[!is.na(fullMetadata_OFFS$peak.RLU.sec.)] = 
  residuals(lm(peak.RLU.sec.~Slfin, data=fullMetadata_OFFS))

fullMetadata_OFFS$Tpeak.min.resSlfin[!is.na(fullMetadata_OFFS$Tpeak.min.)] = 
  residuals(lm(Tpeak.min.~Slfin, data=fullMetadata_OFFS))

fullMetadata_OFFS$L.W.resSlfin[!is.na(fullMetadata_OFFS$L.W)] = 
  residuals(lm(L.W~Slfin, data=fullMetadata_OFFS))

fullMetadata_OFFS$HK.W.resSlfin[!is.na(fullMetadata_OFFS$HK.W)] = 
  residuals(lm(HK.W~Slfin, data=fullMetadata_OFFS))

# 1) Sex (but not trt or father trt) affects respiratory burst activity: 
## lower integration[RLU] for males vs females (p=0.004)
## lower peak[RLU/sec] for males vs females (p=0.002)
## no difference in Tpeak[min]

mod = lmerTest::lmer(integration.RLU.resSlfin ~ Sex * PAT * outcome + (1|brotherPairID), data=fullMetadata_OFFS)
step(mod) # sex signif p=0.004
plot(ggpredict(mod, terms = c("outcome", "Sex", "PAT")), add.data=T)

mod = lmerTest::lmer(peak.RLU.sec.resSlfin ~ Sex * PAT * outcome + (1|brotherPairID), data=fullMetadata_OFFS)
step(mod) # sex signif p=0.002
plot(ggpredict(mod, terms = c("outcome", "Sex", "PAT")), add.data=T)

mod = lmerTest::lmer(Tpeak.min.resSlfin ~ Sex * PAT * outcome + (1|brotherPairID), data=fullMetadata_OFFS)
step(mod) # nothing signif
plot(ggpredict(mod, terms = c("outcome", "Sex", "PAT")), add.data=T)

# 2) liver weight: lower liver weight for male vs female (p<0.001), infected vs control (p<0.001), control father vs infected father (p<0.001)
mod = lmerTest::lmer(L.W.resSlfin ~ Sex * PAT * outcome + (1|brotherPairID), data=fullMetadata_OFFS)
step(mod) # sex (p<0.001), paternal trt (p<0.001) and trt (p<0.001) significant [additive]
plot(ggpredict(mod, terms = c("outcome", "Sex", "PAT")), add.data=T)

# 3) head kidney weight: lower head kidney weight upon infection (p=0.025)
mod = lmerTest::lmer(HK.W.resSlfin ~ Sex * PAT * outcome + (1|brotherPairID), data=fullMetadata_OFFS)
step(mod) # trt (p=0.025) significant 
plot(ggpredict(mod, terms = c("outcome", "Sex", "PAT")), add.data=T)

##################
print("Number of CpGs present in all G2 samples:")
print(nrow(uniteCovALL_G2_woSexAndUnknowChr))
           

# Calculate number of methylated sites, mean coverage, and residuals of methylated sites by covered sites (to account for coverage bias)
mycalcRMS <- function(myUniteCov, myMetaData){
  percMethMat = methylKit::percMethylation(myUniteCov)
  # create a dataframe with all info
  percMethDF = data.frame(SampleID = colnames(percMethMat),
                          Nbr_hypermethCpG = colSums(percMethMat>=70 & !is.na(percMethMat)), ## number of methylated sites
                          Nbr_hypomethCpG = colSums(percMethMat<=30 & !is.na(percMethMat)), ## number of methylated sites
                          Nbr_coveredCpG = colSums(!is.na(percMethMat)), ## number of sites covered in this sample
                          Nbr_NOTcoveredCpG = colSums(is.na(percMethMat)),## number of sites NOT covered in this sample
                          MeanCoverage = colMeans(methylKit::getData(myUniteCov)[,myUniteCov@coverage.index], na.rm = T), ## coverage.index: vector denoting which columns in the data correspond to coverage values
                          OverallPercentageMethylation = colMeans(methylKit::percMethylation(myUniteCov), na.rm = T))
  
  ## RMS in this sample based on covered sites
  percMethDF$RMS_coveredCpG = percMethDF$Nbr_hypermethCpG / percMethDF$Nbr_coveredCpG
  ## merge with original metadata:
  myMetaData = merge(myMetaData, percMethDF)
  # calculate also RMS global, considering CpG covered or not (to compare)
  myMetaData$RMS_allCpG_coveredOrNot = myMetaData$Nbr_hypermethCpG / (myMetaData$M.Seqs_rawReads*10e6)
  # calculate residuals of nbr of methCpG by nbr of covered CpG
  myMetaData$res_Nbr_methCpG_Nbr_coveredCpG = residuals(
    lm(myMetaData$Nbr_hypermethCpG ~ myMetaData$Nbr_coveredCpG))
  ## REORDER myMetaData by sample ID
  myMetaData = myMetaData[order(as.numeric(gsub("S", "", myMetaData$SampleID))),]
  return(myMetaData)
}

fullMetadata_OFFS  <- mycalcRMS(uniteCovALL_G2_woSexAndUnknowChr, fullMetadata_OFFS)

# Calculate number of methylated/unmethylated CpGs
df = percMethylation(uniteCovALL_G2_woSexAndUnknowChr) %>% data.frame()
df$CpG = paste(uniteCovALL_G2_woSexAndUnknowChr$chr, uniteCovALL_G2_woSexAndUnknowChr$start)
df = melt(df)
df = df %>% dplyr::group_by(variable) %>% dplyr::summarise(nCpG=n(),
                                                           nMethCpG=sum(value>70),
                                                           nUnmethCpG=sum(value<30)) %>% data.frame()
names(df)[names(df) %in% "variable"]="SampleID"

df=merge(df, fullMetadata_OFFS)

merge(df[c("SampleID", "nMethCpG", "nCpG")],
      fullMetadata_OFFS[c("SampleID", "Nbr_hypermethCpG", "Nbr_coveredCpG")])

mylm=lmerTest::lmer(nMethCpG ~ Sex * PAT * outcome + (1|brotherPairID), data = df)
step(mylm, reduce.random = F) 
# Model found: nMethCpG ~ Sex + PAT + (1 | brotherPairID)
#                   Eliminated   Sum Sq  Mean Sq NumDF   DenDF F value    Pr(>F)   
# Sex                    0 27651103 27651103     1 103.431 19.6174 2.355e-05 ***
# PAT                    0  7962282  7962282     1 101.150  5.6489   0.01935 * 


mylm=lmerTest::lmer(Nbr_hypermethCpG ~ Sex * PAT * outcome + (1|brotherPairID), data = fullMetadata_OFFS)
step(mylm, reduce.random = F) 

mylm=lmerTest::lmer(nMethCpG ~ Sex + PAT + (1|brotherPairID), data = df)
step(mylm, reduce.random = F) 

P1=plot(ggpredict(mylm, terms = c("PAT", "Sex")), add.data=T)+
  ylab("Number of hyper-methylated CpGs")+
  xlab(NULL)+
  ggtitle("Predicted and true number of hyper-methylated CpGs")

mylm=lmerTest::lmer(nUnmethCpG ~ Sex * PAT * outcome + (1|brotherPairID), data = df)
step(mylm, reduce.random = F) 
P2=plot(ggpredict(mylm, terms = c("PAT", "Sex")), add.data=T)+
  ylab("Number of hypo-methylated CpGs")+
  xlab(NULL)+
  ggtitle("Predicted and true number of hypo-methylated CpGs")
# not significative for unmethylated

pdf(file = "../../dataOut/figures/figB_globalMethDiffSexG1.pdf", width = 10, height = 5)
ggarrange(P1, P2, ncol=2, common.legend = T, legend="right")
dev.off() 

### Males have a lower global methylation than females (residuals of nbr of methylated sites by nbr of sites covered)
mod = lm(res_Nbr_methCpG_Nbr_coveredCpG ~ Sex, data = fullMetadata_OFFS)
summary(step(mod)) # sex is significant p = 0.000157 ***
anova(mod)

plot(ggpredict(mod, terms = c("Sex")), add.data = T) +
  xlab(NULL)+
  ylab("Residuals of N methylated sites on N covered sites") +
  ggtitle("Predicted values of global methylation in offspring")
