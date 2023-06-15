## The raw exome data was processed through the our exome pipeline and the output from this pipeline was a list of somatic mutations that were identified by comparison to the exome data from the blood.

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Exome_data")
mut1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Exome_data/mutations.snvs.indels.csv")
dim(mut1)
# [1] 647  58

## I would like to get all of the mutations that are a real mutation by selecting the mutations that have "ourJudgment" column =yes.
mut1a = mut1[mut1[, 38] == "yes", ]
dim(mut1a)
# [1] 108  58
## There are 108 mutations in total that according to the Costello lab algorithm are genuine mutations.

unique(mut1a$gene)
# [1] ABCC2    ACCS     ANKRD18B BRD7     C2       CALHM1   CSAG1    CTNNA3   FAM193B  GDPD2    GPR173
# [12] HIST1H4K IAPP     IDH1     IMMT     KRTAP5-3 KRTAP5-4 MUC16    MUC21    PASD1    PEX1     PHF8
# [23] PLA2G4D  SCN10A   SFTPA1   SIRT5    SUSD4    TBCC     TEC      TP53     TRIL     CRIPAK   CTSE
# [34] FAM60A   GPR126   NOL7     NTRK1    PABPC4L  PKD2     RECQL    RTN4     SCAPER   SLC20A2  SPEF2
# [45] WASL     XDH      APOBR    ATRX     C14orf39 CCDC168  CTNND2   CYSLTR2  DBC1     IL7R     MORC1
# [56] OR52N2   OR7C2    PCLO     SIGLEC5  SLAMF1   ZFYVE16
## There are 61 unique genes that are mutated in all three of the samples. Many of these mutations are likely to be the same across all of the samples. 
## Some of these mutations are actually silent somatic mutations so they can be removed from the list. Also remove genes whose type is unknown. 
## ANKRD18B is the only mutation with unknown type and I checked this muation and it is actually silent.

mut1a = mut1a[!(mut1a$type == "unknown"), ]
dim(mut1a)
# [1] 86 58

unique(mut1a$gene)
# [1] ABCC2    ACCS     BRD7     C2       CSAG1    CTNNA3   FAM193B  GDPD2
# [9] GPR173   HIST1H4K IAPP     IDH1     KRTAP5-3 MUC16    MUC21    PASD1
# [17] PEX1     PHF8     PLA2G4D  SUSD4    TBCC     TEC      TP53     TRIL
# [25] CRIPAK   FAM60A   GPR126   NOL7     NTRK1    PKD2     RECQL    RTN4
# [33] SCAPER   SLC20A2  SPEF2    WASL     ATRX     C14orf39 CCDC168  CTNND2
# [41] DBC1     IL7R     MORC1    OR52N2   OR7C2    PCLO     SLAMF1   ZFYVE16

## This leaves 48 non-synonamous somatic mutations.
length(unique(mut1a$gene))
# [1] 48

## I renamed samples so from their original section names so that they are sequentially oredered.

## I would like to narrow down the list to just those mutations that have a probe detected above background in at least one of the 69 tumor sections on the microarray:
dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_microarray_data/2013-140 sample probe profile.csv")

## Read in the original SampleInformation file that contains the missing samples (2 references and section 48)
datSample = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor/SF9495_SampleInfo_for_SampleNetwork.csv")

## I need to recalculate the CountDetPval and AvgDetPval to take account of the 3 samples that were removed during SampleNetwork (2 reference samples and section #48:
## Remove these samples from the expression data#.
## Take the average signal for remaining  samples
## Calculate the average detection P value by taking the next column over from avg signal (seq1+1)
## Now we ask what fraction of transcripts are detected at P < .05 in at least one sample:

CountDetPval = apply(dat1[, seq1 + 1], 1, function(x) {
    length(x[x < .05])
})

## Now I have recalculated the CountDetPval and AvgDetPval for the 69 samples

## Clear memory:
rm(list = ls()
   
## Now I can start to ask what the detection p Values are for each of the mutations.
## There are 44 mutated genes that are detected above background in at least 1 section.
## Get the variant frequency for each of the mutations in each of the three exomes...

## The next thing to do is to incorporate the amplicon Sanger sequencing data into the exome analysis. The majority of these 44 mutations were Sanger sequenced and a reasonably high proportion were found to also be present in the blood and were therefore false positives. I manually removed the mutations that were found to be false positives by Sanger sequencing and this file can be used to subset this data down. I also want to create a master table of mutation information based on this validated list of mutations:
## Read in the newly generated data frame with mutation and the number of samples it was called in as well as how many samples it was detected above background in.
## Restrict the list of confident mutations to the 33 which I am very confident about (called in exome data and look positive in IGV and/or Sanger validated).
## Looking at the known variant status, it doesn't make sense because it is saying that IDH1 R132H is novel when it is one of the most common glioma variants. For this reason I should leave it out of this table.

## Write out the comprehensive table:
   
setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Tables_for_paper")
write.table(dat2d, file = "Complete_mutation_information.csv", col.names = TRUE, row.names = FALSE, sep = ",")

