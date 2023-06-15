## I copied all of the raw idat files into a single combined folder:
# cp -r /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_methylation_data_raw/2015-9159_Oldham_Meth_450/UNGC_2015-9159/2015-9159_Oldham_Meth_450/2015-9159_iScans/3999119004/*.idat /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_methylation_data_raw/2015-9159_Oldham_Meth_450/UNGC_2015-9159/2015-9159_Oldham_Meth_450/2015-9159_iScans/Combined

# cp -r /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_methylation_data_raw/2015-9159_Oldham_Meth_450/UNGC_2015-9159/2015-9159_Oldham_Meth_450/2015-9159_iScans/3999119017/*.idat /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_methylation_data_raw/2015-9159_Oldham_Meth_450/UNGC_2015-9159/2015-9159_Oldham_Meth_450/2015-9159_iScans/Combined

# cp -r /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_methylation_data_raw/2015-9159_Oldham_Meth_450/UNGC_2015-9159/2015-9159_Oldham_Meth_450/2015-9159_iScans/3999119030/*.idat /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_methylation_data_raw/2015-9159_Oldham_Meth_450/UNGC_2015-9159/2015-9159_Oldham_Meth_450/2015-9159_iScans/Combined

# cp -r /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_methylation_data_raw/2015-9159_Oldham_Meth_450/UNGC_2015-9159/2015-9159_Oldham_Meth_450/2015-9159_iScans/3999119041/*.idat /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_methylation_data_raw/2015-9159_Oldham_Meth_450/UNGC_2015-9159/2015-9159_Oldham_Meth_450/2015-9159_iScans/Combined

# cp -r /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_methylation_data_raw/2015-9159_Oldham_Meth_450/UNGC_2015-9159/2015-9159_Oldham_Meth_450/2015-9159_iScans/3999119064/*.idat /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_methylation_data_raw/2015-9159_Oldham_Meth_450/UNGC_2015-9159/2015-9159_Oldham_Meth_450/2015-9159_iScans/Combined

# cp -r /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_methylation_data_raw/2015-9159_Oldham_Meth_450/UNGC_2015-9159/2015-9159_Oldham_Meth_450/2015-9159_iScans/3999119096/*.idat /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_methylation_data_raw/2015-9159_Oldham_Meth_450/UNGC_2015-9159/2015-9159_Oldham_Meth_450/2015-9159_iScans/Combined

## Load ChAMP package

.libPaths("/fast-data/R-3.1.1-packages-SJS")
library(ChAMP)

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_methylation_data_raw/2015-9159_Oldham_Meth_450/UNGC_2015-9159/Processed_data")

## datSample1 has the methylation data.
## datSample2 has the sample information for the samples that were analyzed for gene expression and exome sequencing. This includes 2 reference samples and #48 which was removed as an outlier
## datSample3 contains the information about the concentration of the
## Load the well ID info:
datSample4 = read.csv("/big-data/Sam/Sam GBM dataset/SF9495_methylation_data_081215/UNGC 2015-9159/2015-9159 Oldham Meth 450/UNGC_complete_data_sheet.csv")

dim(datSample1)
# [1] 72  3
dim(datSample2)
# [1] 72 14
dim(datSample3)
# [1] 69 20
dim(datSample4)
# [1] 69  7

## To start with I will exclude duplicated samples:
datSample1 = datSample1[1:69, ]

## Clean up datSample2. Remove the two reference samples that are at the end.
datSample2 = datSample2[1:70, ]
datSample2 = datSample2[!datSample2$Original == 48, ]

dim(datSample1)
# [1] 69  3
dim(datSample2)
# [1] 69 14

## Check that they are all in the same order:
all.equal(datSample1$External.Sample.ID, datSample2$Original)
# [1] TRUE
all.equal(datSample2$Original, datSample3$Original)
# [1] TRUE
all.equal(datSample3$Original, datSample4[, 4])
# [1] TRUE

## Now combine the data into a single sample information file. It is not worth incorporating anything related to RNA in this table since it doesn't affect gDNA. Anything to do with DNA should be included.
## Now convert Isolation.Batch and QC.Batch into actual batches as they were done in groups of 12 and there were 86 sections actually cut.
86 / 12
# [1] 7.166667
## There were actually 8 batches in total. I now need to convert these into batches.
## I will use .bincode function to assign samples to a batch. Create break points for the borders. I will go upto 108 which represents 9 batches because it gives me 97 which would be the cutoff for batch 8 otherwise the highest samples in batch 8 would be NA.
Batch = c(seq(1, 108, 12))
Isolation.Batch = .bincode(datSample$Isolation.Batch, Batch, FALSE)
QC.Batch = .bincode(datSample$QC.Batch, Batch, FALSE)

## Now replace these in datSample:
datSample$Isolation.Batch = Isolation.Batch
datSample$QC.Batch = QC.Batch

## Now renames the samples in sequential order and create a SampleLabel to match what I already have for the gene expression.
SampleLabel = paste("SF9495_", c(1:69), sep = "")
datSample = data.frame(Sample_Name = SampleLabel, datSample, New.Sample.No = c(1:69), Group = "ALL")

## It looks like ChAMP needs a column called "Basename" which is essentially the file location and basename. I will create this now"
Basename = paste("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_methylation_data_raw/2015-9159_Oldham_Meth_450/UNGC_2015-9159/2015-9159_Oldham_Meth_450/2015-9159_iScans/Combined", datSample$Corelabel, sep = "/")
datSample = data.frame(datSample, Basename = Basename)

## Write the table out:
write.table(datSample, file = "Complete_SampleInfo_For_Methylation.csv", col.names = TRUE, row.names = FALSE, sep = ",")

## It looks like the sampleInfo needs the following columnnames Sample_Name, Sample_Well, Sample_Plate, Sample_Group, Pool_ID, Slide, Array as phenotype data. It will not take any additional data or any column names apart from these. If I want to add additional covariates then these need to added as studyInfo files.
datSampleChAMP = data.frame(Sample_Name = datSample[, 1], Sample_Well = datSample[, 14], Sample_Plate = "Plate1", Sample_Group = "All", Pool_ID = NA, Slide = datSample[, 10], Array = datSample[, 11], Basename = Basename)

## I would also like to create a studyInfo.txt file that has other batch information.
studyInfo = data.frame(Sample_Name = datSample$Sample_Name, datSample[, 3:9])

## Now write out these tables to disk.
setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_methylation_data_raw/2015-9159_Oldham_Meth_450/UNGC_2015-9159/2015-9159_Oldham_Meth_450/2015-9159_iScans/Combined")
write.table(datSampleChAMP, file = "sampleInfo_For_ChAMP.csv", col.names = TRUE, row.names = FALSE, sep = ",")

## StudyInfo is used for SVD in ChAMP to address batch effects. When using SampleNetwork this file is not required.
write.table(studyInfo, file = "studyInfo.txt", col.names = TRUE, row.names = FALSE, sep = "\t")

rm(list = ls())

## I would now like to try to process the data from the raw IDAT files through probe filtering and BMIQ normalization. I will try using ChAMP for this.

.libPaths("/fast-data/R-3.1.1-packages-SJS")
library(ChAMP)

## Load the data using the champ.load function. This uses the current directory as default. I had to copy the SampleInformation file into the directory containing the idat files. This function used the minfi function to load the data and automatically filters out 65 SNP probes that are included on the chip as internal controls:
myLoad = champ.load()
## champ.load worked and this is the output of the function:

# Filtering probes with a detection p-value above 0.01 in more than one sample has removed 11799 probes from the analysis. If a large number of probes have been removed, ChAMP suggests you look at the failedSample.txt file to identify potentially bad samples.
# Filtering probes with a beadcount <3 in at least 5% of samples, has removed 760 from the analysis.
# Zeros in your dataset have been replaced with 0.000001
# Cluster image is not saved when the number of samples exceeds 60.
# The analysis will proceed with 461797 probes and 69 samples.

## This function created a new folder called resultsChamp and it output some metrics for the data and also indicated the fraction of failed probes in each sample.

## The next thing to do is the normalization of the data. The default setting for ChAMP is BMIQ which accounts for type II bias by normalizing the type I and type II probes seperately. This creates figures for each sample showing the normalization of type I and type II probes and creates a new folder called "Normalization".

myNorm = champ.norm()

## I think that at this stage it makes sense to analyze the data with SampleNetwork to assess technical batch effects.

dat1 = data.frame(PROBE_ID = rownames(myNorm$beta), myNorm$beta)

write.table(dat1, file = "BMIQ_normalized_data.csv", col.names = TRUE, row.names = FALSE, sep = ",")
write.table(dat1, file = "BMIQ_normalized_data.txt", col.names = TRUE, row.names = FALSE, sep = "\t")
