setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor")
## Read in the sample information file:
datSample = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_SampleInfo/SF9495_SampleInfo_for_SampleNetwork.csv")

## Add the array ID and position so I can correct for any batch effects. 

## Now I will load the expression file and relabel the column names. I opened the file 2013-140_Sample probe profile.xlsx and resaved it as csv:
dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_microarray_data/2013-140 sample probe profile.csv")
dim(dat1)
# [1] 47202   302

## Take the average signal for each of these samples:
seq1 = seq(3, 287, 4)

## Calculate the average detection P value by taking the next column over from avg signal (seq1+1):
AvgDetPval = apply(dat1[, seq1 + 1], 1, mean)

## Now we ask what fraction of transcripts are detected at P < .05 in at least one sample:
CountDetPval = apply(dat1[, seq1 + 1], 1, function(x) {
    length(x[x < .05])
})

## Now to run the SampleNetwork Function. I ran this locally due to interactivity needed:

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor")
dat1 = read.csv("SF9495_ExpressionFile_for_SampleNetwork.csv")
datSample = read.csv("SF9495_SampleInfo_for_SampleNetwork.csv")
source("/fast-data/Shared/Code/SJS/SampleNetwork_1.01_cluster_SJS.r")
indexAll = c(7:78)

## Color by QC batch (col3):

SampleNetwork(
    datExprT = dat1,
    method1 = "correlation",
    impute1 = FALSE,
    subset1 = NULL,
    skip1 = 6,
    indices1 = list(indexAll),
    subgroup1 = 3,
    sampleinfo1 = datSample,
    samplelabels1 = 16,
    grouplabels1 = 19,
    fitmodels1 = TRUE,
    whichmodel1 = "univariate",
    whichfit1 = "pc1",
    btrait1 = c(1, 3, 5, 6, 7, 8, 9, 10, 14, 17, 18),
    trait1 = NULL,
    asfactors1 = c(3, 9, 17, 18),
    projectname1 = "SF9495",
    cexlabels = 0.7,
    normalize1 = TRUE,
    replacenegs1 = FALSE,
    exportfigures1 = TRUE,
    verbose = TRUE
)

## After the first round the data looked pretty good and the two reference samples looked clearly different to the rest of the samples. I removed these by using <-2.
## After the second round there was one outlier (section 48) which I removed using <-3:
## No further samples removed, 3 samples in total were removed, the controls and 1 other sample:
## There were no further significant batch effects other than section number but since this represents variation between samples this is what I would expect and it will not be corrected for.
## After removing 1 outlier which was sample 48 this leaves me with 69 samples.
