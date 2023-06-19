library(data.table)
library(dplyr)
library(WGCNA)
library(openblasctl)
library(doFuture)
library(doRNG)
library(flexiblas)
flexiblas_load_backend("OPENBLASPTHREAD")
flexiblas_switch(2)
options(future.globals.maxSize = +Inf)
openblas_set_num_threads(1)

downsampleSelfCov <- function(projectName,
                              readsPerSample,
                              whichAmp,
                              coverageIntervals,
                              nResamples,
                              nThreads,
                              outDir) {
  registerDoFuture()
  plan(multicore, workers = nThreads)
  for (coverage in coverageIntervals) {
    cat("\n--- Starting coverage =", coverage,"---\n")
    downsampledCoverage <- list()
    for (sample in readsPerSample) {
      cat("\nStarting", sample, "\n")
      workingSample <- as.data.frame(fread(sample, nThread = nThreads))
      workingSample <- workingSample[workingSample[, 1] %in% whichAmp,]
      downsample <- foreach(seq(nResamples), .combine = cbind) %dorng% {
        idx <- sample(seq(nrow(workingSample)), coverage*2) 
        idx1 <- sample(idx, coverage)
        idx2 <- idx[!idx %in% idx1]
        downsample1 <- workingSample[idx1,]
        downsample2 <- workingSample[idx2,]
        if(sum(downsample1$type %in% "alt") != 0) {
          downsample1 <- merge(as.data.frame(table(downsample1$gene[downsample1$type == "alt"])), as.data.frame(table(downsample1$gene[downsample1$type == "ref"])), by = 1, all = TRUE)
        } else {
          downsample1 <- merge(data.frame(whichAmp, NA), as.data.frame(table(downsample1$gene[downsample1$type == "ref"])), by = 1, all = TRUE)
        }
        if(sum(downsample2$type %in% "alt") != 0) {
          downsample2 <- merge(as.data.frame(table(downsample2$gene[downsample2$type == "alt"])), as.data.frame(table(downsample2$gene[downsample2$type == "ref"])), by = 1, all = TRUE)
        } else {
          downsample2 <- merge(data.frame(whichAmp, NA), as.data.frame(table(downsample2$gene[downsample2$type == "ref"])), by = 1, all = TRUE)
        }
        downsample1[is.na(downsample1)] <- 0
        downsample2[is.na(downsample2)] <- 0
        colnames(downsample2) <- colnames(downsample1) <- c("Amplicon", "Alt", "Ref")
        return(rbind(downsample1, downsample2))
      } # end downsamples
      downsampledCoverage[[sample]] <- downsample
      rm(downsample)
      collectGarbage()
    } # end samples
    assign(paste0("downsample_self_coverage_", coverage, "x"), downsampledCoverage)
  } ## // end coverage intervals
  save(list = ls()[grep("downsample_self_coverage_", ls())], file = paste0(outDir, "downsampleSelfCov_", projectName, "_", nResamples, "_resamples.RData"))
} ## downsampleSelfCov()
