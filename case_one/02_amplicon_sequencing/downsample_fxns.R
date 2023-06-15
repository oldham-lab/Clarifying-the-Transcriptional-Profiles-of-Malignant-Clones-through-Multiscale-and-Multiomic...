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

downsampleCoverage <- function(projectName,
                               readsPerSample, ## Paths to reads per sample (one file per sample)
                               whichAmps,
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
      downsample <- foreach(seq(nResamples), .combine = rbind) %dorng% {
        ## Initialize downsample1 with first amp
        workingAmp <- workingSample[workingSample[, 1] %in% whichAmps[1],]
        replace <- FALSE
        if (nrow(workingAmp) < coverage)
          replace <- TRUE
        downsample1 <- sample_n(workingAmp, coverage, replace = replace) 
        ## Rest amps
        for (amp in whichAmps[2:length(whichAmps)]) {
          workingAmp <- workingSample[workingSample[, 1] %in% amp,]
          replace <- FALSE
          if(nrow(workingAmp) < coverage)
            replace <- TRUE
          downsample1 <- rbind(downsample1, sample_n(workingAmp, coverage, replace = replace)) ## Sample x reads per amplicon
        }
        ## Summarize number of alt and ref reads per amplicon:
        if(sum(downsample1$type %in% "alt") != 0) {
          downsample1 <- merge(as.data.frame(table(downsample1$gene[downsample1$type == "alt"])), as.data.frame(table(downsample1$gene[downsample1$type == "ref"])), by = 1, all = TRUE)
        } else {
          downsample1 <- merge(data.frame(whichAmps, NA), as.data.frame(table(downsample1$gene[downsample1$type == "ref"])), by = 1, all = TRUE)
        }
        downsample1[is.na(downsample1)] <- 0
        colnames(downsample1) <- c("Amplicon", "Alt", "Ref")
        return(downsample1)
      } # end downsamples
      downsampledCoverage[[sample]] <- downsample
      rm(downsample)
      collectGarbage()
    } # end samples
    assign(paste0("downsample_coverage_", coverage, "x"), downsampledCoverage)

  } ## // end coverage intervals
  save(list = ls()[grep("downsample_", ls())], file = paste0(outDir, "downsampleCoverage_", projectName, "_", nResamples, "_resamples.RData"))

} ## // end fxn

