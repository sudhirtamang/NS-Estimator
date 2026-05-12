# ==============================================================================
# PROJECT:  NS-Estimator
# SCRIPT:   file2.R
# PURPOSE:  BIAS correction on the NS-Estimator
# AUTHOR:   Sudhir T.
# DATE:     12th May, 2026 
# ==============================================================================

rm(list = ls())
library(Tlasso)
library(tensr)
library(glasso)
library(expm)
library(rTensor)
library(doParallel)

source("C:/Users//sudhi//Desktop//Fast and Separable Estimation Replication//replication//Model 1//Separate.fit.R")
source("C://Users//sudhi//Desktop//Fast and Separable Estimation Replication//replication//Model 1//simulation.summary.R")
source("C://Users//sudhi//Desktop//Fast and Separable Estimation Replication//replication//Model 1//Model7.R")
# source("C:/Users/sudhi/Desktop/PhD projects/NS-Estimator/Lfunctions.R")
source("Lfunctions.R")