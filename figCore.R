library(tidyverse)
library(compiler)
library(ggrepel)
library(grid)
library(ggcorrplot)
options(stringsAsFactors = F)
mainFolder <- "D:/Teams"
cwd <- paste0(mainFolder, "/Figures")
randcompiled <- paste0(mainFolder, "/Final_Results/Compiled_data")
randRaw <- "D:\\Teams\\Final_Results\\Raw_Data"
WTcoherence <- "D:\\Teams\\Final_Results\\WT_NPerts"
edgeDel <- paste0(mainFolder, "/Final_Results/edge_deletions")
QC <- "D:\\Teams\\Final_Results\\Quantitative_convergence"
RACIPE <- "D:\\Teams\\Final_Results\\RACIPE"
RACIPE_OE <- paste0(RACIPE, "\\OverExpression")
RACIPE_WT <- paste0(RACIPE, "\\WildType")
RACIPE_DE <- paste0(RACIPE, "\\DownExpression")
RACIPE_QC <- paste0(RACIPE, "\\QuantitativeConvergence")
NodeDel <- "D:\\Teams\\Final_Results\\NodePert"
netList <- c("EMT_RACIPE", "EMT_RACIPE2", "SIL", "SIL2", "melanoma", "SCLC", "EMT_MET_reduced", "drosophila")
EMPNets <- netList[c(1:4, 7)]
groupStrengthCausation <- paste0(mainFolder, "\\Final_Results\\GsMultiPert")

Controls <- c(EMPNets, "drosophila")
labelKeys <- c("maxFrust", "minFrust", "meanFrust", "maxNetFrust", "minNetFrust",
               "meanNetFrust","maxCoh", "minCoh", "meanCoh", "maxNetCoh", "minNetCoh",
               "meanNetCoh", "corFreqFrust", "pFreqFrust", "corFreqCoh", "pFreqCoh", 
               "corFrustCoh", "pFrustCoh"   )
labelvals <- c("Maximum Frustration", "Minimum Frustration", "Mean Frustration", 
               "Max Net Frustration", "Min Net Frustration", "Mean Net Frustration",
               "Maximum Coherence", "Minimum Coherence", "Mean Coherence", 
               "Max Net Coherence", "Min Net Coherence", "Mean Net Coherence",
               "\u03c1 Frequency-frustration", "pVal Frequency-frustration", 
               "\u03c1 Frequency-coherence", "pVal Frequency-coherence",
               "\u03c1 Coherence-frustration", "pVal Coherence-frustration")
names(labelvals) <- labelKeys

labelshorts <- c("Max Frust", "Min Frust", "Mean Frust", 
               "Max Net Frust", "Min Net Frust", "Mean Net Frust",
               "Max Coh", "Min Coh", "Mean Coh", 
               "Max Net Coh", "Min Net Coh", "Mean Net Coh",
               "\u03c1 Freq-Frust", "pVal Freq-Frust", 
               "\u03c1 Freq-Coh", "pVal Freq-Coh",
               "\u03c1 Coh-Frust", "pVal Coh-Frust")
names(labelshorts) <- labelKeys

netNameKey <- c(EMT_RACIPE = "22N 82E", EMT_RACIPE2 = "26N 100E", 
                SIL = "18N 33E", SIL2 = "20N 40E", EMT_MET_reduced = "57N 113E", 
                drosophila = "drosophila")
racEMPNets <- c("EMT_RACIPE", "EMT_RACIPE2", "silviera", "silviera2", "EMT_MET_reduced")
netNameKeyRAC <- c(EMT_RACIPE = "22N 82E", EMT_RACIPE2 = "26N 100E", 
                silviera = "18N 33E", silviera2 = "20N 40E", EMT_MET_reduced = "57N 113E")
phenKey <- c(E = "Epithelial", H = "Hybrid", M = "Mesenchymal")
source(paste0(cwd,"/plot_theme.R"))
source(paste0(cwd,"/utils.R"))
source(paste0(mainFolder, "\\Final_Results\\codes\\inflMat.R"))
source(paste0(mainFolder, "\\Final_Results\\codes\\groupPlotter.R"))
