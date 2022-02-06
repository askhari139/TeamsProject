source("D:\\Teams\\TeamsPaperFigures/figCore.R")
source("D:\\Teams\\Final_Results\\codes\\stateLabeller.R")

library(future)
library(future.apply)
# library(ggcorrplot)


randCorrel <- function(net)
{
    netDir <- paste0(randRaw, "/", net)
    wd <- getwd()
    setwd(netDir)
    topoFiles <- list.files(".", ".topo")
    dists <- sapply(topoFiles, function(topoFile)
    {
        # print(topoFile)
        inflFile <- paste0("Influence/", str_replace(topoFile, ".topo", "_fullMat.csv"))
        if (!file.exists(inflFile))
            return(NA)
        inflMat <- read.csv(inflFile, row.names = 1) %>% as.matrix
        
        corMat <- correlationMatBool(topoFile)
        if (is.na(corMat))
            return(NA)
        inflMat <- inflMat[rownames(corMat), colnames(corMat)]
        
        sqrt(sum((inflMat - corMat)^2))/(2*length(corMat))
    })
    df <- data.frame(Network = topoFiles %>% str_remove(".topo"), Dist = dists)
    write.csv(df, paste0(net,"_correlInfl.csv"), row.names = F)
}
randCorrel <- cmpfun(randCorrel)
plan(multisession, workers= 5)
future_lapply(EMPNets, randCorrel)