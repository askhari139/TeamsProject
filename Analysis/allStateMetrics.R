source("D:/Teams/TeamsPaperFigures/figCore.R")
source(paste0(mainFolder, "/Final_Results/codes/MultiNodePert.R"))
source(paste0(mainFolder, "/Final_Results/codes/stateLabeller.R"))


metricWriterRand <- function(net, nCores = 8)
{
    setwd(paste0(randRaw, "/", net))
    topoFiles <- list.files(".",".topo")
    nets <- topoFiles %>% str_remove(".topo")
    write(net, file = paste0(net, "MetricLog.txt"), append = T)
    netsDone <- readLines(paste0(net, "MetricLog.txt"))
    nets <- nets[-which(nets %in% netsDone)]
    plan(multisession, workers = nCores)
    ### coherence calculations----
    dummy <- future_lapply(nets, function(p){
        topoFile <- paste0(p, ".topo")
        inflFile <- paste0("Influence/", p, ".csv")
        if (!file.exists(inflFile))
            influence_matrix(topoFile)
        dfCoherence <- coherenceSingleNodeAlt(p)
        dfScores <- labeller(topoFile)
        colsOriginal <- colnames(dfScores)[1:11]
        # dfOriginal <- dfScores %>% select(all_of(colsOriginal))
        # dfScores <- dfScores %>% select(-all_of(colsOriginal))
        dfCoherence <- dfCoherence %>% select(all_of(colsOriginal), coherence0)
        df <- merge(dfScores, dfCoherence, by= colsOriginal, all = T)
        write_csv(df, paste0(p, "_finFlagFreq.csv"), quote = "none")
        write(p, file = paste0(net, "MetricLog.txt"), append = T)
    }, future.seed = T)
    print(net)
    future:::ClusterRegistry("stop")
}
sapply(EMPNets[3:5], metricWriterRand)

metricWriterGs <- function(net, nCores = 8)
{
    setwd(groupStrengthCausation)
    setwd(net)
    topoFiles <- list.files(".",".topo")
    nets <- topoFiles %>% str_remove(".topo")
    write(net, file = paste0(net, "MetricLog.txt"), append = T)
    netsDone <- readLines(paste0(net, "MetricLog.txt"))
    netsDone <- readLines(paste0(net, "MetricLog.txt"))
    nets <- nets[-which(nets %in% netsDone)]
    plan(multisession, workers = nCores)
    ### coherence calculations----
    dummy <- future_lapply(nets, function(p){
        topoFile <- paste0(p, ".topo")
        inflFile <- paste0("Influence/", p, ".csv")
        if (!file.exists(inflFile))
            influence_matrix(topoFile)
        dfCoherence <- coherenceSingleNodeAlt(p)
        dfScores <- labeller(topoFile)
        colsOriginal <- colnames(dfScores)[1:11]
        # dfOriginal <- dfScores %>% select(all_of(colsOriginal))
        # dfScores <- dfScores %>% select(-all_of(colsOriginal))
        dfCoherence <- dfCoherence %>% select(all_of(colsOriginal), coherence0)
        df <- merge(dfScores, dfCoherence, by= colsOriginal, all = T)
        write_csv(df, paste0(p, "_finFlagFreq.csv"), quote = "none")
        write(p, file = paste0(net, "MetricLog.txt"), append = T)
    }, future.seed = T)
    print(net)
    future:::ClusterRegistry("stop")
}
metricWriterGs("EMT_RACIPE")

### Bimodality coefficients ----
dfAll <- lapply(EMPNets, function(net){#browser()
    setwd(randRaw)
    
    setwd(net)
    ffFils <- list.files(".", "_finFlagFreq")
    df <- sapply(ffFils, function(x){
        d <- read.csv(x) %>%
            filter(flag == 1, !is.na(Avg0))
        ssBm <- bmCoeff(d$Avg0)
        frustBm <- bmCoeff(d$frust0)
        cohBm <- NA
        if (any(colnames(d) == "coherence0"))
            cohBm <- bmCoeff(d$coherence0 %>% na.omit)
        
        return(c(ssBm, frustBm, cohBm))
    }) %>% t %>% data.frame %>% 
        set_names(c("SSF" ,"Frustration", "Coherence"))
    setwd("..")
    print(net)
    d1 <- df
    colnames(d1) <- paste0("bm", colnames(d1))
    d1 <- d1  %>%
        mutate(Net = str_remove(ffFils, "_finFlagFreq.csv"))
    d2 <- read.csv(paste0(randcompiled, "/", net, "_ALL.csv"))
    d2 <- merge(d2, d1, by= "Net", all = T)
    write.csv(d2, paste0(randcompiled, "/", net, "_ALL.csv"))
    # df
})


### group strengths ----

dfGroups <- lapply(EMPNets, function(net){
    setwd(randRaw)
    setwd(net)
    topoFiles <- list.files(".", ".topo")
    groups <- sapply(topoFiles, groupCalc) %>% t %>% data.frame %>%
        set_names(c("G11", "G22", "G12", "G21", "Net")) %>%
        mutate(across(c("G11", "G22", "G12", "G21"), as.numeric)) %>%
        mutate(Gs = (abs(G11) + abs(G22) + abs(G21) + abs(G12))/4)
    write.csv(groups, paste0(randcompiled, "/", net, "_rand_groups.csv"), row.names = F)
    write.csv(groups, paste0(net, "_rand_groups.csv"), row.names = F)
})
