mainFolder <- "D:/Teams"
plotFolder <- paste0(mainFolder, "/Figures/Fig5")
source(paste0(mainFolder, "/Figures/figCore.R"))
source(paste0(mainFolder, "/Final_Results/codes/stateLabeller.R"))
library(future)
library(future.apply)

### Fig 5A-B -----
coherenceHeatmaps <- function(net)
{#browser()
    setwd(plotFolder)
    if(!dir.exists("coherencePlots"))
        dir.create("coherencePlots")
    setwd("coherencePlots")
    freqFile <- paste0(randRaw, "/", net, "/", net, "_finFlagFreq.csv")
    multiCohFile <- paste0(WTcoherence, "/", net, "_coherence.csv")
    nodeCohFile <- paste0(randcompiled, "/", net, "_NodeStateCoherence.csv")
    dfFreq <- read.csv(freqFile)
    dfCoh <- read.csv(multiCohFile)
    dfSingleCoh <- read.csv(nodeCohFile)
    
    
    rownames(dfCoh) <- dfCoh$X
    colnames(dfCoh) <- colnames(dfCoh) %>% str_remove("X") %>% str_remove_all("\\.")
    dfCoh <- dfCoh[, -1] %>% t %>% data.frame %>% mutate(states = rownames(.))
    dfCoh$states <- paste0("'", dfCoh$states, "'")
    df <- merge(dfCoh, dfFreq, by = "states", all.x = T) %>%
        select(-contains("Avg"), -flag, -contains("frust"), -contains("SD"), -Score, -Partial,
               -Epithelial, -Mesenchymal)
    df$phenotype <- paste0(df$phenotype, 1:nrow(df))
    
    dfGat <- df %>% gather(key = "nPert", value = "Coherence", -states, -phenotype)
    dfGat$nPert <- dfGat$nPert %>% str_remove("X") %>% as.numeric
    dfGat$nPert <- dfGat$nPert/22
    # dfGat$Phenotype <- paste0(dfGat$phenotype, 1:nrow(dfGat))
    
    ggplot(dfGat, aes(x = nPert, y = phenotype, 
                      fill = Coherence)) + geom_tile() +
        scale_fill_viridis_c(limits = c(0,1)) + theme_Publication() +
        scale_x_continuous(expand=c(0,0)) + 
        scale_y_discrete(expand=c(0,0)) + 
        theme(legend.position = "right", legend.direction = "vertical",
              legend.key.height = unit(1, "cm"), 
              legend.key.width = unit(0.6,"cm"),
              legend.title = element_text(hjust = 0.5),
              panel.background = element_blank(), 
              plot.background = element_blank(), panel.grid = element_blank(), 
              axis.text.y = element_blank()) +
        labs(x = "Level of Perturbation", y = "", fill = "Mean\nCoherence")
    ggsave(paste0(net, "_coherenceMultiNode.png"), width = 6.5, height = 5)
    
    
    
    dfCoh <- read.csv(nodeCohFile)
    rownames(dfCoh) <- dfCoh$X
    colnames(dfCoh) <- colnames(dfCoh) %>% str_remove("X") %>% str_remove_all("\\.")
    dfCoh <- dfCoh[, -1] %>% t %>% data.frame %>% mutate(states = rownames(.))
    dfCoh$states <- paste0("'", dfCoh$states, "'")
    df <- merge(dfCoh, dfFreq, by = "states", all.x = T) %>%
        select(-contains("Avg"), -flag, -contains("frust"), -contains("SD"), 
               -Score, -Partial,
               -Epithelial, -Mesenchymal)
    df$phenotype <- paste0(df$phenotype, 1:nrow(df))
    
    dfGat <- df %>% gather(key = "Nodes", value = "Coherence", -states, -phenotype) %>%
        mutate(Nodes = Nodes %>% str_replace_all(regex("\\W+"), ""))
    
    setwd(paste0(randRaw, "/", net))
    l <- groupCalc1(paste0(net, ".topo"))
    nodesC <- l[[1]] %>% unlist
    nodes <- dfGat$Nodes %>% unique
    nodesS <- nodes[!(nodes %in% nodesC)]
    setwd(cwd)
    dfGat$Nodes <- dfGat$Nodes %>% factor(levels = c(nodesC, nodesS))
    ggplot(dfGat, aes(x = Nodes, y = phenotype, 
                      fill = Coherence)) + geom_tile() +
        scale_fill_viridis_c(limits = c(0,1)) + theme_Publication() +
        scale_x_discrete(expand=c(0,0)) + 
        scale_y_discrete(expand=c(0,0)) + 
        theme(legend.position = "right", legend.direction = "vertical",
              legend.key.height = unit(1, "cm"), 
              legend.key.width = unit(0.6,"cm"),
              legend.title = element_text(hjust = 0.5),
              panel.background = element_blank(), 
              plot.background = element_blank(), panel.grid = element_blank(), 
              axis.text.y = element_blank(), 
              axis.text.x = element_text(angle = 90, hjust = 1, size = rel(0.6))) +
        labs(x = "", y = "", fill = "Mean\nCoherence")
    ggsave(paste0(net, "_coherenceSingleNode.png"), width = 6.5, height = 5)
    
    print(net)
}
sapply(EMPNets, coherenceHeatmaps)


### Fig 5C -----
coherenceScatterPlots <- function(net)
{
    print(net)
    setwd(plotFolder)
    setwd("coherencePlots")
    freqFile <- paste0(randRaw, "/", net, "/", net, "_finFlagFreq.csv")
    multiCohFile <- paste0(WTcoherence, "/", net, "_coherence.csv")
    dfFreq <- read.csv(freqFile)
    dfCoh <- read.csv(multiCohFile)
    
    rownames(dfCoh) <- dfCoh$X
    colnames(dfCoh) <- colnames(dfCoh) %>% str_remove("X") %>% str_remove_all("\\.")
    dfCoh <- dfCoh[, -1] %>% t %>% data.frame %>% mutate(states = rownames(.))
    dfCoh$states <- paste0("'", dfCoh$states, "'")
    df <- merge(dfCoh, dfFreq, by = "states", all.x = T) %>%
        select(-contains("Avg"), -flag, -contains("frust"), -contains("SD"), -Score, -Partial,
               -Epithelial, -Mesenchymal)
    df$state <- paste0(df$phenotype, 1:nrow(df))
    df$phenotype <- ifelse(df$phenotype == "H", "Hybrid", "Terminal")
    
    dfGat <- df %>% gather(key = "nPert", value = "Coherence", -states, -phenotype, -state)
    dfGat$nPert <- dfGat$nPert %>% str_remove("X") %>% as.numeric
    dfGat$nPert <- dfGat$nPert/22
    dfGat %>% 
        group_by(state) %>% summarise(HalfMax = nPert[Coherence < 0.5][1]) %>%
        mutate(phenotype = ifelse(str_detect(state, "H"), "Hybrid", "Terminal")) %>%
        group_by(phenotype) %>% summarise(Avg = mean(HalfMax), Std = sd(HalfMax)) %>% 
        mutate(Network = net)
}

df <- lapply(EMPNets, coherenceScatterPlots) %>% reduce(rbind.data.frame) %>%
    mutate(Network = netNameKey[Network])
dfHyb <- df %>% filter(phenotype == "Hybrid") %>% select(-phenotype) %>% 
    set_names(c("Hybrid", "hStd", "Network"))
df1 <- df %>% filter(phenotype == "Terminal") %>% select(-phenotype, -Network) %>% 
    set_names(c("Terminal", "tStd")) %>% cbind.data.frame(dfHyb)

ggplot(df1, aes(x = Terminal, y = Hybrid, shape = Network)) +
    geom_point(size = 5) +
    geom_errorbar(aes(ymin = Hybrid - hStd, ymax = Hybrid + hStd)) +
    geom_errorbarh(aes(xmin = Terminal - tStd, xmax = Terminal + tStd)) +
    # geom_label(aes(label = Network)) +
    xlim(c(0,0.35)) + ylim(c(0, 0.35)) +
    geom_abline(slope = 1, intercept = 0) +
    theme_Publication() +
    theme(legend.position = c(0.2,0.8), legend.direction = "vertical") +
    labs(title = "Half-minimum perturbation")

ggsave("mutliNodeCoherence.png", width = 6.5, height = 5)


### Fig 5D-E - transition and distribution plots ----
meanPlotWT <- function(net, phenotype, df=NULL, nStates = 1)
{
    setwd(WTcoherence)
    
    if (is.null(df))
    {
        l <- list.files(".", paste0(net, "_allNodeCoherence"))
        df <- read_csv(l, show_col_types = F)
    }
    
    df <- df %>% filter(initPhen %in% phenotype) %>% mutate(Fraction = nNode/max(nNode)) %>%
        mutate(EMTScore = score_states(fin, paste0(net, ".topo")))
    inits <- unique(df$init)
    if (!is.null(nStates))
        inits <- inits[1:nStates]
    nStates <- length(inits)
    sapply(1:nStates, function(i){
        x <- inits[i]
        d <- df %>% filter(init == x) %>% group_by(Fraction) %>% 
            summarise(Hamming = sum(Hamming*Freq), 
                      Std = sum(Freq*(Hamming - sum(Hamming*Freq))^2), 
                      EMTScore = sum(EMTScore*Freq), 
                      StdScore = sum(Freq*(EMTScore - sum(EMTScore*Freq))^2), 
                      initPhen =  phenKey[unique(initPhen)])
        ggplot(d, aes(x = Fraction, y = Hamming, color = EMTScore)) + 
            geom_point(size = 3) + facet_wrap(~initPhen, ncol = 1) +
            # geom_errorbar(aes(ymin = Hamming - Std, ymax = Hamming + Std)) + 
            theme_Publication() + 
            scale_color_viridis_c(limits = c(-1,1))+
            xlim(0,1) +
            theme(legend.position = "right", legend.direction = "vertical",
                  legend.key.height = unit(0.8, "cm"), axis.text = element_text(size = rel(1.4)),
                  strip.text = element_text(size = rel(1.5))) +
            labs(x = "Perturbation", y = "Hamming")
        setwd(paste0(plotFolder, "/", "hammings/"))
        if(!dir.exists(net))
            dir.create(net)
        setwd(net)
        # browser()
        ggsave(paste0(net, "_", unique(d$initPhen), "_mean_", i, ".png"), 
               width = 6.5, height = 5.5)
    })
}

levelize <- function(x)
{
    y <- as.character(x)
    y[x<0.25] <- "Low"
    y[x>=0.25 & x<0.75] <- "Med"
    y[x>=0.75] <- "High"
    y
        
}
distributions <- function(net, phenotype, df = NULL, nStates = NULL, groupIt = T)
{
    if (is.null(df))
    {
        l <- list.files(".", paste0(net, "_allNodeCoherence"))
        df <- read_csv(l, show_col_types = F, lazy = F)
    }
    
    df <- df %>% filter(initPhen == phenotype) %>% mutate(Fraction = nNode/max(nNode)) %>%
        mutate(EMTScore = ifelse(finPhen == "H", 0, 1)) %>% 
        mutate(EMTScore = ifelse(finPhen == "E", -1, EMTScore)) %>%
        mutate(Level = levelize(Fraction))
    pertRanges <- c(0, 0.25, 0.75, 1)
    if (!groupIt)
        sapply(2:4, function(x){
            pertRange <- c(pertRanges[x-1], pertRanges[x])
            d <- filter(df , Fraction < pertRange[2], Fraction >= pertRange[1]) %>%
                mutate(Count = Freq*100 %>% round) %>% lapply(rep, .[["Count"]]) %>%
                as.data.frame
            inits <- d$init %>% unique
            if (!is.null(nStates))
                inits <- inits[1:nStates]
            sapply(1:length(inits), function(i){
                initS <- inits[i]
                ggplot(d %>% filter(init == initS) %>%
                           mutate(Phenotype =factor(finPhen, levels = c("E","H", "M"))), 
                       aes(x = Hamming, fill = Phenotype)) +
                    geom_histogram(aes(y = ..count../sum(..count..))) + 
                    theme_Publication() + 
                    scale_fill_manual(breaks = c("E", "H", "M"), 
                                      values = c("red", "green", "blue")) + 
                    labs(x = "Hamming", y = "Frequency") +
                    theme(legend.position = c(0.5, 0.9))+ 
                    ylim(0,0.7)
                setwd(paste0(plotFolder, "/", "hammings/"))
                if(!dir.exists(net))
                    dir.create(net)
                ggsave(paste0(plotFolder, "/hammings/", net, "/",net, "_", phenotype, "_", i, "_pert_", 
                              pertRange[2] %>% round(2), ".png"), width = 5.5, height = 5)
                setwd(WTcoherence)
            })
            
        })
    else
    {
        d <- df %>%
            mutate(Count = Freq*100 %>% round) %>% lapply(rep, .[["Count"]]) %>%
            as.data.frame
        inits <- d$init %>% unique
        if (!is.null(nStates))
            inits <- inits[1:nStates]
        sapply(1:length(inits), function(i){
            initS <- inits[i]
            d1 <- d %>% filter(init == initS) %>%
                mutate(Level = factor(Level, levels = c("Low", "Med", "High")))
            # %>%
            #     mutate(finPhen = phenKey[finPhen]) %>%
            #     mutate(Phenotype = factor(finPhen, 
            #                               levels = c("Epithelial", "Hybrid", "Mesenchymal")))
            
            ggplot(d1, aes(x = Hamming, fill = finPhen)) +
                facet_wrap(~Level, ncol = 3) +
                geom_histogram(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + 
                theme_Publication() + 
                scale_fill_manual(values = c("red", "green", "blue")) +
                labs(x = "Hamming", y = "Frequency") +
                ylim(c(0, 0.7))+
                theme(legend.position = "none", strip.text = element_text(size = rel(1.5)))
            setwd(paste0(plotFolder, "/", "hammings/"))
            if(!dir.exists(net))
                dir.create(net)
            ggsave(paste0(plotFolder, "/hammings/", net, "/",net, "_", 
                          phenotype, "_", i, ".png"), width = 13, height = 5)
            setwd(WTcoherence)
        })
    }
}

phenotypes <- c("E", "H", "M")
sapply(phenotypes, function(x){
    setwd(WTcoherence)
    # sapply(c(EMPNets), meanPlot, phenotype = x, nStates = NULL)
    sapply(c(EMPNets), distributions, phenotype = x)
})

net <- "EMT_RACIPE_rand_102"
setwd(WTcoherence)
df <- read_csv("EMT_RACIPE/EMT_RACIPE_rand_102_allNodeCoherence_nPert100_nIter10.csv")
sapply(phenotypes, function(x){
    distributions(net, df, phenotype = x, nStates = 3, groupIt = T)
})


### Fig 5E demo ----

meanPlotOverlap <- function(WT, rand, net, phenotype, df=NULL, nStates = 1)
{
    setwd(WTcoherence)
    df1 <- read_csv(WT, show_col_types = F)
    df2 <- read_csv(rand, show_col_types = F)
    randName <- str_extract(rand, "rand_\\d+?_all") %>% str_remove("_all") %>% 
        paste0(net, "_", .)
    WTfrust <- paste0(randRaw, "/", net, "/", net, "_finFlagFreq.csv") %>% read_csv(show_col_types = F)
    WTfrust <- WTfrust$frust0 %>% set_names(WTfrust$states)
    randss <- paste0(randRaw, "/", net, "/", randName, "_finFlagFreq.csv") %>% 
        read_csv(show_col_types = F)
    randfrust <- randss$frust0 %>% set_names(randss$states)
    
    
    df1 <- df1 %>% mutate(Frustration = WTfrust[fin]) %>%
        filter(initPhen %in% phenotype) %>% mutate(Fraction = nNode/max(nNode)) 
    inits1 <- unique(df1$init)
    df2 <- df2 %>% mutate(Frustration = randfrust[fin]) %>%
        filter(initPhen %in% phenotype) %>% mutate(Fraction = nNode/max(nNode))
    inits2 <- unique(df2$init)
    sapply(1:nStates, function(i){
        x1 <- inits1[i]
        x2 <- inits2[i]
        d1 <-  df1 %>% filter(init == x1) %>% group_by(Fraction) %>%
            summarise(Hamming_WT = sum(Hamming*Freq), 
                      Frustration_WT = max(Frustration))
        d2 <- df2 %>% filter(init == x2) %>% group_by(Fraction) %>%
            summarise(Hamming_rand = sum(Hamming*Freq),
                      Frustration_rand = max(Frustration))
        d <- merge(d1, d2, by = "Fraction", all = T)
        HammingDf <- d %>% select(Fraction, contains("Hamming")) %>% 
            gather(key = "Net", value = "Hamming", -Fraction) %>% 
            mutate(Net = str_remove(Net, "Hamming_"))
        Frustdf <- d %>% select(Fraction, contains("Frustration")) %>% 
            gather(key = "Net", value = "Frustration", -Fraction) %>% 
            mutate(Net = str_remove(Net, "Frustration_"))
        d <- merge(HammingDf, Frustdf, by = c("Fraction", "Net"), all = T) %>%
            mutate(initPhen = phenKey[phenotype])
        
        ggplot(d, aes(x = Fraction, y = Hamming, color = Frustration, shape = Net)) + 
            geom_point(size = 3) +
            facet_wrap(~initPhen, ncol = 1) +
            # geom_errorbar(aes(ymin = Hamming - Std, ymax = Hamming + Std)) + 
            theme_Publication() + 
            scale_color_viridis_c(limits = c(0, 0.37))+
            xlim(0,1) +
            theme(legend.position = "right", legend.direction = "vertical",
                  legend.key.height = unit(1, "cm"), 
                  axis.text = element_text(size = rel(1.4)),
                  strip.text = element_text(size = rel(1.5))) +
            labs(x = "Perturbation", y = "Hamming", shape = "")
        setwd(paste0(plotFolder, "/", "hammings/"))
        if(!dir.exists(net))
            dir.create(net)
        # browser()
        ggsave(paste0( net, "/",net, "_", unique(d$initPhen), "_randMean_", i, ".png"), 
               width = 6.5, height = 5.8)
        setwd(WTcoherence)
    })
}

meanPlotOverlap("EMT_RACIPE_allNodeCoherence_nPert100_nIter10_reps1.csv",
                "EMT_RACIPE/EMT_RACIPE_rand_102_allNodeCoherence_nPert100_nIter10.csv", 
                "EMT_RACIPE", "H", nStates = 3)


### Fig 5 F,G AUC ----
aocCalc <- function(df)
{
    if (ncol(df)<3)
        return()
    df <- df %>% select(initPhen, finPhen, Freq, nNode, init) %>% 
        mutate(Fraction = nNode/max(nNode)) #%>%
        # mutate(initPhen = ifelse(initPhen == "H", "Hybrid", "Terminal"),
               # finPhen = ifelse(finPhen == "H", "Hybrid", "Terminal"))
    phenotypes <- unique(df$initPhen)
    FracRanges <- c(0, 0.25,0.75,1)
    FracLabs <- c("", "Low", "Med", "High")
    df$class <- "Low"
    df$class <- ifelse(df$Fraction >= 0.25 & df$Fraction<0.75, "Med", df$class)
    df$class <- ifelse(df$Fraction >= 0.75, "High", df$class)
    
    lapply(phenotypes, function(x){
        d <- df %>% filter(initPhen == x) %>% group_by(class, init, finPhen) %>% 
            summarise(Frequency = sum(Freq), .groups = "drop") 
        d1 <- d %>% group_by(class, init) %>% 
            summarise(Sum = sum(Frequency), .groups = "drop")
        d <- merge(d,d1,by = c("class", "init"), all = T) %>% mutate(Frequency = Frequency/Sum) 
        
        d %>% 
            group_by(class, finPhen) %>% 
            summarise(Avg = mean(Frequency), Std = sd(Frequency), .groups = "drop") %>%
            mutate(initPhen = x)
        
    }) %>% reduce(rbind.data.frame)
}

aocData <- function(net)
{
    setwd(WTcoherence)
    WT <- read_csv(paste0(net , "_allNodeCoherence_nPert100_nIter10_reps1.csv")) %>% aocCalc %>%
        mutate(Net = net)
    setwd(net)
    filz <- list.files(".", "allNodeCoherence")
    df <- lapply(filz, function(x){
        n <- str_remove(x, "_allNode.*")
        if (n == net)
            return()
        d <- x %>% read_csv(show_col_types = F)
        if(ncol(d)< 3)
            return()
        d %>% aocCalc %>% mutate(Net = n)
    }) %>% reduce(rbind.data.frame)
    df <- rbind.data.frame(WT, df)
    setwd("..")
    write.csv(df, paste0(net, "_aocDat.csv"), row.names = F)
}
sapply(EMPNets, aocData)

aocPlots <- function(net)
{
    setwd(WTcoherence)
    df <- read_csv(paste0(net, "_aocDat.csv"))
    dfGroups <- read_csv(paste0(randcompiled, "/", net, "_rand_groups.csv")) %>% 
        select(Net, Gs)
    Gs <- dfGroups$Gs
    names(Gs) <- dfGroups$Net
    df$Gs <- Gs[df$Net]
    phenotypes <- unique(df$initPhen)
    dat <- lapply(phenotypes, function(x){
        d <- df %>% filter(initPhen == x)
        corD <- lapply(unique(df$class), function(i){
            d1 <- d %>% filter(class == i) 
            multiFactorCorrelation(d1, "finPhen", "Gs", "Avg", label = F) %>% 
                mutate(class = i)
        }) %>% reduce(rbind.data.frame) %>% 
            mutate(significance = ifelse(pValue < 0.05, "", "X"))
        d2 <- d %>% filter(finPhen != "M") %>% select(finPhen, Avg, Gs, class) %>%
            spread(key = finPhen, value = Avg)
        setwd(plotFolder)
        setwd("AUCs")
        d2$class <- factor(d2$class, levels = c("Low", "Med", "High"))
        ggplot(d2, aes(x = E, y = H, color = Gs)) + 
            geom_point(size = 2) + theme_Publication() + 
            scale_color_viridis_c() + 
            theme(legend.position = "right", legend.direction = "vertical", 
                  legend.key.height = unit(0.8, "cm"),
                  axis.text = element_text(size = rel(1.3)),
                  strip.text = element_text(size = rel(1.5))) +
            facet_wrap(~class, nrow = 1) +
            labs(x = "AUC(Epithelial)", y = "AUC(Hybrid)")
        ggsave(paste0(net, "_", x,"_AucScatter.png"), width = 15, height = 5)
        corD$class <- factor(corD$class, levels = c("Low", "Med", "High"))
        ggplot(corD, aes(x = class, y = Factors, fill = Correlation)) + 
            geom_tile() + scale_fill_gradient2(low ="red", high = "blue", limits = c(-1,1)) +
            geom_text(aes(label = significance))+
            theme_Publication() + theme(legend.position = "right", legend.direction = "vertical", 
                                        legend.key.height = unit(0.8, "cm"),
                                        axis.text = element_text(size = rel(1.3))) +
            labs(x = "", y = "")
        ggsave(paste0(net, "_", x, "_heatmap.png"), width = 5.5, height = 5)
        corD %>% mutate(Net = net, Phenotype = x)
    }) %>% reduce(rbind.data.frame)
    return(dat)
}

df <- lapply(EMPNets, aocPlots) %>% reduce(rbind.data.frame)
df$Network <- netNameKey[df$Net]
ggplot(df %>% filter(Network != "22N 82E"), aes(x = class, y = Factors, fill = Correlation)) +
    geom_tile() + scale_fill_gradient2(low ="red", high = "blue", limits = c(-1,1)) +
    geom_text(aes(label = significance))+
    theme_Publication() + theme(legend.position = "right", legend.direction = "vertical", 
                                legend.key.height = unit(0.8, "cm"),
                                axis.text = element_text(size = rel(1.3))) +
    facet_grid(Phenotype ~ Network) +
    labs(x = "", y = "")
setwd(paste0(plotFolder, "/AUCs"))
ggsave("allNetHeat.png", width = 14, height = 9)
### Fig 5H - sigmoidal Fits----

sigmFit <- function(x)
{
    df <- read_csv(x, show_col_types = F)
    if (ncol(df)<3)
        return()
    df <- df %>% select(init, initPhen, nNode, Hamming, Freq) %>% 
        mutate(Fraction = nNode/max(nNode)) %>%
        group_by(init, Fraction) %>% 
        summarise(Hamming = sum(Hamming*Freq), initPhen = unique(initPhen), .groups = "drop")
    inits <- unique(df$init)
    
    fits <- sapply(inits, function(x){
        d <- df %>% filter(init == x)
        xDat <- log(0.5/d$Fraction)
        yDat <- log((1+0.002)/(d$Hamming + 0.001) - 1)
        fit <- lm(yDat~xDat)
        c(fit$coefficients[2], fit$coefficients[1], unique(d$initPhen))
    }) %>% t %>% data.frame %>% set_names(c("Cooperativity", "Intercept", "Phenotype")) %>%
        mutate(Net = x %>% str_remove("_allNodeCoherence.*"), State = inits)
    return(fits)
}
sigmDat <- function(net)
{
    setwd(WTcoherence)
    fitsWT <- paste0(net, "_allNodeCoherence_nPert100_nIter10_reps1.csv") %>% 
        sigmFit %>% mutate(Net = net)
    setwd(net)
    filz <- list.files(".", "allNodeCoherence")
    plan(multisession, workers = 8)
    fittedDat <- future_lapply(filz, sigmFit) %>% reduce(rbind.data.frame)
    future:::ClusterRegistry("stop")
    setwd("..")
    df <- rbind.data.frame(fitsWT, fittedDat)
    write_csv(df, paste0(net, "_sigmoidalFits.csv"))
    print(net)
}
sapply(EMPNets, sigmDat)
sigmPlots <- function(net)
{
    setwd(WTcoherence)
    fitDat <- read_csv(paste0(net, "_sigmoidalFits.csv")) %>%
        group_by(Net, Phenotype) %>% 
        summarise(Avg = mean(Cooperativity), Std = sd(Cooperativity))
    
    groupDat <- read_csv(paste0(randcompiled, "/", net, "_rand_groups.csv")) %>% 
        select(Net, Gs)
    df <- merge(fitDat, groupDat, by = "Net")
    WTDat <- df %>% filter(!str_detect(Net, "rand")) %>% mutate(Net = "WT")
    phenKey <- c("Terminal", "Hybrid", "Terminal")
    names(phenKey) <- c("E", "H", "M")
    df$Phenotype <- phenKey[df$Phenotype]
    WTDat$Phenotype <- phenKey[WTDat$Phenotype]
    ggplot(df, aes(x = Avg)) + 
        geom_histogram(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + 
        theme_Publication() + facet_wrap(~Phenotype, nrow = 1) +
        geom_vline(data = WTDat, aes(xintercept = Avg), color = "red", size = 1.5) +
        labs(x = "Cooperativity", y = "Frequency")
    setwd(plotFolder)
    if(!dir.exists("Sigmoidals"))
        dir.create("Sigmoidals")
    setwd("Sigmoidals")
    ggsave(paste0(net, "_cooperativityDist.png"), width = 9, height = 5)
    corLabs <- multiFactorCorrelation(df, "Phenotype", "Gs", "Avg", method = "spearman")
    names(corLabs) <- unique(df$Phenotype)
    df$grob <- corLabs[df$Phenotype]
    annDat <- df %>% mutate(top = Avg+Std/2) %>% 
        group_by(Phenotype, grob) %>% 
        summarise(x = min(Gs) + 0.5*(max(Gs)-min(Gs)), 
                  y = min(top, na.rm = T) + 1.1*(max(top, na.rm = T)-min(top, na.rm = T)))
                  # y = Inf)
    ggplot(df, aes(x = Gs, y = Avg)) + geom_point() +
        geom_errorbar(aes(ymin = Avg-Std, ymax = Avg+Std)) +
        facet_wrap(~Phenotype) + 
        geom_text(data = annDat, aes(x = x, y = y, label = grob), size = 5, face = "bold") +
        geom_label_repel(data = WTDat, aes(label = Net), color = "red") +
        theme_Publication() + theme(legend.position = "top", 
                                    axis.text = element_text(size = rel(1.3)),
                                    strip.text = element_text(size = rel(1.5))) + 
        labs(x = "Mean Group Strength", y = "Cooperativity", color = "")
    ggsave(paste0(net, "_cooperativityVGs.png"), width = 9, height = 5)
    multiFactorCorrelation(df, "Phenotype", "Gs", "Avg", 
                                     method = "spearman", label = F) %>%
        mutate(Network = netNameKey[net], Sig = ifelse(pValue < 0.05, "", "X"))
    
}
dat <- lapply(EMPNets, sigmPlots) %>% 
    reduce(rbind.data.frame)
df <- dat %>% mutate(Factors = str_extract(Factors, ".")) %>% 
    mutate(Factors = ifelse(Factors == "H", "Hybrid", "Terminal"))
ggplot(df, aes(y = Network, x = Factors, fill = Correlation)) +
    geom_tile() + geom_text(aes(label = Sig)) +
    theme_Publication() + labs(x = "", y = "") + 
    theme(
        # axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.position = "right", 
          legend.direction = "vertical", legend.key.height = unit(0.8, "cm")) +
    scale_fill_gradient2(low = "red", high = "blue", limits = c(-1,1))
ggsave("SigmoidalCorrel.png", width = 5.5, height = 5)

### Fig 5H schematic -----

meanPlotOverlapFits <- function(net, df = NULL)
{
    print(net)
    setwd(WTcoherence)
    if (is.null(df))
    {
        l <- list.files(".", paste0(net, "_allNodeCoherence"))
        df <- read_csv(l, show_col_types = F)
    }
    phenotype <- unique(df$initPhen)
    df <- df %>% mutate(Fraction = nNode/max(nNode)) %>%
        mutate(EMTScore = score_states(fin, paste0(net, ".topo")))
    inits <- c(df %>% group_by(init, initPhen) %>% summarise(Count = n()) %>% 
                   split(.[["initPhen"]]) %>% sapply(function(x){
                       x$init[1]
                   }))
    phens <- df %>% select(init, initPhen) %>% unique %>% filter(init %in% inits)
    
    
    d <- df %>% filter(init %in% inits) %>% group_by(init, Fraction) %>% 
        summarise(Hamming = sum(Hamming*Freq), 
                  Std = sum(Freq*(Hamming - sum(Hamming*Freq))^2), 
                  EMTScore = sum(EMTScore*Freq), 
                  StdScore = sum(Freq*(EMTScore - sum(EMTScore*Freq))^2), 
                  initPhen =  phenKey[unique(initPhen)], .groups = "drop")
    inits <- unique(d$init)
    
    fits <- lapply(inits, function(x){
        d1 <- d %>% filter(init == x)
        xDat <- log(0.5/d1$Fraction)
        yDat <- log((1+0.002)/(d1$Hamming + 0.001) - 1)
        fit <- lm(yDat~xDat)
        n <- fit$coefficients[2]
        d1 %>% select(Fraction) %>% mutate(HammingFit = Fraction^n/((0.5^n) + Fraction^n), 
                                          init = x, Cooperativity = round(n, 2))
    }) %>% reduce(rbind.data.frame)
    
    d1 <- merge(d, fits, by = c("init", "Fraction"), all = T) %>% 
        filter(initPhen != "Mesenchymal") %>% 
        mutate(initPhen = ifelse(initPhen == "Hybrid", "Hybrid", "Terminal")) %>%
        mutate(Phenotype = paste0(initPhen, " : ", Cooperativity))
    ggplot(d1, aes(x = Fraction)) + geom_point(aes(y = Hamming, shape = Phenotype)) +
        geom_line(aes(y = HammingFit, shape = Phenotype)) +
        theme_Publication() + theme(legend.position = c(0.7, 0.2),
                                    legend.direction = "vertical",
                                    axis.text = element_text(size = rel(1.3))) +
        labs(shape = "Phenotype : Cooperativity")

    setwd(plotFolder)
    ggsave(paste0(net, "_fitDemos.png"), width = 5.5, height = 5)
}
sapply(EMPNets, meanPlotOverlapFits)


### Unused plots ----

meanPlotOverlapPhen <- function(net, df = NULL)
{
    if (is.null(df))
    {
        l <- list.files(".", paste0(net, "_allNodeCoherence"))
        df <- read_csv(l, show_col_types = F)
    }
    phenotype <- unique(df$initPhen)
    df <- df %>% mutate(Fraction = nNode/max(nNode)) %>%
        mutate(EMTScore = score_states(fin, paste0(net, ".topo")))
    inits <- c(df %>% group_by(init, initPhen) %>% summarise(Count = n()) %>% 
                   split(.[["initPhen"]]) %>% sapply(function(x){
                       x$init[1]
                   }))
    
    
    d <- df %>% filter(init %in% inits) %>% group_by(init, Fraction) %>% 
        summarise(Hamming = sum(Hamming*Freq), 
                  Std = sum(Freq*(Hamming - sum(Hamming*Freq))^2), 
                  EMTScore = sum(EMTScore*Freq), 
                  StdScore = sum(Freq*(EMTScore - sum(EMTScore*Freq))^2), 
                  initPhen =  phenKey[unique(initPhen)])
    ggplot(d, aes(x = Fraction, y = Hamming, color = EMTScore, shape = initPhen)) + 
        geom_point(size = 3) + #facet_wrap(~initPhen, ncol = 1) +
        # geom_errorbar(aes(ymin = Hamming - Std, ymax = Hamming + Std)) + 
        theme_Publication() + 
        scale_color_viridis_c(limits = c(-1,1))+
        xlim(0,1) +
        theme(legend.position = "right", legend.direction = "vertical",
              legend.key.height = unit(0.8, "cm"), axis.text = element_text(size = rel(1.4))) +
        labs(x = "Perturbation", y = "Hamming")
    setwd(paste0(cwd, "/", "hammings/"))
    if(!dir.exists(net))
        dir.create(net)
    # browser()
    ggsave(paste0(cwd, "/", "hammings/", net, "/",net, "_all.png"), 
           width = 6.5, height = 5)
    setwd(WTcoherence)
}
sapply(EMPNets, meanPlotOverlapPhen)
### rand nets

halfMaxPlots <- function(net)
{
    setwd(randcompiled)
    print(net)
    df <- read_csv(paste0(net, "_phenotypeCoherence.csv"))
    Phenotypes <- df$Phenotype %>% unique
    lapply(Phenotypes, function(x){
        d <- df %>% filter(Phenotype == x) %>% drop_na()
        metrics <- c("IC50", "StateCollapse", "PhenotypeCollapse", 
                     "Slope1", "Slope2", "Slope3")
        metricNames <- c("IC50", "Perturbation for State collapse",
                         "Perturbation for Phenotype Collapse", 
                         "Slope Hamming vs Perturbation", 
                         "Slope Hamming vs Perturbation", 
                         "Slope Hamming vs Perturbation")
        names(metricNames) <- metrics
        sapply(c("G11", "G12", "G22", "G21", "Gs"), function(g)
        {
            grobs <- lapply(metrics, function(m){
                correlGrob(d, g, paste0(m, "_mean"), method = "spearman")
            })
            names(grobs) <- metrics
            sapply(metrics[4:6], function(m){
                ggplot(d, aes_string(x = g, y = paste0(m, "_mean")))+
                    geom_point(size = 2) + 
                    geom_errorbar(aes_string(ymin = paste0(m, "_mean - ", m,"_sd"),
                                             ymax = paste0(m, "_mean + ", m,"_sd"))) +
                    annotation_custom(grobs[[m]]) + theme_Publication() + 
                    labs(x = "Mean Group Strength", y = metricNames[m])
                ggsave(paste0(plotFolder, "/randHammings/", 
                              net, "_", x, "_",g,"_", m, ".png"), width = 5.5, height = 5)
            })
        })
        grobs <- lapply(metrics, function(m){
            correlGrob(d, "Gs", paste0(m, "_mean"), method = "spearman")
        })
        names(grobs) <- metrics
        sapply(metrics, function(m){
            ggplot(d, aes_string(x = "Gs", y = paste0(m, "_mean")))+
                geom_point(size = 2) + 
                geom_errorbar(aes_string(ymin = paste0(m, "_mean - ", m,"_sd"),
                                         ymax = paste0(m, "_mean + ", m,"_sd"))) +
                annotation_custom(grobs[[m]]) + theme_Publication() + 
                labs(x = "Mean Group Strength", y = metricNames[m])
            ggsave(paste0(cwd, "/randHammings/", net, "_", x, "_", m, ".png"), width = 5.5, height = 5)
        })
        d <- d %>% select(all_of(paste0(c("Slope1", "Slope2", "Slope3"), "_mean")), Gs) %>%
            gather(key = "Metric", value = "Value", -Gs) 
        # browser()
        corrs <- multiFactorCorrelationAnova(d, "Metric", "Value", "Gs", 
                                             label = F, method = "spearman") %>% 
            mutate(Phenotype = x, Net = net)
    }) %>% reduce(rbind.data.frame)
    
    
}
df <- lapply(EMPNets, halfMaxPlots)
df <- df %>% reduce(rbind.data.frame)
df$Phenotype <- phenKey[df$Phenotype]
phen <- unique(df$Phenotype)
df$Metric <- df$Factors %>% str_remove("_mean")
df$Net <- netNameKey[df$Net]
df$Label <- ifelse(df$pValue < 0.05, "", "X")
sapply(phen, function(x){
    d <- df %>% filter(Phenotype == x)
    ggplot(d, aes(x = Net, y = Metric, fill = Correlation)) + geom_tile() +
        scale_fill_gradient2(low = "blue", high = "red", limits = c(-1,1)) +
        theme_Publication() + 
        theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90),
              legend.position = "right", legend.direction = "vertical", 
              legend.key.height = unit(0.8, "cm")) +
        facet_wrap(~Phenotype, ncol = 1) +
        labs(x = "", y = "")
    ggsave(paste0(cwd, "/randHammings/", x, "_heatmap.png"), width = 7, height = 4)
})

ggplot(df, aes(x = Net, y = Metric, fill = Correlation)) + geom_tile() +
    geom_text(aes(label = Label)) +
    scale_fill_gradient2(low = "blue", high = "red", limits = c(-1,1)) +
    theme_Publication() + 
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "right", legend.direction = "vertical", 
          legend.key.height = unit(0.8, "cm")) +
    facet_wrap(~Phenotype, nrow = 1) +
    labs(x = "", y = "")
ggsave(paste0(cwd, "/randHammings/all_heatmap.png"), width = 10, height = 4)

