mainFolder <- "D:/Teams"
plotFolder <- paste0(mainFolder, "/Figures/Fig3")
source(paste0(mainFolder, "/Figures/figCore.R"))
source(paste0(mainFolder, "/Final_Results/codes/stateLabeller.R"))

### Fig 3B - coherence bimodality ----

WTcohViolins <- function(netList)
{
    dfs <- lapply(netList, function(x){
        df <- read.csv(paste0(randRaw, "\\", x, "\\", x, "_finFlagFreq.csv"))
        df <- df %>% filter(flag == 1) %>% select(coherence0) %>% set_names("Coherence") %>%
            drop_na() %>% mutate(Network = netNameKey[x] %>% str_replace(" ", "\n"))
        df
    }) %>% reduce(rbind.data.frame) %>% mutate(Network = str_replace_all(Network, "_", "\n"))
    ggplot(dfs, aes(x = Network, y = Coherence)) + geom_violin() +
        theme_Publication() + theme(axis.title.x = element_blank())
    setwd(plotFolder)
    
    if(!dir.exists("violins"))
    {
        dir.create("violins")
    }
    setwd("violins")
    ggsave(paste0(paste0(netList, collapse = "_"), "_coherencePlot.png"), width = 6.5, height = 5)
    
}
WTcohViolins(EMPNets)

### Fig 3C - min, max and mean coherence ----

cohViolinNpPlotter <- function(net, plot = F)
{
    # browser()
    print(net)
    setwd(randcompiled)
    dfWT <- read.csv(paste0(net, "_WT_Metrics.csv")) %>%
        select(Coherence) %>% set_names("Coherence") %>% drop_na() %>% 
        mutate(Class = "WT") %>%
        select(Class, Coherence)
    wtDat <- c(min(dfWT$Coherence), max(dfWT$Coherence), mean(dfWT$Coherence))
    dfrand1 <- read.csv(paste0(net, "_ALL.csv")) %>% select(Net, minCoh, maxCoh, meanCoh) %>%
        drop_na()
    # wtDat <- dfrand1 %>% filter(!str_detect(Net, "rand")) %>% select(-Net) %>% unlist
    dfrand <- dfrand1 %>%
        select(minCoh, maxCoh) %>% drop_na() %>% 
        gather(key = "Class", value = "Coherence") %>% 
        mutate(Class = paste0("Rand\n", labelvals[Class] %>% str_replace(" ", "\n")))
    # dfrand <- dfrand1 %>%
    #     gather(key = "Class", value = "Coherence") %>% 
    #     mutate(Class = paste0("Rand\n", labelvals[Class] %>% str_replace(" ", "\n")))
    setwd(plotFolder)
    df <- rbind.data.frame(dfWT, dfrand) %>% 
        mutate(Class = factor(Class, levels = c("Rand\nMinimum\nCoherence", 
                                                "WT", "Rand\nMaximum\nCoherence")))
    if(!dir.exists("violins"))
    {
        dir.create("violins")
    }
    setwd("violins")
    if (plot)
    {
        ggplot(df, aes(x = Class, y = Coherence)) + geom_violin() + 
            theme_Publication() + 
            theme(axis.title.x = element_blank(), 
                  axis.text.x = element_text(size = rel(1)))
        ggsave(paste0(net, "_cohViolin.png"), width = 6.5, height = 5)
    }
    
    res <- c(sum(dfrand1$minCoh > wtDat[1]), 
             sum(dfrand1$maxCoh > wtDat[2]), 
             sum(dfrand1$meanCoh > wtDat[3]))/nrow(dfrand1)
    return(res)
}
dat <- sapply(EMPNets, cohViolinNpPlotter, plot = T) %>% data.frame %>% 
    mutate(Metric = c("Min\n Coherence", "Max\n Coherence", "Mean\n Coherence")) %>% 
    gather(key = "Network", value = "Fraction", -Metric) %>%
    mutate(Network = netNameKey[Network])
setwd(plotFolder)
ggplot(dat, aes(x = Metric, y = Network, fill = Fraction)) + geom_tile() +
    scale_x_discrete(expand = c(0,0)) +
    theme_Publication() + scale_fill_viridis_c() + 
    theme(legend.key.height = unit(0.8, "cm"), 
          legend.position = "right", legend.direction = "vertical",
          axis.text = element_text(size = rel(1)))
ggsave("cohpValueMatrix.png")


### Fig 3D - Coherence vs Gs ----

GsCohplots <- function(net)
{
    setwd(randcompiled)
    df <- paste0(net, "_ALL.csv") %>% read_csv %>% filter(str_detect(Net, "rand")) %>%
        select(minCoh, maxCoh, meanCoh, Net)
    key <- c("Min\nCoherence", "Max\nCoherence", "Mean\nCoherence")
    names(key) <- c("minCoh", "maxCoh", "meanCoh")
    dfGroup <- paste0(net, "_rand_groups.csv") %>% read_csv %>% 
        mutate(Gs = (abs(G11) + abs(G22) + abs(G12) + abs(G21))/4) %>%
        select(Net, Gs) %>% filter(str_detect(Net, "rand"))
    df <- merge(df, dfGroup, by = "Net", all = T) 
    WT <- df %>% filter(!str_detect(Net, "rand"))
    grobLow <- correlGrob(df, "minCoh", "Gs")
    grobHigh <- correlGrob(df, "maxCoh", "Gs")
    
    lab <- "WT"
    setwd(plotFolder)
    if (!dir.exists("CohVsGsCorrels"))
        dir.create("CohVsGsCorrels")
    setwd("CohVsGsCorrels")    
    ggplot(df, aes(x = Gs, y = minCoh)) + 
        geom_point() + 
        geom_label_repel(data = WT, 
                         color = "red", label = lab) + 
        annotation_custom(grobLow)+
        theme_Publication()+
        labs(x = "Mean Group Strength", y = "Minimum Coherence")
    ggsave(paste0(net, "_lowCoh.png"), width = 5.5, height = 5)
    ggplot(df, aes(x = Gs, y = maxCoh)) + 
        geom_point() + 
        geom_label_repel(data = WT,
                         color = "red", label = lab) + 
        annotation_custom(grobHigh)+
        theme_Publication()+
        labs(x = "Mean Group Strength", y = "Maximum Coherence")
    ggsave(paste0(net, "_highCoh.png"), width = 5.5, height = 5)
    
    df <- df %>% 
        gather(key = "Metric", value = "Coherence", -Gs, -Net) %>% 
        mutate(Metric = key[Metric])
    multiFactorCorrelation(df, "Metric", "Coherence", "Gs", label = F) %>%
        mutate(Network = netNameKey[net])
    
}

df <- lapply(EMPNets, GsCohplots) %>% 
    reduce(rbind.data.frame) %>% mutate(Label = ifelse(pValue < 0.05, "", "X"))

ggplot(df, aes(x = Factors, y = Network, fill = Correlation)) + 
    geom_tile() + geom_text(aes(label = Label)) +
    theme_Publication() + 
    scale_x_discrete(expand = c(0,0)) +
    theme(axis.text = element_text(size = rel(1)),
        legend.position = "right", legend.direction = "vertical", 
        legend.key.height = unit(0.8, "cm")) +
    scale_fill_gradient2(low = "blue", high = "red", limits = c(-1,1))+
    labs(x = "", y = "")
ggsave(paste0(plotFolder, "/", net, "_cohGsCorrel.png"))

### Fig 3E - coherence vs frustration and frequency ----

cohFreqAll <- function(netList, compdir = randcompiled, plotDir = cwdm, key)
{
    setwd(compdir)
    cohDf <- list.files(".", "cohMnStd") %>% lapply(function(x){
        dfCoh <- read.csv(x)
        net <- x %>% str_remove("cohMnStd.csv")
        i <- which(netList == net)
        dfCoh$Network <- key[i]
        cr <- cor.test(dfCoh$Freq, dfCoh$Mean, method = "spearman")
        dfCoh$Network1 <- paste0(dfCoh$Network, ", \u03c1 : ", round(cr$estimate,2), 
                                 ifelse(cr$p.value < 0.05, "*", ""))
        cr <- cor.test(dfCoh$Frust, dfCoh$Mean, method = "spearman")
        dfCoh$Network2 <- paste0(dfCoh$Network, ", \u03c1 : ", round(cr$estimate,2), 
                                 ifelse(cr$p.value < 0.05, "*", ""))
        
        dfCoh
    }) %>% reduce(rbind.data.frame)
    ggplot(cohDf, aes(x = Freq, y = Mean, color = Network1)) + geom_point(size = 2) + 
        geom_errorbar(aes(ymin = Mean - Stdev, ymax = Mean + Stdev)) +
        # geom_smooth(se = F, method = "lm")+
        theme_Publication() + theme(legend.position = c(0.7,0.6), 
                                    legend.direction = "vertical", 
                                    legend.key.height = unit(0.6, "cm"),
                                    legend.text = element_text(size = rel(1.1)),
                                    legend.title = element_text(size = rel(1.1))) +
        labs(x = "Frequency", y = "Coherence", color = "Network")
    setwd(plotFolder)
    ggsave("allNetsCohVfreq.png", width = 5.5, height = 5)
    ggplot(cohDf, aes(x = Frust, y = Mean, color = Network2)) + geom_point(size = 2) + 
        geom_errorbar(aes(ymin = Mean - Stdev, ymax = Mean + Stdev)) +
        # geom_smooth(se = F, method = "lm")+
        theme_Publication() + theme(legend.position = c(0.3,0.3), 
                                    legend.direction = "vertical", 
                                    legend.key.height = unit(0.6, "cm"),
                                    legend.text = element_text(size = rel(1.1)),
                                    legend.title = element_text(size = rel(1.1))) +
        labs(x = "Frustration", y = "Coherence", color = "Network")
    setwd(plotFolder)
    ggsave("allNetsCohVfrust.png", width = 5.5, height = 5)
}
cohFreqAll(netList = EMPNets, key = netNameKey)
