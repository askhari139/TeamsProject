mainFolder <- "D:/Teams"
plotFolder <- paste0(mainFolder, "/Figures/Fig2")
source(paste0(mainFolder, "/Figures/figCore.R"))
source(paste0(mainFolder, "/Final_Results/codes/stateLabeller.R"))

### Fig 2A - Steady state Frequency distributions ----
setwd(randRaw)
df <- lapply(EMPNets, function(x){
    d <- read.csv(paste0(x, "/", x, "_finFlagFreq.csv")) %>%
        filter(flag == 1)
    nedge <- read.delim(paste0(x, "/", x, ".topo"), sep = " ") %>% nrow
    n <- d$states[1] %>% str_split("") %>% unlist %>% length - 2
    d %>% select(Avg0) %>% drop_na() %>% mutate(Network = paste0(n, "N ", nedge, "E"))
}) %>% reduce(rbind.data.frame)

ggplot(df, aes(x = Avg0, fill = Network,y = ..count../sum(..count..))) + 
    geom_histogram(alpha = 0.5, color = "black") + 
    # geom_density() + 
    scale_x_log10() +
    theme_Publication() + 
    labs(x = "Frequency", y = "Frequency") +
    theme(legend.position = "top")
setwd(plotFolder)
ggsave("SteadystateFreq.png", width = 7, height = 5)

### Fig 2C - Frustration distributions ----

WTfrustViolins <- function(netList)
{
    dfs <- lapply(netList, function(x){
        df <- read.csv(paste0(randRaw, "\\", x, "\\", x, "_finFlagFreq.csv"))
        df <- df %>% filter(flag == 1) %>% select(frust0) %>% set_names("Frustration") %>%
            drop_na() %>% mutate(Network = netNameKey[x] %>% str_replace(" ", "\n"))
        df
    }) %>% reduce(rbind.data.frame) %>% mutate(Network = str_replace_all(Network, "_", "\n"))
    ggplot(dfs, aes(x = Network, y = Frustration)) + geom_violin() +
        theme_Publication() + theme(axis.title.x = element_blank())
    setwd(plotFolder)
    
    if(!dir.exists("violins"))
    {
        dir.create("violins")
    }
    setwd("violins")
    ggsave(paste0(paste0(netList, collapse = "_"), "_frustrationPlot.png"), width = 5.5, height = 5)
    
}
WTfrustViolins(EMPNets)


### Fig 2D - WT vs rand frustration ----

frustViolinNpPlotter <- function(net)
{
    # browser()
    print(net)
    setwd(paste0(randRaw, "/", net))
    dfWT <- read.csv(paste0(net, "_finFlagFreq.csv")) %>% filter(flag == 1) %>% 
        select(frust0) %>% set_names("Frustration") %>% drop_na() %>% 
        mutate(Class = "WT") %>%
        select(Class, Frustration)
    wtDat <- c(min(dfWT$Frustration), max(dfWT$Frustration), mean(dfWT$Frustration))
    setwd(randcompiled)
    dfrand1 <- read.csv(paste0(net, "_ALL.csv")) %>% select(minFrust, maxFrust, meanFrust) %>%
        drop_na()
    dfrand <- dfrand1 %>%
        gather(key = "Class", value = "Frustration") %>% 
        mutate(Class = paste0("Rand\n", labelvals[Class] %>% str_replace(" ", "\n")))
    setwd(cwd)
    df <- rbind.data.frame(dfWT, dfrand) %>% 
        filter(!str_detect(Class, "Mean")) %>% 
        mutate(Class = factor(Class, levels = c("Rand\nMinimum\nFrustration", 
                                                "WT", "Rand\nMaximum\nFrustration")))
    
    ggplot(df, aes(x = Class, y = Frustration)) + geom_violin() + 
        theme_Publication() + 
        theme(axis.title.x = element_blank(), 
              axis.text.x = element_text(size = rel(1)))
    setwd(plotFolder)
    
    if(!dir.exists("violins"))
    {
        dir.create("violins")
    }
    setwd("violins")
    
    ggsave(paste0(net, "_frustViolin.png"), width = 7.5, height = 5)
}

### Fg 2E ----

### Fig 2F - Frequency vs Frustration correlation ----

frustFreqplotAll <- function(netList, plotDir = cwd, key)
{
    # print(dir)
    setwd(randRaw)
    # net <- str_remove(topoFile, ".topo")
    freqdf <- lapply(netList, function(net){
        i <- which(netList == net)
        d <- read.csv(paste0(net, "/", net, "_finFlagFreq.csv")) %>%
            filter(flag == 1) %>% select(Avg0, SD0, frust0) %>%
            drop_na() %>% mutate(Network = key[i])
        cr <- cor.test(d$Avg0, d$frust0, method = "spearman")
        d$Network <- paste0(d$Network, ", \u03c1 : ", round(cr$estimate,2), ifelse(cr$p.value < 0.05, "*", ""))
        d
    }) %>% reduce(rbind.data.frame)
    setwd(plotDir)
    # corr <- cor(freqdf$Avg0, freqdf$frust0)
    # grob <- grobTree(textGrob(paste0("\u03c1 : ", round(corr, 2)), x=0.5,  y=0.9, hjust=0,
    #                           gp=gpar(col="black", fontsize=18, fontface="bold")))
    ggplot(freqdf, aes(x = Avg0, y = frust0, color = Network)) + geom_point(size = 2) + 
        geom_errorbarh(aes(xmin = Avg0-SD0, xmax = Avg0+SD0)) +
        # geom_smooth(method = "lm") + 
        # annotation_custom(grob) + 
        theme_Publication() +
        theme(legend.position = c(0.7,0.6), 
              legend.direction = "vertical", 
              legend.key.height = unit(0.6, "cm"),
              legend.text = element_text(size = rel(1.1)),
              legend.title = element_text(size = rel(1.1))) +
        labs(x = "Frequency", y = "Frustration")
    setwd(plotFolder)
    ggsave(paste0("allNet_freqVfrustScatter.png"), width = 5.5, height = 5)
}
frustFreqplotAll(EMPNets, key = netNameKey)

### Fig 2G, S3D - RACIPE frust vs freq ----

frustPlot <- function(net)
{
    setwd(cwd)
    ssDat <- read_csv(paste0(RACIPE_WT, "/", net, "/", net,"_discreteStates.csv"))
    grob <- correlGrob(ssDat, "Frequency", "Frustration", method = "spearman")
    ggplot(ssDat, aes(x = Frequency, y = Frustration)) + geom_point(size = 2) + 
        # geom_smooth(method = "lm") + 
        annotation_custom(grob) + theme_Publication() +
        labs(x = "Frequency", y = "Frustration")
    setwd(plotFolder)
    dir.create("RACIPE_frustFreq")
    setwd("RACIPE_frustFreq")
    ggsave(paste0(net, "_freqVfrustScatter.png"), width = 5.5, height = 5)
}

sapply(EMPNets, frustPlot)


### Fig 2H, I - Gs vs frustration correlation ----

plotGsFrustFreq <- function(net)
{
    setwd(paste0(randRaw, "/", net))
    frustDat <- read_csv(paste0(randcompiled, "/", net, "_ALL.csv")) %>% 
        select(minFrust, maxFrust, Net) %>% set_names(c("LowFrust", "HighFrust", "Net"))
    groups <- read_csv(paste0(randcompiled, "/", net, "_rand_groups.csv")) %>%
        mutate(Gs = (abs(G11) + abs(G12) + abs(G21) + abs(G22))/4)
    
    df <- merge(groups %>% select(Net, Gs), frustDat, by = "Net", all = T)
    WT <- df %>% filter(Net == net)
    grobLow <- correlGrob(df, "LowFrust", "Gs")
    grobHigh <- correlGrob(df, "HighFrust", "Gs", yPos = 0.2)
    
    lab <- "WT"
    setwd(plotFolder)
    if (!dir.exists("FrustVGsCorrels"))
        dir.create("FrustVGsCorrels")
    setwd("FrustVGsCorrels")    
    ggplot(df, aes(x = Gs, y = LowFrust)) + 
        geom_point() + 
        geom_label_repel(data = WT, 
                         color = "red", label = lab) + 
        annotation_custom(grobLow)+
        theme_Publication()+
        labs(x = "Mean Group Strength", y = "Minimum Frustration")
    ggsave(paste0(net, "_lowFrust.png"), width = 5.5, height = 5)
    ggplot(df, aes(x = Gs, y = HighFrust)) + 
        geom_point() + 
        geom_label_repel(data = WT,
                         color = "red", label = lab) + 
        annotation_custom(grobHigh)+
        theme_Publication()+
        labs(x = "Mean Group Strength", y = "Maximum Frustration")
    ggsave(paste0(net, "_highFrust.png"), width = 5.5, height = 5)
}
sapply(EMPNets, plotGsFrustFreq)

### Fig 2J - heatmap of frustration vs Gs ----

GsFrustHeatmap <- function(net)
{
    setwd(randcompiled)
    df <- paste0(net, "_ALL.csv") %>% read_csv %>% filter(str_detect(Net, "rand")) %>%
        select(minFrust, maxFrust, meanFrust, Net)
    key <- c("Min\nFrust", "Max\nFrust", "Mean\nFrust")
    names(key) <- c("minFrust", "maxFrust", "meanFrust")
    dfGroup <- paste0(net, "_rand_groups.csv") %>% read_csv %>% 
        mutate(Gs = (abs(G11) + abs(G22) + abs(G12) + abs(G21))/4) %>%
        select(Net, Gs) %>% filter(str_detect(Net, "rand"))
    df <- merge(df, dfGroup, by = "Net", all = T) %>% 
        gather(key = "Metric", value = "Frustration", -Gs, -Net) %>% 
        mutate(Metric = key[Metric])
    
    multiFactorCorrelation(df, "Metric", "Frustration", "Gs", label = F) %>%
        mutate(Network = netNameKey[net])
    
}

df <- lapply(EMPNets, GsFrustHeatmap) %>% 
    reduce(rbind.data.frame) %>% mutate(Label = ifelse(pValue < 0.05, "", "X"))

ggplot(df, aes(x = Factors, y = Network, fill = Correlation)) + 
    geom_tile() + geom_text(aes(label = Label)) +
    theme_Publication() + 
    scale_x_discrete(expand = c(0,0)) +
    theme(#axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "right", legend.direction = "vertical", 
        legend.key.height = unit(0.8, "cm")) +
    scale_fill_gradient2(low = "blue", high = "red", limits = c(-1,1))+
    labs(x = "", y = "")
setwd(plotFolder)
ggsave(paste0("all_frustGsCorrel.png"), width = 5.5, height = 5)





