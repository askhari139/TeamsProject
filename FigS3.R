mainFolder <- "D:/Teams"
plotFolder <- paste0(mainFolder, "/Figures/FigS3")
source(paste0(mainFolder, "/Figures/figCore.R"))
source(paste0(mainFolder, "/Final_Results/codes/stateLabeller.R"))

### Bimodality coefficients ----

dfAll <- lapply(EMPNets, function(net){
    df_metrics <- read.csv(paste0(randcompiled, "/", net, "_ALL.csv"))
    dfGroups <- read.csv(paste0(randcompiled, "/", net, "_rand_groups.csv"))
    df <- merge(dfGroups, df_metrics, by = "Net", all = T) %>%
        mutate(Network = netNameKey[net])
    dfWT <- df %>% filter(!str_detect(Net, "_rand"))
    setwd(plotFolder)
    if (!dir.exists("Bimodality"))
        dir.create("Bimodality")
    setwd("Bimodality")
    ggplot(df, aes(x = bmSSF)) + geom_histogram(aes(y = ..count../sum(..count..))) +
        geom_vline(data = dfWT, aes(xintercept = bmSSF), color = "red", size = 1.5) +
        theme_Publication() + facet_wrap(~Network) +
        labs(x = "SSF Bimodality Coefficient", y = "Frequency")
    ggsave(paste0(net, "_bmSSFDist.png"), width = 5.5, height = 5)
    grob <- correlGrob(df, "bmSSF", "Gs", xPos = 0.7, yPos = 0.1)
    ggplot(df, aes(x = Gs, y = bmSSF)) + geom_point() +
        annotation_custom(grob)+
        geom_label_repel(data = dfWT, label = "WT", color = "red") +
        theme_Publication() + facet_wrap(~Network) +
        labs(y = "SSF Bimodality Coefficient", x = "Mean Group Strength")
    ggsave(paste0(net, "_bmSSFGs.png"), width = 5.5, height = 5)
    
    ggplot(df, aes(x = bmFrustration)) + geom_histogram(aes(y = ..count../sum(..count..))) +
        geom_vline(data = dfWT, aes(xintercept = bmFrustration), color = "red", size = 1.5) +
        theme_Publication() + facet_wrap(~Network) +
        labs(x = "Frustraion Bimodality Coefficient", y = "Frequency")
    ggsave(paste0(net, "_bmFrustDist.png"), width = 5.5, height = 5)
    grob <- correlGrob(df, "bmFrustration", "Gs", xPos = 0.7, yPos = 0.1)
    ggplot(df, aes(x = Gs, y = bmFrustration)) + geom_point() +
        annotation_custom(grob)+
        geom_label_repel(data = dfWT, label = "WT", color = "red") +
        theme_Publication() + facet_wrap(~Network) +
        labs(y = "Frustration Bimodality Coefficient", x = "Mean Group Strength")
    ggsave(paste0(net, "_bmFrustGs.png"), width = 5.5, height = 5)
    
    ggplot(df, aes(x = corFreqFrust)) + geom_histogram(aes(y = ..count../sum(..count..))) +
        geom_vline(data = dfWT, aes(xintercept = corFreqFrust), color = "red", size = 1.5) +
        theme_Publication() + facet_wrap(~Network) +
        labs(x = "\u03c1(Frustration, SSF)", y = "Frequency")
    ggsave(paste0(net, "_corFreqFrust.png"), width = 5.5, height = 5)
    grob <- correlGrob(df, "corFreqFrust", "Gs")
    ggplot(df, aes(x = Gs, y = corFreqFrust)) + geom_point() +
        annotation_custom(grob)+
        geom_label_repel(data = dfWT, label = "WT", color = "red") +
        theme_Publication() + facet_wrap(~Network) +
        labs(y = "\u03c1(Frustration, SSF)", x = "Mean Group Strength")
    ggsave(paste0(net, "_corFreqFrustGs.png"), width = 5.5, height = 5)
    
    
    df
}) %>% reduce(rbind.data.frame)


df <- dfAll
dfWT <- df %>% filter(!str_detect(Net, "rand"))

ggplot(df, aes(x = bmSSF)) + geom_histogram(aes(y = ..count../sum(..count..))) +
    geom_vline(data = dfWT, aes(xintercept = bmSSF), color = "red", size = 1.5) +
    theme_Publication() + facet_wrap(~Network, nrow = 1) +
    labs(x = "SSF Bimodality Coefficient", y = "Frequency")
ggsave(paste0("all", "_bmSSFDist.png"), width = 14, height = 5)


ggplot(df, aes(x = bmFrustration)) + geom_histogram(aes(y = ..count../sum(..count..))) +
    geom_vline(data = dfWT, aes(xintercept = bmSSF), color = "red", size = 1.5) +
    theme_Publication() + facet_wrap(~Network, nrow = 1) +
    labs(x = "Frustraion Bimodality Coefficient", y = "Frequency")
ggsave(paste0("all", "_bmFrustDist.png"), width = 14, height = 5)

ggplot(df, aes(x = corFreqFrust)) + geom_histogram(aes(y = ..count../sum(..count..))) +
    geom_vline(data = dfWT, aes(xintercept = corFreqFrust), color = "red", size = 1.5) +
    theme_Publication() + facet_wrap(~Network, nrow = 1) +
    labs(x = "\u03c1(Frustration, SSF)", y = "Frequency")
ggsave(paste0("all", "_corFreqFrust.png"), width = 14, height = 5)


### Fraction Heatmap----

dfFrac <- lapply(EMPNets, function(net){
    d <- df %>% filter(Network == net)
    
})
