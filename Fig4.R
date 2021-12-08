mainFolder <- "D:/Teams"
plotFolder <- paste0(mainFolder, "/Figures/Fig4")
source(paste0(mainFolder, "/Figures/figCore.R"))
source(paste0(mainFolder, "/Final_Results/codes/stateLabeller.R"))

### Fig 4A i - 2D plots ----
cohFreqplot <- function(net, dir, plotDir = cwd)
{
    print(dir)
    setwd(dir)
    # net <- str_remove(topoFile, ".topo")
    freqdf <- read.csv(paste0(net, "_finFlagFreq.csv")) %>%
        filter(flag == 1) %>% select(states,Avg0, SD0, frust0) %>%
        drop_na()
    setwd(randcompiled)
    cohDf <- read.csv(paste0(net, "cohMnStd.csv"))
    scoreDf <- read.csv(paste0(net, "_WT_Metrics.csv")) %>% select(States, Score) %>% 
        set_names(c("states", "EMTScore"))
    colnames(cohDf)[1] <- "states"
    freqdf <- merge(freqdf, cohDf, by = "states", all = T)
    freqdf <- merge(freqdf, scoreDf, by = "states", all = T)
    setwd(plotDir)
    
    grob1 <- correlGrob(freqdf, "Avg0", "Mean")
    grob2 <- correlGrob(freqdf, "frust0", "Mean")
    
    ggplot(freqdf, aes(x = Avg0, y = Mean)) + geom_point(size = 2) + 
        geom_errorbar(aes(ymin = Mean-Stdev, ymax = Mean+Stdev)) +
        geom_errorbarh(aes(xmin = Avg0-SD0, xmax = Avg0+SD0)) +
        # geom_smooth(method = "lm") + 
        annotation_custom(grob1) + theme_Publication() +
        labs(x = "Frequency", y = "Coherence")
    ggsave(paste0(net, "_freqVCohScatter.png"), width = 6.5, height = 5)
    ggplot(freqdf, aes(x = frust0, y = Mean, color = Avg0)) + geom_point(size = 2) + 
        geom_errorbar(aes(ymin = Mean-Stdev, ymax = Mean+Stdev), width = 0.01) +
        # geom_errorbarh(aes(xmin = Avg0-SD0, xmax = Avg0+SD0)) +
        # geom_line() + 
        annotation_custom(grob2) + theme_Publication() + scale_color_viridis_c()+
        theme(legend.position = "right", legend.direction = "vertical", 
              legend.key.size = unit(1, "cm")) +
        labs(x = "Frustration", y = "Coherence", color = "Frequency")
    ggsave(paste0(net, "_frustVCohScatter.png"), width = 5.5, height = 5)
    ggplot(freqdf, aes(x = EMTScore, y = Mean, color = Avg0)) + geom_point(size = 2) + 
        geom_errorbar(aes(ymin = Mean-Stdev, ymax = Mean+Stdev), width = 0.01) +
        # geom_errorbarh(aes(xmin = Avg0-SD0, xmax = Avg0+SD0)) +
        # geom_line() + 
        # annotation_custom(grob2) + 
        theme_Publication() + scale_color_viridis_c()+
        theme(legend.position = "right", legend.direction = "vertical", 
              legend.key.size = unit(1, "cm")) +
        labs(x = "EMT Score", y = "Coherence", color = "Frequency")
    ggsave(paste0(net, "_scoreVCohScatter.png"), width = 6.5, height = 5)
}

### Fig 4A ii - Top bottom Gs networks ----

topBottom <- function(net)
{
    netsList <- paste0(randcompiled, "/", net, "_rand_groups.csv") %>% read_csv %>%
        mutate(across(c(G11, G22, G12, G21), as.numeric)) %>%
        mutate(Gs = (abs(G11) + abs(G12) + abs(G22) + abs(G21))/4) %>%
        select(Net, Gs) %>% filter(str_detect(Net, "rand")) %>%
        arrange(Gs) %>% slice(c(1:10, (nrow(.) - 9): (nrow(.)))) %>%
        select(Net) %>% unlist
    class <- c(rep("Low", 10), rep("High", 10))
    names(class) <- netsList
    setwd(paste0(randRaw, "/", net))
    df <- lapply(netsList, function(x){
        print(x)
        # browser()
        d <- read_csv(paste0(x, "_finFlagFreq.csv"), show_col_types = F) %>%
            filter(!is.na(Avg0), flag == 1) 
        if (!any(colnames(d) == "coherence0"))
            return()
        d <- d %>% 
            select(Avg0, frust0, coherence0, phenotype) %>%
            mutate(Phenotype = ifelse(phenotype == "H", "Hybrid", "Terminal")) %>% 
            drop_na() %>% mutate(Gs = class[x])
    }) %>% reduce(rbind.data.frame)
    
    ggplot(df, aes(x = frust0, y = coherence0, color = Gs, shape = Phenotype)) +
        geom_point() + theme_Publication() + 
        theme(legend.position = c(0.4, 0.15), 
              legend.text = element_text(size = rel(1.1)),
              legend.title = element_text(size = rel(1.1)), 
              legend.spacing = unit(0.1, "cm")) +
        labs(x = "Frustration", y = "Coherence", color = "Mean Group Strength")
    setwd(plotFolder)
    if(!dir.exists("topBottom"))
        dir.create("topBottom")
    setwd("topBottom")
    ggsave(paste0( net, "_topBottomRand.png"), width = 5.5, height = 3.5)
    
    ggplot(df, aes(x = frust0 + coherence0, fill = Gs)) + 
        geom_density() + theme_Publication() + 
        theme(legend.position = "none") +
        labs(x = "Frustration + Coherence", y = "")
    ggsave(paste0(net, "_topBottomDist.png"), width = 5.5, height = 2)
    
    ggplot(df, aes(x = frust0, fill = Gs)) + 
        geom_density() + theme_Publication() + 
        theme(legend.position = "none") +
        labs(x = "Frustration", y = "")
    ggsave(paste0(net, "_topBottomFrustDist.png"), width = 5.5, height = 2)
    
    ggplot(df, aes(x = coherence0, fill = Gs)) + 
        geom_density() + theme_Publication() + 
        theme(legend.position = "none") +
        labs(x = "Coherence", y = "")
    ggsave(paste0(net, "_topBottomCohDist.png"), width = 5.5, height = 2)
    
}
sapply(EMPNets, topBottom)


topBottomMeans <- function(net)
{
    netsList <- paste0(randcompiled, "/", net, "_rand_groups.csv") %>% read_csv %>%
        mutate(across(c(G11, G22, G12, G21), as.numeric)) %>%
        mutate(Gs = (abs(G11) + abs(G12) + abs(G22) + abs(G21))/4) %>%
        select(Net, Gs) %>% filter(str_detect(Net, "rand")) %>%
        arrange(Gs) %>% slice(c(1:10, (nrow(.) - 9): (nrow(.)))) %>%
        select(Net) %>% unlist
    class <- c(rep("Low", 10), rep("High", 10))
    names(class) <- netsList
    setwd(paste0(randRaw, "/", net))
    df <- lapply(netsList, function(x){
        read_csv(paste0(x, "_finFlagFreq.csv"), show_col_types = F) %>%
            filter(!is.na(Avg0), flag == 1) %>% 
            select(Avg0, frust0, coherence0) %>%
            drop_na() %>% mutate(Gs = class[x])
    }) %>% reduce(rbind.data.frame)
    
    ggplot(df, aes(x = frust0, y = coherence0, color = Gs)) +
        geom_point() + theme_Publication() + 
        theme(legend.position = c(0.4, 0.1), 
              legend.text = element_text(size = rel(1.1)),
              legend.title = element_text(size = rel(1.1))) +
        labs(x = "Frustration", y = "Coherence", color = "Mean Group Strength")
    ggsave(paste0(cwd, "/", net, "_topBottomMeans.png"), width = 5.5, height = 5)
    
}
sapply(EMPNets, topBottomMeans)

### Fig 4B, S5B - state matrices ----

matrixPlot <- function(net)
{
    setwd(cwd)
    if(!dir.exists("StateMatrix"))
        dir.create("StateMatrix")
    setwd("StateMatrix")
    freqFile <- paste0(randRaw, "/", net, "/", net, "_finFlagFreq.csv")
    nodes <- readLines(paste0(randRaw, "/", net, "/", net, "_nodes.txt")) %>% 
        str_replace_all(regex("\\W+"), "")
    freqData <- read_csv(freqFile, col_types = cols()) %>% filter(flag == 1, !is.na(Avg0)) %>% 
        select(states, phenotype) %>%
        mutate(phenotype = factor(phenotype, levels = rev(c("H", "M", "E")))) %>%
        mutate(states = str_remove_all(states, "'")) %>%
        arrange(phenotype) %>% 
        mutate(states = str_split(states, "") %>% sapply(function(x){paste0(x, collapse = "_")})) %>%
        separate(states, nodes, sep = "_")
    topoFile <- paste0(randRaw, "/", net, "/", net, ".topo")
    file.copy(topoFile, paste0(net, ".topo"))
    topoFile <- paste0(net, ".topo")
    size <- 0.8
    if (length(nodes) > 30)
        size <- 0.5
    gr <- groupCalc1(topoFile)[[1]] %>% unlist
    sig <- nodes[!(nodes %in% gr)]
    nds <- c(gr, sig)
    freqData <- freqData %>% mutate(num = paste0(as.character(phenotype), 1:nrow(.)))
    freqGat <- freqData %>%
        gather(key = "Nodes", value = "Level", -phenotype, -num) %>%
        mutate(Nodes = factor(Nodes, levels = nds))
    breaks <- freqData %>% 
        mutate(n = 1:nrow(.)) %>%
        split(freqData$phenotype) %>% sapply(function(x){
            median(x$n) %>% round
        })
    breaks <- freqData$num[breaks]
    labels <- c("Epithelial", "Mesenchymal", "Hybrid")
    
    ggplot(freqGat, aes(x = Nodes, y = reorder(num, -as.numeric(phenotype)), fill = Level)) + 
        geom_tile() + 
        theme_Publication() + 
        theme(axis.title.x = element_blank(),
              # axis.text.y = element_text(angle = 90, hjust = 0.5),
              axis.text.x = element_text(angle = 90, 
                                         size = rel(size), vjust = 0.5, hjust = 1),
              legend.position = "right", 
              legend.direction = "vertical") + 
        scale_fill_manual(values = c("white", "black")) +
        scale_y_discrete(breaks= breaks,
                         labels = labels)+
        # scale_y_reverse()+
        labs(y = "")
    setwd(plotFolder)
    if (!dir.exists("StateMatrices"))
        dir.create("StateMatrices")
    setwd("StateMatrices")
    ggsave(filename = paste0(net, "_stateFreq.png"),
           width = 10, height = 6)
}

sapply(EMPNets, matrixPlot)

### Fig 4C,D - RACIPE frequency and state strength ----

distMakeRAC <- function(net, plot = T)
{
    net1 <- net
    if(str_detect(net1, "rand"))
        net1 <- str_remove(net1, "_rand.*")
    freqFile <- paste0(RACIPE_WT, "\\", net1, ".csv")
    freqDf <- read_csv(freqFile, col_types = cols())
    if (!any(colnames(freqDf) == "Strength"))
        return()
    df <- freqDf %>% select(Phenotype, Strength, Frequency, EMTScore, Frustration) %>% 
        set_names(c("phenotype", "Strength", "Frequency", "Score", "Frustration"))
    if (plot)
    {
        df1 <- df %>% mutate(Phenotype = ifelse(Score == 1, "Epithelial", "Hybrid")) %>%
            mutate(Phenotype = ifelse(Score == -1, "Mesenchymal", Phenotype))
        setwd(plotFolder)
        if (!dir.exists("RACIPE"))
            dir.create("RACIPE")
        setwd("RACIPE")
        ggplot(df1, aes(y = Strength, x = Phenotype)) + 
            geom_violin() +
            theme_Publication() + 
            # scale_fill_manual(values = c("Epithelial", "Hybrid", "Mesenchymal")) + 
            labs(y = "Stength (Sm)", x = "", color = "") +
            theme(legend.position = c(0.7, 0.9), legend.direction = "vertical")
        ggsave(paste0(net, "_scoreDistViolin.png"), width = 5.5, height = 5)
        ggplot(df1, aes(y = Frequency, x = Phenotype)) + geom_violin() + theme_Publication() + 
            labs(y = "SSF", x = "") + scale_y_log10()
        ggsave(paste0(net, "_FreqDistViolin.png"), width = 5.5, height = 5)
    }
    return(df %>% mutate(Net = net, Network = ifelse(str_detect(net, "rand"), "Rand", "WT")))
}
sapply(EMPNets, distMakeRAC)

### Fig 4E, S5E - State strength distribution ----
distMake <- function(net, plot = T)
{
    net1 <- net
    if(str_detect(net1, "rand"))
        net1 <- str_remove(net1, "_rand.*")
    freqFile <- paste0(randRaw, "\\", net1, "\\", net, "_finFlagFreq.csv")
    namez <- readLines(freqFile, n = 1) %>% str_split(",")%>% unlist
    if (!any(str_detect(namez,"Strength")))
        return()
    freqDf <- read_csv(freqFile, col_types = cols())
    
    df <- freqDf %>% select(phenotype, Strength)
    if (plot)
    {
        df1 <- df %>% mutate(Phenotype = ifelse(phenotype == "E", "Epithelial", "Hybrid")) %>%
            mutate(Phenotype = ifelse(phenotype == "M", "Mesenchymal", Phenotype))
        setwd(plotFolder)
        if (!dir.exists("StateStrength"))
            dir.create("StateStrength")
        setwd("StateStrength")
        ggplot(df1, aes(x = Strength, fill = Phenotype)) + 
            geom_histogram(aes(y = ..count../sum(..count..)), color = NA) +
            theme_Publication() + 
            # scale_fill_manual(values = c("Epithelial", "Hybrid", "Mesenchymal")) + 
            labs(x = "Stength (Sm)", y = "Frequency", color = "") +
            theme(legend.position = c(0.7, 0.9), legend.direction = "vertical")
        ggsave(paste0(net, "_strengthDist.png"), width = 5.5, height = 5)
    }
    return(df %>% mutate(Net = net, Network = ifelse(str_detect(net, "rand"), "Rand", "WT")))
}

sapply(EMPNets, distMake)

randStrengthDist <- function(net)
{
    setwd(randRaw)
    setwd(net)
    topoFiles <- list.files(".", ".topo")
    networks <- topoFiles %>% str_remove(".topo")
    networks <- c(sample(networks, 100), net)
    dfList <- lapply(networks, distMake, plot = F) %>% 
        reduce(rbind.data.frame)
    setwd(plotFolder)
    setwd("StateStrength")
    ggplot(dfList, aes(x = phenotype, y = Strength, fill  = Network)) +
        geom_violin() + 
        theme_Publication() + 
        labs(x = "Phenotype") + 
        theme(legend.position = c(0.5, 0.9))
    ggsave(paste0(net, "_randStrengthdist.png"), width = 6, height = 5)
    # write_csv(dfList, "randScoreDisct.csv", quote = "none")
}
sapply(EMPNets, randStrengthDist)

### Fig 4E - state strength vs stability metrics ----

sctterPlots <- function(net)
{
    freqDf <- read.csv(paste0(randRaw, "/", net, "/", net, "_finFlagFreq.csv")) %>%
        filter(flag == 1, !is.na(Avg0))
    cohDf <- read.csv(paste0(randcompiled, "/", net, "cohMnStd.csv"))
    colnames(cohDf)[1] <- "states"
    df <- merge(freqDf, cohDf, by = "states", all = T) %>% select(-contains("1")) %>%
        select(-contains("2")) %>% select(-flag)
    dfUseful <- df %>% select(Freq, Frust, Mean, Strength, Partial) %>% 
        set_names(c("Frequency", "Frustration", "Coherence", "Strength", "Partial"))
    setwd(plotFolder)
    if(!dir.exists("StateStrength"))
        dir.create("StateStrength")
    setwd("StateStrength")
    grob <- correlGrob(dfUseful, "Strength", "Frequency", xPos = 0.2)
    ggplot(dfUseful, aes(x = Strength, y = Frequency)) + geom_point() + 
        annotation_custom(grob) + labs(x = "Strength (Sm)", y = "SSF") +
        # geom_smooth(method = "lm") +
        theme_Publication()
    ggsave(paste0(net, "_frequency.png"), width = 5.5, height = 5)
    grobFreq <- grob$children[[1]]$label %>% str_replace(" :", "SSF :")
    grobFrust <- correlGrob(dfUseful, "Strength", "Frustration")
    grobFrust <- grobFrust$children[[1]]$label %>% str_replace(" :", "Frust :")
    grob <- grobTree(textGrob(paste0(grobFreq, "\n", grobFrust), 
                              x=0.55,  y=0.3, hjust=0,
                              gp=gpar(col="black", fontsize=18, fontface="bold")))
    
    ggplot(dfUseful, aes(x = Strength, y = Frequency, color = Frustration)) + 
        geom_point() + 
        annotation_custom(grob) + labs(x = "Strength (Sm)", y = "SSF") +
        # geom_smooth(method = "lm") +
        theme_Publication() +
        scale_color_viridis_c()+
        theme(legend.position = c(0.2,0.6), legend.direction = "vertical",
              legend.key.height = unit(0.6, "cm"))
    ggsave(paste0(net, "_frequencyFrustration.png"), width = 5.5, height = 5)
    
    grob <- correlGrob(dfUseful, "Strength", "Frustration")
    ggplot(dfUseful, aes(x = Strength, y = Frustration)) + geom_point() + 
        annotation_custom(grob) + labs(x = "Strength (Sm)") +
        # geom_smooth(method = "lm") +
        theme_Publication()
    ggsave(paste0(net, "_frustraion.png"), width = 5.5, height = 5)
    
    
    grob <- correlGrob(dfUseful, "Strength", "Coherence", xPos = 0.2)
    ggplot(dfUseful, aes(x = Strength, y = Coherence)) + geom_point() + 
        annotation_custom(grob) + labs(x = "Strength (Sm)") +
        # geom_smooth() +
        theme_Publication()
    ggsave(paste0(net, "_coherence.png"), width = 5.5, height = 5)
    
    # setwd(cwd)
    dfUseful$Network <- net
    print(net)
    return(dfUseful)
}
allNetsDf <- lapply(EMPNets, sctterPlots) %>% reduce(rbind.data.frame)
corDf <- sapply(EMPNets, function(x){
    df <- allNetsDf %>% filter(Network == x) %>% select(-Network) %>% cor %>%
        data.frame
    metrs <- df$Strength
    names(metrs) <- colnames(df)
    metrs <- metrs[-4]
    setwd(paste0(randRaw, "/", x))
    topoFile <- paste0(x, ".topo")
    g <- groupCalc(topoFile)[-5] %>% as.numeric %>% abs %>% mean
    setwd(cwd)
    return(c(metrs, GroupStrength = g))
}) %>% t %>% data.frame %>% mutate(Network = netNameKey[EMPNets])

cor2 <- corDf %>% gather(key = "Metric", value = "Correlation", -Network) %>% 
    mutate(Network = factor(Network, levels = c("57N 113E", "20N 40E", "18N 33E", 
                                                "26N 100E", "22N 82E"))) %>%
    mutate(Metric = ifelse(Metric == "Frequency", "SSF", Metric))
ggplot(cor2 %>% filter(!(Metric %in% c("Partial", "GroupStrength"))), 
       aes(y = Metric, x = Network, fill = Correlation)) + geom_tile() + 
    theme_Publication() + scale_fill_gradient2(low = "red", high = "blue", limits = c(-1,1)) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1),
          axis.text = element_text(face = "bold", size = rel(1.2)),
          legend.position = "right", legend.direction = "vertical",
          legend.key.height = unit(1, "cm"))
setwd(plotFolder)
setwd("StateStrength")
ggsave("StrengthCorrel.png", width = 6, height = 3.5)

ggplot(cor2 %>% filter(Metric == "GroupStrength") %>% mutate(Metric = "Mean\nGroup Strength"), 
       aes(y = Metric, x = Network, fill = Correlation)) + geom_tile() + 
    theme_Publication() + scale_fill_viridis_c() + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1),
          axis.text = element_text(face = "bold", size = rel(1.3)),
          axis.text.y = element_text(hjust = 0.5),
          legend.position = "right", legend.direction = "vertical",
          legend.key.height = unit(1, "cm")) +
    labs(x = "", y = "", fill = "Mean\nGroup\nStrength")
ggsave("Gs.png", height = 2.3, width = 6)


### Fig 4E - Strength correlations random networks ----

correlations <- function(net)
{
    setwd(randRaw)
    setwd(net)
    topoFiles <- list.files(".", ".topo")
    freqFiles <- str_replace(topoFiles, ".topo", "_finFlagFreq.csv")
    dat <- lapply(freqFiles, function(x){
        print(x %>% str_remove("_finFlagFreq.csv"))
        d <- read_csv(x, show_col_types = F, col_types = "c") %>% 
            filter(!is.na(Avg0), flag ==1)
        if (nrow(d) < 3 || (!any(colnames(d) == "Strength")))
            return()
        df <- d %>% select(Avg0, frust0, Strength)
        if (any(colnames(d) == "coherence0"))
            df <- d %>% select(Avg0, frust0, Strength, coherence0)
        if (any(colnames(d) == "coherence0.x"))
            df <- d %>% select(Avg0, frust0, Strength, contains("coherence")) %>%
            set_names(c("Avg0", "frust0", "Strength", "coherence0", "dummy")) %>%
            select(-dummy)
        if(all(is.na(df$coherence0)))
            df <- df %>% select(-coherence0)
        df <- df %>% gather(key = "Metric", value = "Value", -Strength) %>%
            mutate(Value = as.numeric(Value))
        metricKey <- c("SSF", "Frustration", "Strength", "Coherence")
        names(metricKey) <- c("Avg0", "frust0", "Strength", "coherence0")
        df$Metric <- metricKey[df$Metric]
        df$Value[any(is.logical(df$Value))] <- NA
        multiFactorCorrelation(df, "Metric", "Strength", "Value", label = F, method = "spearman") %>%
            mutate(Net = str_remove(x, "_finFlagFreq.csv"), 
                   Significance = ifelse(pValue < 0.05, "", "*"))
    }) %>% reduce(rbind.data.frame)
    groupDat <- read.csv(paste0(net, "_rand_groups.csv"))
    dat <- merge(dat, groupDat %>% select(Net, Gs), by = "Net", all.x = T)
    write.csv(dat, paste0(randRaw, "\\", net, "_stateStrengthCorrel.csv"), row.names = F)
}
sapply(EMPNets[5], correlations)

correlationPlots <- function(net)
{
    dat <- read_csv(paste0(randRaw, "/", net, "_stateStrengthCorrel.csv"), show_col_types = F)
    Metrics <- unique(dat$Factors)
    setwd(plotFolder)
    setwd("StateStrength")
    sapply(Metrics, function(x){
        print(x)
        d <- dat %>% filter(Factors == x) %>% select(Gs, Correlation) %>% drop_na
        cVal <- cor(d$Gs, d$Correlation)
        ypos <- NULL
        if (cVal>0)
            ypos <- 0.1
        grob <- correlGrob(d, "Gs", "Correlation", yPos = ypos, method = "spearman")
        ggplot(d, aes(x = Gs, y = Correlation)) +
            geom_point() + annotation_custom(grob) +
            theme_Publication() + labs(x = "Mean Group Strength", y = x)
        ggsave(paste0(net, "_rand_", x, ".png"), width = 5.5, height = 5)
    })
    multiFactorCorrelation(dat, "Factors", "Gs", "Correlation", method = "spearman", label = F) %>% 
        mutate(Network = netNameKey[net], Sig = ifelse(pValue < 0.05, "", "X"))
}
df <- lapply(EMPNets, correlationPlots) %>% reduce(rbind.data.frame)
ggplot(df, aes(x = Network, y = Factors, fill = Correlation)) + 
    geom_tile() + geom_text(aes(label = Sig)) +
    theme_Publication() + scale_fill_gradient2(low = "red", high = "blue", limits = c(-1,1)) +
    labs(x = "", y= "") + 
    theme(legend.position = "right", legend.direction = "vertical", legend.key.height = unit(0.8,"cm"), 
          axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1),
          axis.text = element_text(face = "bold", size = rel(1.4)))
ggsave("StateStrengthRand.png", width = 6.5, height = 5)
