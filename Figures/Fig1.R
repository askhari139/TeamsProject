mainFolder <- "D:/Teams"
plotFolder <- paste0(mainFolder, "/Figures/Fig1")
source(paste0(mainFolder, "/Figures/figCore.R"))
source(paste0(mainFolder, "/Final_Results/codes/stateLabeller.R"))

library(igraph)

### 1A, S1 - preparing networks for cytoscape ----

GMLcreate <- function(net)
{
    topoFile <- paste0(mainFolder, "/Final_Results/WT_topoFiles/", net, ".topo")
    topoDf <- read.delim(topoFile, sep = "") %>%
        mutate(Source = Source %>% str_replace_all(regex("\\W+"), ""),
               Target  = Target %>% str_replace_all(regex("\\W+"), ""))
    gr <- graph_from_data_frame(topoDf)
    write_graph(gr, paste0("D:\\Teams\\Final_Results\\WT_topoFiles\\", net, ".GML"),
                format ="gml")
    write.csv(topoDf, paste0("D:\\Teams\\Final_Results\\WT_topoFiles\\", net, "_net.csv"), 
              row.names = F, quote = F)
    nodeChar <- getEMSONodes(net)
    df <- data.frame(Node = nodeChar, Nature = names(nodeChar) %>% str_remove_all("\\d"))
    write.csv(df, paste0("D:\\Teams\\Final_Results\\WT_topoFiles\\", net, "_nodes.csv"), 
              row.names = F, quote = F)
}
sapply(EMPNets, GMLcreate)

### 1B, S2A: Interaction and influence matrix plot ----

matPlot <- function(net)
{
    net1 <- str_remove(net, "_rand.*")
    setwd(paste0(randRaw, "/",net1))
    topoFile <- paste0(net,".topo")
    influenceMat <- read.csv(paste0("Influence/", net, "_fullMat.csv"))
    nodes <- influenceMat[[1]] %>% str_replace_all(regex("\\W+"), "")
    influenceMat[[1]] <- nodes
    colnames(influenceMat) <- c("Nodes", nodes)
    ls <- topo_to_int_mat(topoFile)
    intMat <- ls[[1]]
    nodes <- ls[[2]] %>% str_replace_all(regex("\\W+"), "")
    colnames(intMat) <- rownames(intMat) <- nodes
    intMat <- data.frame(intMat) %>% mutate(Nodes = nodes)
    nodes <- getEMSONodes(net)
    inflDf <- influenceMat %>% gather(key = "Nodes1", value = "Influence", -Nodes) %>%
        mutate(Nodes1 = factor(Nodes1, levels = nodes), Nodes = factor(Nodes, levels = nodes))
    intDf <- intMat %>% gather(key = "Nodes1", value = "Influence", -Nodes) %>%
        mutate(Nodes1 = factor(Nodes1, levels = nodes), Nodes = factor(Nodes, levels = nodes))
    if (length(nodes) > 30)
    {
        inflDf <- inflDf %>% mutate(Nodes = as.integer(Nodes), Nodes1 = as.integer(Nodes1))
        intDf <- intDf %>% mutate(Nodes = as.integer(Nodes), Nodes1 = as.integer(Nodes1))
    }
    setwd(plotFolder)
    if (!dir.exists("Matrices"))
        dir.create("Matrices")
    setwd("Matrices")
    ggplot(inflDf, aes(x = Nodes1, y = Nodes, fill = Influence)) + geom_tile() +
        theme_Publication() + scale_fill_gradient2(low = "blue", high = "red", limits = c(-1,1)) + 
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 90), legend.position = "right", 
              legend.direction = "vertical", legend.key.height = unit(0.5, "in"))
    ggsave(paste0(net, "_influence.png"), width = 7.5, height = 6)
    ggplot(intDf, aes(x = Nodes1, y = Nodes, fill = Influence)) + geom_tile() +
        theme_Publication() + scale_fill_gradient2(low = "blue", high = "red", limits = c(-1,1)) + 
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 90), legend.position = "right", 
              legend.direction = "vertical", legend.key.height = unit(0.5, "in"))
    ggsave(paste0(net, "_interaction.png"), width = 7.5, height = 6)
}
sapply(EMPNets, matPlot)

### 1C, S2B - Gs distributions ----
densPlot <- function(net)
{
    setwd(randcompiled)
    df <- read.csv(paste0(net, "_rand_groups.csv"))
    df <- df %>% mutate(GrpStr = df %>% select(-Net) %>% 
                            apply(1, function(x){x %>% abs %>% mean}))
    wtNet <- which(!str_detect(df$Net, "rand"))
    wtStr <- df$GrpStr[wtNet]
    ggplot(df, aes(x = GrpStr)) + geom_histogram(aes(y = ..count../sum(..count..))) + 
        geom_vline(xintercept = wtStr, color = "red", size = 1.5) + 
        theme_Publication() + 
        labs(x = "Mean Group Strength", y = "Frequency")
    setwd(plotFolder)
    if (!dir.exists("GsDist"))
        dir.create("GsDist")
    setwd("GsDist")
    ggsave(paste0(net, "_rand_Gs.png"), width = 5.5, height = 5)
}

sapply(EMPNets, densPlot)

### S2C,D ----

### 1D, S2E,F - correlation matrices ----

correlationMatrixRAC <- function(net)
{
    setwd(paste0(RACIPE_WT, "/", net))
    topoFile <- paste0(net, ".topo")
    
    
    nodes <- read.delim(paste0(net, ".prs")) %>% 
        filter(str_detect(Parameter, "Prod")) %>% 
        select(Parameter) %>% unlist %>% str_remove("Prod_of_")
    df <- read_delim(paste0(net, "_solution.dat"), delim = "\t", col_names = F) %>%
        set_names(c("ParIndex", "nStates", "Count", nodes)) %>%
        select(-ParIndex, -nStates, -Count)
    setwd(cwd)
    corDf <- cor(df) 
    nodes <- getEMSONodes(net)
    corDf <- corDf %>% data.frame %>% mutate(Nodes = rownames(.) %>% 
                                                 str_replace_all(regex("\\W+"), "")) %>%
        gather(key = "Nodes1", value = "Correlation", -Nodes) %>% 
        mutate(Nodes = factor(Nodes, levels = nodes), 
               Nodes1 = factor(Nodes1 %>% str_replace_all(regex("\\W+"), ""), 
                               levels = nodes))
    
    ggplot(corDf, aes(x = Nodes, y = Nodes1, fill = Correlation)) + geom_tile() +
        theme_Publication() + scale_fill_gradient2(low = "blue", high = "red", 
                                                   limits = c(-1,1)) + 
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "right", 
              legend.direction = "vertical", legend.key.height = unit(0.5, "in"))
    setwd(plotFolder)
    if (!dir.exists("CorMats"))
        dir.create("CorMats")
    ggsave(paste0("CorMats/", net, "_RAC_correlation.png"), 
           width = 7.5, height = 6)
}
sapply(EMPNets, correlationMatrixRAC)

correlationMatBool <- function(topoFile)
{#browser()
    print(topoFile)
    nodes <- readLines(topoFile %>% str_replace(".topo", "_nodes.txt"))
    corMat <- read_csv(topoFile %>% str_replace(".topo", "_finFlagFreq.csv"), 
                       show_col_types = F) %>% 
        filter(flag == 1) %>% select(states, Avg0) %>% drop_na %>%
        mutate(Avg0 = Avg0*10000 %>% round) %>% apply(1, function(x){
            s <- x[1] %>% str_remove_all("'")
            n <- x[2] %>% as.integer
            rep(s, n)
        }) %>% unlist 
    if (length(corMat) < 3)
        return(NA)
    corMat <- corMat %>% lapply(function(x){
        str_split(x, "") %>% unlist %>% as.integer
    }) %>% reduce(rbind.data.frame) %>% set_names(nodes) %>% cor
    colnames(corMat) <- colnames(corMat) %>% str_replace_all(regex("\\W+"), "")
    rownames(corMat) <- rownames(corMat) %>% str_replace_all(regex("\\W+"), "")
    return(corMat)
}
correlationMatBool <- cmpfun(correlationMatBool)
corMatPlot <- function(net)
{
    netDir <- paste0(randRaw, "/", net)
    setwd(netDir)
    topoFile <- paste0(net, ".topo")
    corMat <- correlationMatBool(topoFile)
    nodes <- getEMSONodes(net)
    setwd(cwd)
    corDf <- corMat %>% data.frame %>% mutate(Nodes = rownames(.) %>% 
                                                  str_replace_all(regex("\\W+"), "")) %>%
        gather(key = "Nodes1", value = "Correlation", -Nodes) %>% 
        mutate(Nodes = factor(Nodes, levels = nodes), 
               Nodes1 = factor(Nodes1 %>% str_replace_all(regex("\\W+"), ""), 
                               levels = nodes))
    
    ggplot(corDf, aes(x = Nodes, y = Nodes1, fill = Correlation)) + geom_tile() +
        theme_Publication() + scale_fill_gradient2(low = "blue", high = "red", 
                                                   limits = c(-1,1)) + 
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5), 
              legend.position = "right", 
              legend.direction = "vertical", legend.key.height = unit(0.5, "in"))
    setwd(plotFolder)
    if (!dir.exists("CorMats"))
        dir.create("CorMats")
    ggsave(paste0("CorMats/", net, "_bool_correlation.png"), 
           width = 7.5, height = 6)
}
sapply(EMPNets, corMatPlot)

### 1F, S2G  - Distance between correlation and influence ----

randCorrelPlot <- function(net)
{
    setwd(paste0(randRaw, "/", net))
    distDat <- read.csv(paste0(net, "_correlInfl.csv"))
    groupDat <- read.csv(paste0("../", net, "_rand_groups.csv")) %>%
        mutate(Gs = (abs(G11) + abs(G22) + abs(G12) + abs(G21))/4) %>%
        select(Net, Gs) %>% set_names(c("Network", "Gs"))
    df <- merge(distDat, groupDat, by = "Network", all = T) %>%
        mutate(Net = netNameKey[net])
    grob <- correlGrob(df, "Dist", "Gs")
    WT <- df %>% filter(Network == net)
    ggplot(df, aes(x = Gs, y = Dist*100)) + geom_point(size = 2) + 
        # geom_smooth(method = "lm") + 
        geom_label_repel(data = WT, color = "red", label = "WT") + 
        facet_wrap(~Net, ncol = 1) + 
        annotation_custom(grob) + theme_Publication() +
        theme(legend.position = "none") +
        labs(x = "Mean Group Strength", y = "Percent difference between\nCorrelation and Influence")
    setwd(plotFolder)
    if (!dir.exists("CorrelVInfl"))
        dir.create("CorrelVInfl")
    setwd("CorrelVInfl")
    ggsave(paste0(net, "_correlVInfl.png"), width = 5.5, height = 5.1)
    
    ggplot(df, aes(x = Dist)) + 
        geom_histogram(aes(y = ..count../sum(..count..))) +
        geom_vline(xintercept = WT$Dist, size = 1.5, color = "red") + 
        theme_Publication() + labs(x = "Distance between\nCorrelation and Influence", 
                                   y = "Frequency")
    ggsave(paste0(net, "_correlVInfldist.png"))
}
sapply(EMPNets, randCorrelPlot)
