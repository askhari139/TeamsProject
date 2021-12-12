library(tidyverse)
source("D:/Teams/TeamsPaperFigures/plot_theme.R")
options(stringsAsFactors = F)
library(compiler)

topo_to_int_mat <- function(topo_file) {
    # print(topo_file)
    df <- read.delim(topo_file, sep = " ", stringsAsFactors = F)
    if (ncol(df) != 3) {
        df <- read.delim(topo_file, stringsAsFactors = F)
    }
    # browser()
    colnames(df) <- c("Source", "Target", "Type")
    df <- df %>% 
        mutate(Source = str_remove_all(Source, "\\s")) %>%
        mutate(Target = str_remove_all(Target, "\\s")) %>%
        mutate(Type = ifelse(Type == 2, -1, 1))
    
    nodes <- unique(c(df$Source, df$Target)) %>% 
        sort(decreasing = T)
    n_nodes <- length(nodes)
    intmat <- rep(0, n_nodes * n_nodes) %>% 
        matrix(ncol = n_nodes)
    df1 <- df %>% 
        mutate(Source = sapply(Source, function(x) {which(nodes == x)})) %>% 
        mutate(Target = sapply(Target, function(x) {which(nodes == x)}))
    # browser()
    dummy <- apply(df1, 1, function(x) {
        # browser()
        i <- x[1]
        j <- x[2]
        k <- x[3]
        intmat[i,j] <<- k
    })
    return(list(intmat, nodes))
}

compute_power_matrix <- function(mat, power) {
    res <- mat
    if (power == 1)
    {
        return(res)
    }
    for (i in 2:power) {
        res <- res %*% mat
    }
    return(res)
}

influence_matrix <- function(topoFile, lmax = 10) {
    ls <- topo_to_int_mat(topoFile)
    intmat <- ls[[1]]
    intmax <- intmat
    intmax[which(intmax == -1)] <- 1
    res <- 0
    for (l in 1:lmax) {
        intM <- compute_power_matrix(intmat, l)
        maxM <- compute_power_matrix(intmax, l)
        r1 <- intM / maxM
        r1[is.nan(r1)] <- intM[is.nan(r1)]
        res <- res + r1
    }
    res <- res / lmax
    
    nodes <- ls[[2]]
    
    influence_mat <- res
    colnames(influence_mat) <- rownames(influence_mat) <- nodes
    signal <- which(apply(intmat, 2, function(x){all(x==0)}))
    output <- which(apply(intmat, 1, function(x){all(x==0)}))
    secondary_signal <- which(apply(intmat, 2, function(x){
        if (length(signal) !=0)
            all(x[-signal] == 0)
        else
            F
    }))
    secondary_output <- which(apply(intmat, 1, function(x){
        if (length(output) != 0)
            all(x[-output] == 0)
        else
            F
    }))
    nonEssentials <- c(signal, output, secondary_signal, secondary_output)
    if(length(nonEssentials))
    {
        influence_reduced <- influence_mat[-nonEssentials, 
                                           -nonEssentials]
        nodes_reduced <- nodes[-nonEssentials] %>% str_replace_all(regex("\\W+"), "")
    }
    else
    {
        influence_reduced <- influence_mat
        nodes_reduced <- nodes %>% str_replace_all(regex("\\W+"), "")
    }
    if (length(nodes_reduced) < 2)
        return()
    rownames(influence_reduced) <- colnames(influence_reduced) <- nodes_reduced
    
    
    influence_reduced <- influence_reduced[order(influence_reduced[,1]), order(influence_reduced[1,])]
    net <- str_remove(topoFile, ".topo") %>% paste0("Influence/",.)
    if (!dir.exists("Influence"))
        dir.create("Influence")
    write.csv(influence_reduced, paste0(net, ".csv"))
    write.csv(influence_mat, paste0(net, "_fullMat.csv"))
}
influence_matrix <- cmpfun(influence_matrix)

groupPlotter <- function(topoFile, remNodes = NULL)
{
    # print(topoFile)
    inflFile <- paste0("Influence/", str_replace(topoFile, ".topo", "_fullMat.csv"))
    if (!file.exists(inflFile))
        return()
    inflMat <- read.csv(inflFile)
    nodes <- inflMat[[1]]
    rownames(inflMat) <- nodes
    df <- as.matrix(inflMat[,-1])
    if (!is.null(remNodes))
    {
        ids <- which(nodes %in% remNodes)
        nodes <- nodes[-ids]
        df <- df[-ids, -ids]
    }
    df1 <- apply(df, 2, function(x){
        ifelse(x > 0, 1, -1)
    })
    df1 <- cbind(df1, t(df1))
    hc <- hclust(dist(df1))
    clust <- cutree(hc, 2)
    g1 <- nodes[clust == 1]
    g2 <- nodes[clust == 2]
    nOrder <- c(g1,g2)
    df2 <- data.frame(df) %>% mutate(nodes1 = nodes) %>%
        gather(key = "Nodes", value = "Influence", -nodes1) %>%
        mutate(nodes1 = factor(nodes1, levels = nOrder), Nodes = factor(Nodes, levels = nOrder))
    ggplot(df2, aes(x = Nodes, y = nodes1, fill = Influence)) + geom_tile() +
        theme_Publication() + scale_fill_gradient2(low = "blue", high = "red", limits = c(-1,1)) + 
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 90), legend.position = "right", 
              legend.direction = "vertical", legend.key.height = unit(0.5, "in"))
    ggsave(paste0("Influence/", str_replace(topoFile, ".topo", "_group.png")), width = 7.5, height = 6)
}
groupPlotter <- cmpfun(groupPlotter)

groupCalc <- function(topoFile, remNodes = NULL)
{#browser()
    inflFile <- paste0("Influence/", str_replace(topoFile, ".topo", ".csv"))
    if (!file.exists(inflFile))
        return(c(NA, NA, NA, NA, str_remove(topoFile , ".topo")))
    inflMat <- read.csv(inflFile)
    
    nodes <- inflMat[[1]]
    rownames(inflMat) <- nodes
    df <- as.matrix(inflMat[,-1])
    if (!is.null(remNodes))
    {
        ids <- which(nodes %in% remNodes)
        nodes <- nodes[-ids]
        df <- df[-ids, -ids]
    }
    df1 <- apply(df, 2, function(x){
        ifelse(x > 0, 1, -1)
    })
    df1 <- cbind(df1, t(df1))
    hc <- hclust(dist(df1))
    clust <- cutree(hc, 2)
    g1 <- nodes[clust == 1] %>% sort
    g2 <- nodes[clust == 2] %>% sort
    if(g1[length(g1)]> g2[1])
    {
        g0 <- g1
        g1 <- g2
        g2 <- g0
    }
    nOrder <- c(g1,g2)
    df2 <- data.frame(df) %>% mutate(nodes1 = nodes) %>%
        gather(key = "Nodes", value = "Influence", -nodes1) %>%
        mutate(nodes1 = factor(nodes1, levels = nOrder), Nodes = factor(Nodes, levels = nOrder))
    
    g11 <- df[g1,g1] %>% as.vector %>% mean
    g22 <- df[g2,g2] %>% as.vector %>% mean
    g12 <- df[g1,g2] %>% as.vector %>% mean
    g21 <- df[g2,g1] %>% as.vector %>% mean
    return(c(g11, g22, g12, g21, str_remove(topoFile , ".topo")))
}
groupCalc <- cmpfun(groupCalc)

# topoFiles <- list.files(".", ".topo$")
# dummy <- sapply(topoFiles, influence_matrix)
# # sapply(topoFiles, groupPlotter)
# df <- lapply(topoFiles, groupCalc) %>% reduce(rbind.data.frame) %>% set_names(c("G11", "G22", "G12", "G21", "Net"))
# write.csv(df, "EMT_RACIPE_edel_groups.csv", row.names = F)
