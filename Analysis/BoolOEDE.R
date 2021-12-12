library(tidyverse)
library(stringr)
library(future)
library(future.apply)
library(data.table)
library(compiler)
source("D:\\Teams\\TeamsPaperFigures\\figCore.R")
source("D:\\Teams\\Final_Results\\codes\\inflMat.R")
source("D:\\Teams\\Final_Results\\codes\\stateLabeller.R")


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
hamming <- function(x)
{
    x <- str_split(x, "")
    sum(x[[1]]!=x[[2]])/(length(x[[1]])-2)
} 
hamming <- cmpfun(hamming)

dec2bin <- function(x, l=32) 
{
    paste(rev(rev(as.integer(rev(intToBits(x))))[1:l]), collapse = "")
} %>% cmpfun


inwards <- function(x,df){
    colnames(df) <- c("State", "Step")
    df %>% filter(Step == x) %>% select(State) %>% unlist
} %>% cmpfun

asyncInit_OE <- function(sInit, update_matrix, maxT = 1000, init = NULL,
                         OE = T, ID = NULL)
{#browser()
    f <- 0 #flag
    N_nodes <- length(sInit)
    # update_matrix <- 2*update_matrix + diag(N_nodes)
    indices <- 1:N_nodes
    up_ind_list <- sample(indices[-ID], maxT, replace = T) # generate the update indices 
    s <- sInit
    for (k in 1:maxT){
        s_dummy <- s%*%update_matrix %>% sign # update the state
        if (OE)
            s_dummy[ID] <- 1
        else
            s_dummy[ID] <- -1
        if (all(s_dummy == s)) 
            f <- 1 # flag 1 means it is a proper steady state
        
        up_ind <- up_ind_list[k]
        s[up_ind] <- s_dummy[up_ind]
        
        if (f == 1) 
            break
    }
    fin <- paste0("'", 
                  paste0(
                      replace(s, which(s==-1),0), 
                      collapse = ""),
                  "'")
    return(c(fin, f))
}

asyncInit_OE <- cmpfun(asyncInit_OE)

asyncSim <- function(net, nInit = 10000, maxT = 1000, OE = T)
{
    # browser()
    setwd(paste0(randRaw, "/", net))
    if (!dir.exists("OEDE"))
        dir.create("OEDE")
    nodeOrder <- readLines(paste0(net, "_nodes.txt"))
    nam <- paste0("OEDE/", net , "_", ifelse(OE, "OE", "DE"))
    ls <- topo_to_int_mat(paste0(net, ".topo"))
    intMat <- ls[[1]]
    nodes <- ls[[2]]
    colnames(intMat) <- rownames(intMat) <- nodes
    intMat <- intMat[nodeOrder, nodeOrder]
    nNodes <- length(nodes)
    plan(multisession, workers = 8)
    dummyOE <- future_sapply(nodeOrder, function(x){
        upm <- 2*intMat + diag(nNodes)
        id <- which(nodeOrder == x)
        dfList <- lapply(1:3, function(i){
            sInitList <- sample(c(-1, 1), nInit*nNodes, replace = T) %>%
                matrix(ncol = nNodes)
            sInitList[, id] <- ifelse(OE, 1, -1)
            apply(sInitList, 1, function(y){
                asyncInit_OE(y, upm,maxT, OE = OE, ID = id)
            }) %>% t %>% data.frame %>% set_names(c("states", "flag")) %>% 
                group_by(states, flag) %>% summarise(Freq = n(), .groups = "drop") %>% 
                mutate(Freq = Freq/sum(Freq))
        }) %>% reduce(merge, by = c("states", "flag"), all = T)
        dfList[is.na(dfList)] <- 0
        dfList <- dfList %>% mutate(Avg0 = dfList %>% select(contains("Freq")) %>% apply(1,mean),
                                    SD0 = dfList %>% select(contains("Freq")) %>% apply(1,sd)) %>%
            select(states, flag, Avg0, SD0)
        write_csv(dfList, paste0(nam, "_", x, ".csv"), quote = "none")
        
    })
    future:::ClusterRegistry("stop")
    
}
asyncSim <- cmpfun(asyncSim)
