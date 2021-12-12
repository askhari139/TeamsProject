source("D:\\Teams\\Final_Results\\codes\\inflMat.R")

groupCalc1 <- function(topoFile, remNodes = NULL)
{#browser()
    # print(topoFile)
    inflFile <- paste0("Influence/", str_replace(topoFile, ".topo", ".csv"))
    if (!file.exists(inflFile))
        influence_matrix(topoFile)
    if (!file.exists(inflFile))
        return(list(NA, NA))
    inflMat <- read.csv(inflFile, stringsAsFactors = F)
    
    # nodes <- readLines(str_replace(topoFile, ".topo", "_nodes.txt"))
    nodes <- rownames(inflMat) <- colnames(inflMat)[-1]
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
    l <- list(g1, g2)
    mirdetect <- c(sum(str_detect(g1 %>% str_to_upper, "MIR")), 
                   sum(str_detect(g2 %>% str_to_upper, "MIR")))
    egroup <- which.max(mirdetect)
    mgroup <- ifelse(egroup == 1, 2, 1)
    names(l)[c(egroup, mgroup)] <- c("E", "M")
    return(list(l, list(nodes)))
}
groupCalc1 <- cmpfun(groupCalc1)

scoreCalc <- function(state, nodes, group, inflMat)
{
    # browser()
    state <- state %>% str_split("") %>% unlist
    if (state[1] == "'")
      state <- state[-c(1, length(state))]
    names(state) <- nodes %>% unlist
    stategroups <- state[group %>% unlist]
    sOn <- nodes[state == "1"]
    sgOn <- unlist(group)[stategroups == "1"]
    EOn <- nodes[which(nodes %in% group$E[state[group$E] == "1"])]
    MOn <- nodes[which(nodes %in% group$M[state[group$M] == "1"])]
    rownames(inflMat) <- inflMat[,1] %>% str_replace_all(regex("\\W+"), "")
    inflMat <- as.matrix(inflMat[,-1])
    colnames(inflMat) <- colnames(inflMat) %>% str_replace_all(regex("\\W+"), "")
    c(Strength = sum(inflMat[sOn, sOn]), Partial = sum(inflMat[sgOn, sgOn]),
                  Epithelial = sum(inflMat[EOn, EOn]), 
      Mesenchymal = sum(inflMat[MOn, MOn]),
      EMTScore = (length(EOn)/ length(group$E)) - length(MOn)/length(group$M))
}
scoreCalc <- cmpfun(scoreCalc)

stateLabel <- function(state, nodes, stateNodes, groups)
{#browser()
    state <- str_split(state, "") %>% unlist
    if (state[1] == "'")
      state <- state[-c(1, length(state))]
    s <- stateNodes[which(state == "1")]
    nodeGroups <- s
    if (all(groups[[1]] %in% nodeGroups))
        return(names(groups)[1])
    else if (all(groups[[2]] %in% nodeGroups))
        return(names(groups)[2])
    else
        return("H")
}
stateLabel <- cmpfun(stateLabel)
stateLabel <- Vectorize(stateLabel, vectorize.args = "state")

labeller <- function(topoFile, rand = T, write = F){#browser()
    print(topoFile)
    net <- topoFile %>% str_remove(".topo")
    
    freqFile <- paste0(net, "_finFlagFreq.csv")
    #nodes <- readLines(paste0(net, "_nodes.txt")) %>% str_replace_all(regex("\\W+"), "")
    x  <- groupCalc1(topoFile)
    groups <- x[[1]]
    nodes <- x[[2]]
    stateNodes <- readLines(paste0(net, "_nodes.txt")) %>% str_replace_all(regex("\\W+"), "")
    freqDf <- read_csv(freqFile, show_col_types = F, lazy = F)
    if (ncol(freqDf) > 11)
        freqDf <- freqDf[, 1:11]
    freqDf$phenotype <- stateLabel(freqDf$states, nodes,stateNodes, groups)
    inflfile <- paste0("Influence/", net, "_fullMat.csv")
    # browser()
    cond <- file.exists(inflfile) && !str_detect(topoFile, "rand")
    if (rand)
        cond <- file.exists(inflfile)
    if (cond)
    {
        inflMat <- read.csv(inflfile)
        freqDf <- cbind.data.frame(freqDf, 
                                   sapply(freqDf$states, function(x){scoreCalc(x, 
                                       nodes = stateNodes, group = groups, 
                                       inflMat = inflMat)}) %>% t)
    }
    if (write)
      write.csv(freqDf, freqFile, row.names = F)
    else
      return(freqDf)
    df <- freqDf %>% select(states, Avg0, phenotype)
    df <- df[complete.cases(df),]
    hybridNess <- sum(df$phenotype == "H")/nrow(df)
    hybridFreq <- sum(df$Avg0[df$phenotype == "H"])/sum(df$Avg0)
    print(topoFile)
    return(c(net, hybridNess, hybridFreq))
}
labeller <- cmpfun(labeller)

labeller_states <- function(states, topoFile,rand = F){#browser()
    print(topoFile)
    net <- topoFile %>% str_remove(".topo")
    if(!file.exists(topoFile))
      setwd(paste0(randRaw, "/",net))
    
    # freqFile <- paste0(net, "_finFlagFreq.csv")
    #nodes <- readLines(paste0(net, "_nodes.txt")) %>% str_replace_all(regex("\\W+"), "")
    x  <- groupCalc1(topoFile)
    groups <- x[[1]]
    nodes <- x[[2]]
    stateNodes <- readLines(paste0(net, "_nodes.txt")) %>% str_replace_all(regex("\\W+"), "")
    
    phenotype <- stateLabel(states, nodes,stateNodes, groups)
    return(phenotype)
}

labeller_states_RACIPE <- function(states, topoFile,rand = F){#browser()
  print(topoFile)
  net <- topoFile %>% str_remove(".topo")
  if (!file.exists(topoFile))
    file.copy(paste0(RACIPE_WT, "/", net, "/", net, ".topo"), topoFile)
  # setwd(paste0(randRaw, "/",net))
  
  # freqFile <- paste0(net, "_finFlagFreq.csv")
  #nodes <- readLines(paste0(net, "_nodes.txt")) %>% str_replace_all(regex("\\W+"), "")
  x  <- groupCalc1(topoFile)
  groups <- x[[1]]
  nodes <- x[[2]]
  prsFile <- paste0(net, ".prs")
  if (!file.exists(prsFile))
    prsFile <- paste0(RACIPE_WT, "/", net, "/", net, ".prs")
  stateNodes <- read.delim(prsFile) %>% 
    filter(str_detect(Parameter, "Prod_")) %>% select(Parameter) %>% 
    unlist %>% str_remove("Prod_of_") %>%
    str_replace_all(regex("\\W+"), "")
  
  inflFile <- paste0("Influence/",net, "_fullMat.csv")
  if (!file.exists(inflFile))
    inflFile <- paste0(randRaw, "/", net, "/Influence/", net, "_fullMat.csv")
  inflMat <- read.csv(inflFile)
  
  
  score <- sapply(states, scoreCalc, nodes = stateNodes, group = groups, 
                  inflMat = inflMat) %>% t %>% data.frame %>%
    mutate(Phenotype = "H", State = states) %>%
    mutate(Phenotype = ifelse(EMTScore > 0.5, "E", Phenotype)) %>% 
    mutate(Phenotype = ifelse(EMTScore < -0.5, "M", Phenotype)) 
  
  return(score)
}
labeller_states <- cmpfun(labeller_states)
labeller_states_RACIPE <- cmpfun(labeller_states_RACIPE)


score_states <- function(states, topoFile)
{
  print(topoFile)
  net <- topoFile %>% str_remove(".topo")
  pwd <- getwd()
  setwd(paste0(randRaw, "/", net))
  freqFile <- paste0(net, "_finFlagFreq.csv")
  #nodes <- readLines(paste0(net, "_nodes.txt")) %>% str_replace_all(regex("\\W+"), "")
  x  <- groupCalc1(topoFile)
  group <- x[[1]]
  nodes <- x[[2]]
  stateNodes <- readLines(paste0(net, "_nodes.txt")) %>% str_replace_all(regex("\\W+"), "")
  s <- sapply(states, function(state){
    state <- state %>% str_split("") %>% unlist
    if (state[1] == "'")
      state <- state[-c(1, length(state))]
    names(state) <- stateNodes %>% unlist
    EOn <- stateNodes[which(stateNodes %in% group$E[state[group$E] == "1"])]
    MOn <- stateNodes[which(stateNodes %in% group$M[state[group$M] == "1"])]
    (length(EOn)/ length(group$E)) - (length(MOn)/length(group$M))
  }) 
  setwd(pwd)
  return(s)
}
score_states <- cmpfun(score_states)

hybridMetric <- function(dir, nam = NULL, rand = F){
    cwd <- getwd()
    setwd(dir)
    topoFiles <- list.files(".", ".topo$")
    df <- sapply(topoFiles, labeller, rand = rand) %>% t%>%
        data.frame %>% set_names(c("Net", "hybridness", "hybridFreq"))
    if(is.null(nam)){
        d <- getwd()
        nam <- d %>% str_remove(".*/")
    }
    write.csv(df, paste0(nam, "_hybridData.csv"), row.names = F)
    setwd(cwd)
    print(dir)
}
hybridMetric <- cmpfun(hybridMetric)


