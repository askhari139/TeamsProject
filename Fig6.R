mainFolder <- "D:/Teams"
plotFolder <- paste0(mainFolder, "/Figures/Fig6")
source(paste0(mainFolder, "/Figures/figCore.R"))
source(paste0(mainFolder, "/Final_Results/codes/stateLabeller.R"))

### AUC data and Plots ----
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
    setwd(groupStrengthCausation)
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
    setwd("..")
    write.csv(df, paste0(net, "_aocDat.csv"), row.names = F)
}
sapply(EMPNets, aocData)

aocPlots <- function(net)
{
    setwd(groupStrengthCausation)
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
        if (!dir.exists("AUCs"))
            dir.create("AUCs")
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

### Sigmoidal fits ----
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
    setwd(groupStrengthCausation)
    
    setwd(net)
    filz <- list.files(".", "allNodeCoherence")
    plan(multisession, workers = 8)
    fittedDat <- future_lapply(filz, sigmFit) %>% reduce(rbind.data.frame)
    future:::ClusterRegistry("stop")
    setwd("..")
    df <- fittedDat
    write_csv(df, paste0(net, "_sigmoidalFits.csv"))
    print(net)
}
sapply(EMPNets, sigmDat)
sigmPlots <- function(net)
{
    setwd(groupStrengthCausation)
    fitDat <- read_csv(paste0(net, "_sigmoidalFits.csv")) %>%
        group_by(Net, Phenotype) %>% 
        summarise(Avg = mean(Cooperativity), Std = sd(Cooperativity))
    
    groupDat <- read_csv(paste0(net, "/", net, "_groups.csv")) %>% 
        select(Net, Gs)
    df <- merge(fitDat, groupDat, by = "Net")
    WTDat <- df %>% filter(Net == net) %>% mutate(Net = "WT")
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
        geom_text(data = annDat, aes(x = x, y = y, label = grob), size = 5) +
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

### Heatmap ----

metricPlots <- function(net, plot = F)
{
    setwd(groupStrengthCausation)
    setwd(net)
    csvFiles <- list.files(".", "_finFlagFreq")
    df <- sapply(csvFiles, function(x){
        print(x)
        df <- read_csv(x, show_col_types = F) %>% filter(!is.na(Avg0), flag == 1)
        if (nrow(df) == 0)
            return(c(NA,NA,NA,NA,NA))
        c(max(df$frust0), min(df$frust0), max(df$coherence0, na.rm = T), min(df$coherence0, na.rm = T),
          df %>% filter(phenotype != "H") %>% select(Avg0) %>% unlist %>% sum)
    }) %>% t %>% data.frame %>% set_names(c("maxFrust", "minFrust", "maxCoh", "minCoh",
                                            "terminalFreq")) %>%
        mutate(Net = str_remove(csvFiles, "_finFlagFreq.csv"))
    write_csv(df, paste0(net, "_Metrics.csv"))
    dfGroup <- read_csv(paste0(net, "_groups.csv")) %>% select(Net, Gs) %>%
        mutate(nEdges = Net %>% str_remove(net) %>% str_count("_"))
    dfSigm <- read_csv(paste0("../",net, "_sigmoidalFits.csv")) %>% 
        select(Net, Cooperativity, Phenotype) %>%
        mutate(Phenotype = ifelse(Phenotype == "H", "Hybrid", "Terminal")) %>%
        group_by(Net, Phenotype) %>% summarise(Cooperativity = mean(Cooperativity))%>%
        spread(key = Phenotype, value = Cooperativity)
    setwd(paste0(plotFolder, "/GsCausation"))
    if (!dir.exists(net))
        dir.create(net)
    setwd(net)
    ggplot(dfGroup, aes(x = nEdges, y = Gs)) + 
        geom_point(size = 2) + geom_line() + 
        theme_Publication() +
        labs(x = "Number of edges deleted", y = "Mean Group Strength")
    ggsave(paste0(net, "_GroupStrengthChange.png"), width = 5.5, height = 5)
    
    df <- merge(df, dfGroup, by = "Net", all = T)
    df <- merge(df, dfSigm, by = "Net", all = T)
    metrics <- c("maxFrust", "minFrust", "maxCoh", "minCoh",
                 "terminalFreq", "Hybrid", "Terminal")
    MetricNames <- c("Maximum\nFrustration", "Minimum\nFrustration", "Maximum\nCoherence",
                     "Minimum\nCoherence", "Terminal State\nFrequency", "Cooperativity\nHybrid",
                     "Cooperativity\nTerminal")
    names(MetricNames) <- metrics
    if(plot)
    {
        sapply(metrics, function(x){
            grob <- correlGrob(df, "Gs", x)
            ggplot(df, aes_string(x = "Gs", y = x)) + 
                geom_point(size = 2) + 
                annotation_custom(grob) +
                theme_Publication() +
                labs(x = "Mean Group Strength", y = MetricNames[x] %>% str_replace("\n", " "))
            ggsave(paste0(net, "_", x, "vsGs.png"), width = 5.5, height = 5)
        })
    }
    
    df <- df %>% select(Gs, all_of(metrics)) %>% gather(key = "Metrics", value = "Value", -Gs) %>%
        mutate(Metrics = MetricNames[Metrics])
    multiFactorCorrelation(df, 'Metrics', "Value", "Gs", label = F) %>% mutate(Net = net)
}
df <- lapply(EMPNets, metricPlots, plot = T) %>% reduce(rbind.data.frame) %>%
    mutate(Label = ifelse(pValue < 0.05, "", "X")) %>%
    mutate(Network = netNameKey[Net]) %>%
    mutate(Metrics = Factors)

ggplot(df, aes(y = Network, x = Factors, fill = Correlation)) +
    geom_tile() +
    geom_text(aes(label = Label))+
    scale_fill_gradient2(low = "blue", high = "red", limits = c(-1,1)) +
    theme_Publication() + 
    theme(axis.title = element_blank(), legend.position = "right", 
          legend.direction = "vertical", legend.key.height = unit(0.8, "cm"), 
          axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 1)) +
    labs(x = "", y = "")
ggsave(paste0(plotFolder, "/GsCausation/allNetCorrs.png"), width = 10, height = 6)
