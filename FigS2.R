mainFolder <- "D:/Teams"
plotFolder <- paste0(mainFolder, "/Figures/FigS2")
source(paste0(mainFolder, "/Figures/figCore.R"))


### Quantitative Convergence ---

QCPlotterBool <- function(){
    setwd(plotFolder)
    if(!dir.exists("QCplots"))
        dir.create("QCplots")
    
    setwd(QC)
    dirs <- list.dirs(".", recursive = F) %>% str_remove("./")
    dfList <- lapply(EMPNets, function(x){
        freqFiles <- paste0(dirs, "/", x, "_finFlagFreq.csv")
        df <- lapply(freqFiles, function(y){
            n <- str_extract(y, ".*/") %>% str_remove("/")
            read_csv(y, col_types = cols()) %>% filter(!is.na(Avg0), flag == 1) %>%
                select(states, Avg0, SD0) %>% set_names("states", paste0(n, c("_Avg", "_Std")))
        }) %>% reduce(merge, by = "states", all = T)
        df[is.na(df)] <- 0
        dfAvg <- df %>% select(states, contains("Avg")) %>% 
            gather(key = "nInit", value = "Frequency", -states) %>% 
            mutate(nInit = str_remove(nInit, "_Avg"))
        dfStd <- df %>% select(states, contains("Std")) %>% 
            gather(key = "nInit", value = "Std", -states) %>% 
            mutate(nInit = str_remove(nInit, "_Std"))
        df <- merge(dfAvg, dfStd, by = c("states", "nInit"), all = T) %>% mutate(Net = netNameKey[x])
        df$CV <- df$Std/(df$Frequency + 1e-5)
        write.csv(df, paste0(x, "_QCDat.csv"), row.names = F)
        setwd(paste0(plotFolder, "/QCPlots"))
        ggplot(df %>% filter(Frequency > 0.005*max(Frequency)), 
               aes(x = reorder(states, -Frequency), y = Frequency, color = nInit, shape = nInit)) + 
            geom_point(size = 2.5, position = position_dodge2(width = 0.9)) + 
            geom_errorbar(aes(ymin = Frequency - Std, ymax = Frequency + Std), 
                          position = position_dodge2(width = 0.9)) +
            theme_Publication() + 
            theme(axis.text.x = element_blank(),
                  legend.position = c(0.8, 0.8), legend.direction = "vertical",
                  legend.text = element_text(size = rel(1)), 
                  legend.title = element_text(size = rel(1))) +
            facet_wrap(~Net) +
            labs(x = "States", color = "# Initial Conditions", shape = "# Initial Conditions")
        ggsave(paste0(x, "_Bool_QCAvgStd.png"), width = 5.5, height = 5)
        
        ggplot(df %>% filter(Frequency > 0.005*max(Frequency)), 
               aes(x = reorder(states, -Frequency), y = CV, color = nInit, shape = nInit, 
                   group = nInit)) + 
            geom_point(size = 2.5, position = position_dodge2(width = 0.9)) + 
            geom_line()+
            scale_y_log10() + 
            theme_Publication() + 
            theme(axis.text.x = element_blank(),
                  legend.position = "top") +
            facet_wrap(~Net) +
            labs(x = "States", y = "Std/Mean",
                 color = "# Initial Conditions", shape = "# Initial Conditions")
        ggsave(paste0(x, "_Bool_QCCV.png"), width = 5.5, height = 5)
        
        ggplot(df , 
               aes(x = CV, fill = nInit)) + 
            geom_density(alpha =  0.5)+
            theme_Publication() + 
            theme(legend.position = c(0.8, 0.8), legend.direction = "vertical", 
                  legend.text = element_text(size = rel(1)), 
                  legend.title = element_text(size = rel(1))) +
            facet_wrap(~Net) +xlim(c(-0.05, 1)) +
            labs(x = "Std/Mean", y = "Density", fill = "# Initial Conditions")
        ggsave(paste0(x, "_Bool_QCCVdist.png"), width = 5.5, height = 5)
        
        setwd(QC)
    })
}
QCPlotterBool()

QCPlotterRAC <- function()
{
    setwd(plotFolder)
    if(!dir.exists("QCplots"))
        dir.create("QCplots")
    
    setwd(RACIPE_QC)
    dirs <- list.dirs(".", recursive = F) %>% str_remove("./")
    nPars <- c("1000", "10000", "100000")
    dfList <- lapply(racEMPNets, function(x){
        df <- lapply(nPars, function(n){
            y <- paste0(x, "/", x, "_", n, ".csv")
            read_csv(y, col_types = cols()) %>% set_names("states", paste0(n, c("_Avg", "_Std")))
        }) %>% reduce(merge, by = "states", all = T)
        df[is.na(df)] <- 0
        dfAvg <- df %>% select(states, contains("Avg")) %>% 
            gather(key = "nInit", value = "Frequency", -states) %>% 
            mutate(nInit = str_remove(nInit, "_Avg"))
        dfStd <- df %>% select(states, contains("Std")) %>% 
            gather(key = "nInit", value = "Std", -states) %>% 
            mutate(nInit = str_remove(nInit, "_Std"))
        df <- merge(dfAvg, dfStd, by = c("states", "nInit"), all = T) %>% 
            mutate(Net = netNameKeyRAC[x])
        df$CV <- df$Std/(df$Frequency + 1e-5)
        write.csv(df, paste0(x, "_QCDat.csv"), row.names = F)
        setwd(paste0(plotFolder, "/QCPlots"))
        ggplot(df %>% filter(Frequency > 0.005*max(Frequency)), 
               aes(x = reorder(states, -Frequency), y = Frequency, color = nInit, shape = nInit)) + 
            geom_point(size = 2.5, position = position_dodge2(width = 0.9)) + 
            geom_errorbar(aes(ymin = Frequency - Std, ymax = Frequency + Std), 
                          position = position_dodge2(width = 0.9)) +
            theme_Publication() + 
            theme(axis.text.x = element_blank(),
                  legend.position = c(0.8, 0.8), legend.direction = "vertical",
                  legend.text = element_text(size = rel(1)),
                  legend.title = element_text(size = rel(1))) +
            facet_wrap(~Net) +
            labs(x = "States", color = "# Parameter Sets", shape = "# Parameter Sets")
            
        ggsave(paste0(x, "_RAC_QCAvgStd.png"), width = 5.5, height = 5)
        
        ggplot(df %>% filter(Frequency > 0.005*max(Frequency)), 
               aes(x = reorder(states, -Frequency), y = CV, color = nInit, shape = nInit, 
                   group = nInit)) + 
            geom_point(size = 2.5, position = position_dodge2(width = 0.9)) + 
            geom_line()+
            scale_y_log10() + 
            theme_Publication() + 
            theme(axis.text.x = element_blank(),
                  legend.position = c(0.8, 0.8), legend.direction = "vertical",
                  legend.text = element_text(size = rel(1)), 
                  legend.title = element_text(size = rel(1))) +
            facet_wrap(~Net) +
            labs(x = "States", y = "Std/Mean",
                 color = "# Parameter Sets", shape = "# Parameter Sets")
        ggsave(paste0(x, "_RAC_QCCV.png"), width = 5.5, height = 5)
        
        ggplot(df , 
               aes(x = CV, fill = nInit)) + 
            geom_density(alpha =  0.5)+
            theme_Publication() + 
            theme(legend.position = c(0.8, 0.8), legend.direction = "vertical",
                  legend.text = element_text(size = rel(1)), 
                  legend.title = element_text(size = rel(1))) +
            facet_wrap(~Net) +xlim(c(-0.05, 1)) +
            labs(x = "Std/Mean", y = "Density", fill = "# Parameter Sets")
        ggsave(paste0(x, "_RAC_QCCVdist.png"), width = 5.5, height = 5)
        
        setwd(RACIPE_QC)
    })
}
QCPlotterRAC()
