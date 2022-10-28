
library(RColorBrewer)
library(ggpubr)
library(pheatmap)
library(pathview)
library(graphite)
library(msigdbr)
library(singscore)
library(GSEABase)
library(openxlsx)
library(dplyr)
library(stringr)


#
Dados <- read.xlsx("./Dados.xlsx")




w <- 1
Dados$LampDose <- sapply(1:dim(Dados)[1], function(w){
  paste0(Dados$Wave[w], "_", Dados$Dose[w])
})



ScoreData <- Dados[Dados$Plot == "CPD",]



# Create mean table
df_mean <- ScoreData[,c("Wave", "Time", "Dose", "LampDose", "Value")] %>% 
  group_by(Time, Wave, Dose, LampDose) %>% 
  summarize(average = mean( Value )) %>%
  ungroup()



png(file=paste0("./", "", "", "",
                "CPD", ".png"),
    width=1800,
    height=600,
    res=300)

print(
  
  # Plot
  ggplot(mapping = aes(x = Time, y = Value, fill = LampDose),
         data = ScoreData) + 
    geom_boxplot(alpha=0.6, 
                 outlier.shape = NA) +
    ylab(paste0(unique(ScoreData$Plot))) +
    theme_classic() + 
    geom_jitter(width=0.2,size=0.8,
                alpha=0.3) +
    xlab("Lamp") + 
    
    facet_wrap(~ LampDose,
               ncol = 5) +
    
    geom_point(mapping = aes(x = Time, y = average),
               data = df_mean,
               alpha=0) +
    geom_line(mapping = aes(x = Time, y = average, group=LampDose,
                            color= LampDose),
              data = df_mean
    )   +
    
    scale_fill_manual( values=c("#FFFFFF", "#008435", "#00441B",  "#EA1C26", "#A50F15"    ) ) +
    
    scale_color_manual( values=c("#000000", "#000000", "#000000", "#000000", "#000000") ) +
    
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     color = "black"),
          axis.text.y = element_text(angle = 45,
                                     hjust = 1,
                                     color = "black")
    )
  
)



dev.off()















ScoreData <- Dados[Dados$Plot == "6-4-PP",]



# Create mean table
df_mean <- ScoreData[,c("Wave", "Time", "Dose", "LampDose", "Value")] %>% 
  group_by(Time, Wave, Dose, LampDose) %>% 
  summarize(average = mean( Value )) %>%
  ungroup()



png(file=paste0("./", "", "", "",
                "6-4-PP", ".png"),
    width=1800,
    height=600,
    res=300)

print(
  
  # Plot
  ggplot(mapping = aes(x = Time, y = Value, fill = LampDose),
         data = ScoreData) + 
    geom_boxplot(alpha=0.6, 
                 outlier.shape = NA) +
    ylab(paste0(unique(ScoreData$Plot))) +
    theme_classic() + 
    geom_jitter(width=0.2,size=0.8,
                alpha=0.3) +
    xlab("Lamp") + 
    
    facet_wrap(~ LampDose,
               ncol = 5) +
    
    geom_point(mapping = aes(x = Time, y = average),
               data = df_mean,
               alpha=0) +
    geom_line(mapping = aes(x = Time, y = average, group=LampDose,
                            color= LampDose),
              data = df_mean
    )   +
    
    scale_fill_manual( values=c("#FFFFFF", "#008435", "#00441B",  "#EA1C26", "#A50F15"    ) ) +
    
    scale_color_manual( values=c("#000000", "#000000", "#000000", "#000000", "#000000") ) +
    
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     color = "black"),
          axis.text.y = element_text(angle = 45,
                                     hjust = 1,
                                     color = "black")
    )
  
)



dev.off()










