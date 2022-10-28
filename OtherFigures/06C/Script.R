
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



ScoreData <- Dados


# Create mean table
df_mean <- ScoreData[,c("Wave", "Time", "Dose", "LampDose", "Value")] %>% 
  group_by(Time, Wave, Dose, LampDose) %>% 
  summarize(average = mean( Value )) %>%
  ungroup()



png(file=paste0("./", "", "", "",
                "DCF", ".png"),
    width=2840,
    height=1100,
    res=400)

print(
  
  # Plot
  ggplot(mapping = aes(x = LampDose, y = Value, fill = LampDose),
         data = ScoreData) + 
    geom_boxplot(alpha=0.6, 
                 outlier.shape = NA) +
    ylab(paste0(unique(ScoreData$Plot))) +
    theme_classic() + 
    geom_jitter(width=0.2,size=1.2,
                alpha=0.4) +
    xlab("Lamp") + 
    
    #facet_wrap(~ LampDose,
    #           ncol = 7) +
    
    geom_point(mapping = aes(x = LampDose, y = average),
               data = df_mean,
               alpha=0) +
    scale_fill_manual( values=c("#FFFFFF", "#00BF4D", "#008435", "#00441B",  "#EE4851", "#EA1C26", "#A50F15"    ) ) +
    
    scale_color_manual( values=c("#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000") ) +
    
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     color = "black"),
          axis.text.y = element_text(angle = 45,
                                     hjust = 1,
                                     color = "black")
    )
  
)



dev.off()








