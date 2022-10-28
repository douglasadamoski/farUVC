


#
library(openxlsx)
library(sp)
library(raster)
library(fields)

#
Dados <- read.xlsx(
  xlsxFile="Dados.xlsx",
  sheet="Planilha1",
  startRow = 1,
  colNames = TRUE,
  rowNames = FALSE,
  detectDates = FALSE,
  skipEmptyRows = TRUE,
  skipEmptyCols = TRUE,
  # rows = NULL,
  cols = 1:7,
  check.names = FALSE,
  sep.names = ".",
  namedRegion = NULL,
  na.strings = "NA",
  fillMergedCells = FALSE
)



unique(Dados$Distancia)

#


# Plot agora
ondaNow <- "222nm"
distanciaNow <- "Box_01"


for(ondaNow in unique(Dados$Onda)){
  for(distanciaNow in unique(Dados$Distancia)){
    
    
    #
    TempTable <- Dados[which(Dados$Distancia == distanciaNow &
                               Dados$Onda == ondaNow),  ]
    if(dim(TempTable)[1] > 0){
      
      
      
      # Adiciona extrapolações da borda
      linhasAgora <- (dim(TempTable)[1] +1 ):(dim(TempTable)[1] +4 )
      
      TempTable[linhasAgora, "Distancia"] <- distanciaNow
      TempTable[linhasAgora, "Onda"] <- ondaNow
      
      
      # x min, y min
      ValuesNow <- c(TempTable[which(TempTable[, "x"] == min(TempTable[, "x"], na.rm=TRUE) &
                        TempTable[, "y"] == min(TempTable[, "y"], na.rm=TRUE)), "mW_cm2"],
      # 
      TempTable[which(TempTable[, "x"] == max(TempTable[, "x"], na.rm=TRUE) &
                        TempTable[, "y"] == min(TempTable[, "y"], na.rm=TRUE)), "mW_cm2"],
      
      # 
      TempTable[which(TempTable[, "x"] == min(TempTable[, "x"], na.rm=TRUE) &
                        TempTable[, "y"] == max(TempTable[, "y"], na.rm=TRUE)), "mW_cm2"],
      
      # 
      TempTable[which(TempTable[, "x"] == max(TempTable[, "x"], na.rm=TRUE) &
                        TempTable[, "y"] == max(TempTable[, "y"], na.rm=TRUE)), "mW_cm2"])

      #
      TempTable[linhasAgora, "x"] <- c(0, 16.8, 0, 16.8)
      TempTable[linhasAgora, "y"] <- c(0, 0, 12.6, 12.6)
      TempTable[linhasAgora, "mW_cm2"] <- ValuesNow
      TempTable[linhasAgora, "Inside.area"] <- "Não"
  
      
      #
      coords <- TempTable[, c("x", "y")]
      s      <- SpatialPointsDataFrame(coords = coords, data = TempTable)
      
      # interpolate with a thin plate spline 
      # (or another interpolation method: kriging, inverse distance weighting). 
      
      tps <- Tps(coordinates(s), as.vector(TempTable$mW_cm2))
      
      p   <- raster(s)
      p   <- interpolate(p, tps)
      
      # Stdize values
      p@data@values <- (( max(TempTable$mW_cm2) - min(TempTable$mW_cm2) ) * ((p@data@values - min(p@data@values))/( max(p@data@values) - min(p@data@values) ))) + min(TempTable$mW_cm2)
      
      #
      # print(paste0("For", "minimum: ", min(p@data@values), "; maximum: ", max(p@data@values)))

      
      # plot raster, points, and contour lines
      brks <- seq(0,
                  1.05*max(TempTable$mW_cm2),
                  by=0.0001) 
      
      #color_now <- heat.colors(length(brks))[length(brks):1]
      color_now <- colorRampPalette(c("royalblue","springgreen","yellow","red"))(length(brks))
      
      ####
      png(file=paste0(ondaNow, "_", distanciaNow, ".png"),
          width=2500,
          height=1600,
          res=300)
      
      plot(p,
           xlim=c(0,130),
           ylim=c(0,85),
           xaxs="i",
           yaxs="i",
           breaks=brks,
           lab.breaks= brks,
           #zlim=c(0,1)
           col=color_now,
           asp=1.03
      )
      
      plot(s, add=T)
      
      dev.off()
      #
      
    }
    
  }
  
}



#





# contour(p, add=T) 




