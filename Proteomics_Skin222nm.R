
## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

##### Library loading #####
#https://datascienceplus.com/proteomics-data-analysis-1-3-data-acquisition-and-cleaning/#:~:text=In%20our%20case%2C%20Q.,values%20are%20below%20the%20threshold.
library(data.table)
library(ggplot2)
library(ggfortify)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(ggpmisc)
library(gridExtra)
library(dplyr)
library(httr)
library(jsonlite)
# Bioconductor packages
library(limma)
library(vsn)
library(affy)
library(biomaRt)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)

#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install(c("limma", "vsn", "affy", "biomaRt",
#                       "ReactomePA", "clusterProfiler", "org.Hs.eg.db"))

# Read proteinGroup table
DadosBrutos <- read.table(file=paste0("proteinGroups.txt"),
                          header=TRUE,
                          stringsAsFactors = FALSE,
                          check.names = FALSE,
                          sep="\t")

# https://pubmed.ncbi.nlm.nih.gov/24942700/
# Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ
# LFQ
#
# A systematic evaluation of normalization methods in quantitative label-free proteomics 
# https://academic.oup.com/bib/article/19/1/1/2562889
#

##### Change mislabelled sample #####
colnames(DadosBrutos) <- gsub("Lampada_2_48h_3", "TROCANDONOME", colnames(DadosBrutos))
colnames(DadosBrutos) <- gsub("Controle_neg_24h_3", "Lampada_2_48h_3", colnames(DadosBrutos))
colnames(DadosBrutos) <- gsub("TROCANDONOME", "Controle_neg_24h_3", colnames(DadosBrutos))

#
nameModifier <- "_00_RawData"

####
AmostrasDisponiveis <- gsub("Identification type ", "", colnames(DadosBrutos)[104:130])

#
dir.create(path = paste0("./", "Results", nameModifier, "/"), showWarnings = FALSE)

##### General description #####
## Sequence Coverage

png(paste0("./", "Results", nameModifier, "/",
           "01_CoberturaMedia.png"),
    width=1500, height=1000, res=300)
hist(rowMeans(DadosBrutos[, paste0("Sequence coverage ", AmostrasDisponiveis, " [%]")]),
     breaks=200,
     main="Protein-wise sequence Coverage (%)",
     xlab="Mean sequence coverage (%)",
     ylab="Amount",
     #ylim=c(0,1.5),
     xlim=c(0, 100),
     freq=FALSE,
     border="white", 
     col=rgb(1, 0, 0, 1/4))
dev.off()



#### Melted coverage
MeltedCoverage <- data.table::melt(setDT(DadosBrutos[, c("Majority protein IDs", paste0("Sequence coverage ", AmostrasDisponiveis, " [%]"))]),
                                   id.vars = c("Majority protein IDs"),
                                   variable.name = "SequenceCoverage")

MeltedCoverage <- data.frame(MeltedCoverage, stringsAsFactors = FALSE, check.names=FALSE)
MeltedCoverage$SequenceCoverage <- gsub(" \\[%\\]", "", MeltedCoverage$SequenceCoverage)
MeltedCoverage$SequenceCoverage <- gsub("Sequence coverage ", "", MeltedCoverage$SequenceCoverage)
#

p <- ggplot(data=MeltedCoverage,
            aes(x=SequenceCoverage,
                y=value)) + 
  geom_boxplot() +
  theme_classic() +
  ylab("Sequence Coverage (%)") +
  geom_jitter(shape=16,
              position=position_jitter(0.2),
              alpha=0.01) +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))


png(paste0("./", "Results", nameModifier, "/",
           "02_CoberturaPorAmostra.png"),
    width=1800, height=1200, res=300)
print(p)
dev.off()

##### Filters ######
## Remove potential contaminants
DadosBrutos[,"QValue_fail"] <- ""
DadosBrutos[which(DadosBrutos[,"Q-value"] > 0.01),"QValue_fail"] <- "+"

plotDF <- data.frame(table(DadosBrutos[, c("Potential contaminant",
                                 "Reverse",
                                 "Only identified by site") ]))

colnames(plotDF) <- c("Potential contaminant", "Reverse", "Only identified by site", "Proteins" )

#
png(paste0("./", "Results", nameModifier, "/",
           "03_ProteinQuality.png"), height=900, width=1000,
    res=120)
p<-tableGrob(plotDF)
grid.arrange(p)
dev.off()


#
# for data manipulation

# Filter false hits
DadosBrutos[which(DadosBrutos[,"Potential contaminant"] != "+" &
                    DadosBrutos[,"Reverse"] != "+" & 
                    DadosBrutos[,"Only identified by site"] != "+" & 
                    DadosBrutos[,"QValue_fail"] != "+"
)  , ] -> DadosBrutos


#



# Only LFQ table

####
DadosLFQ <- DadosBrutos[, c("Protein IDs",
                            "Majority protein IDs",
                            "Fasta headers",
                            "Gene names",
                            "Protein names",
                            "Number of proteins",
                            paste0("LFQ intensity ", AmostrasDisponiveis))]


colnames(DadosLFQ) <- gsub("LFQ intensity ", "", colnames(DadosLFQ))

#


png(paste0("./", "Results", nameModifier, "/",
           "04_TotalLFQCounts_perSample.png"),
    width=2000, height=1500, res=200)

#bottom, left, top and right
# par(mar=c(5.1,20.1,4.1,2.1))
print(
  ggbarplot(data.frame(Condition=AmostrasDisponiveis,
                       values=as.numeric(colSums(DadosLFQ[,AmostrasDisponiveis])) ),
            x = "Condition",
            y = "values", 
            add = c("mean_se", "jitter")) +
    rotate_x_text(90) +
    
    stat_compare_means() # comparisons=myComparisons
)
dev.off()
####


png(paste0("./", "Results", nameModifier, "/",
           "05_TotalLFQCounts_summary.png"),
    width=1500, height=1500, res=300)

#bottom, left, top and right
# par(mar=c(5.1,20.1,4.1,2.1))
print(
  ggbarplot(data.frame(Condition=gsub("h_3", "h", gsub("h_2", "h", gsub("h_1", "h", AmostrasDisponiveis))),
                       values=as.numeric(colSums(DadosLFQ[,AmostrasDisponiveis])) ),
            x = "Condition",
            y = "values", 
            ylab = "Total LFQ Counts",
            add = c("mean_se", "jitter")) +
    rotate_x_text(90) +
    
    stat_compare_means() # comparisons=myComparisons
)
dev.off()
####




#
MeltedLFQ <- data.table::melt(setDT(DadosLFQ[,c("Majority protein IDs", AmostrasDisponiveis)]),
                              id.vars = c("Majority protein IDs")
)


##
png(paste0("./", "Results", nameModifier, "/",
           "06_TotalLFCCounts_histogram.png"),
    width=1500, height=800, res=300)

print(
  ggplot(MeltedLFQ, aes(x=value) ) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.2, fill="#FF6666") +
    theme_classic() +
    ylab("Density") +
    xlab("LFQ")
  
)

dev.off()


##
png(paste0("./", "Results", nameModifier, "/",
           "07_TotalLogCounts_histogram.png"),
    width=1500, height=800, res=300)

print(
  ggplot(MeltedLFQ, aes(x=log2(value))) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.2, fill="#FF6666") +
    theme_classic() +
    ylab("Density") +
    xlab("log2(LFQ)")
  
)

dev.off()


##



### REMOVE EMPTY GENE NAMES
### SOLVE DUPLICATES WITH SUMS
# Extract IDs
# Isolate the first entry
DadosLFQ[,"UniprotID"] <- sub(";.*", "", DadosLFQ[,"Protein IDs"])
DadosLFQ[,"Fasta headers"] = sub(";.*", "", DadosLFQ[,"Fasta headers"])

#Remove isoform
DadosLFQ[ ,"UniprotID"] <- sapply(1:dim(DadosLFQ)[1], function(w){
  strsplit(DadosLFQ[w,"UniprotID"], "-")[[1]][1]
})


# Extract Protein name
regex = regexpr("(?<=_HUMAN.).*(?=.OS)", DadosLFQ[,"Fasta headers"], perl = TRUE)
DadosLFQ$Protein.name = regmatches(DadosLFQ[,"Fasta headers"], regex)

# Extract UniProtID
#regex = regexpr("(?<=\\|).*(?=\\|)", df$Protein.IDs, perl = TRUE)
#df$Protein = regmatches(df$Protein.IDs, regex)

# Extract Gene ID
#regex = regexpr("((?<=\\|[[:alnum:]]{6}\\|).*(?=_HUMAN)|(?<=\\|[[:alnum:]]{10}\\|).*(?=_HUMAN))",
#                df$Protein.IDs, perl = TRUE)
#df$Gene = regmatches(df$Protein.IDs, regex)

UniprotLocalBaseJSON <- "D:/OneDrive/Databases/Uniprot/UniprotIDs_tsv/"

#
for(idNow in unique(DadosLFQ$UniprotID) ){
  if(file.exists( paste(UniprotLocalBaseJSON, idNow, ".tsv", sep="") )){
    if(file.info(paste(UniprotLocalBaseJSON, idNow, ".tsv", sep="") )$size>50  ){
      DownloadFile <- FALSE
    } else {
      DownloadFile <- TRUE
    }
  } else {
    DownloadFile <- TRUE
  }
  
  if(DownloadFile){
    # Create URL 
    URL_now <- paste0("https://rest.uniprot.org/uniprotkb/", idNow, "?format=tsv&fields=accession,id,protein_name,gene_names,xref_ensembl,xref_geneid")
    # JSON download
    writeBin(content(GET(URL_now, progress()),
                     "text"), paste(UniprotLocalBaseJSON, idNow, ".tsv", sep=""))
  }
}

##### read uniprot tables #####
if(exists("UniprotUpdateTable")){
  remove(UniprotUpdateTable)
}
myColNames <- c("Entry", "Entry Name", "Protein names", "Gene Names", "Ensembl", "GeneID")

for(idNow in unique(DadosLFQ$UniprotID) ){
  tempTable <-  read.table(file=paste(UniprotLocalBaseJSON, idNow, ".tsv", sep=""),
                           sep="\t", header=TRUE, stringsAsFactors = FALSE,
                           check.names = FALSE,
                           comment.char = "#")
  if(sum(myColNames %in% colnames(tempTable)) != 6){
    tempTable[,myColNames[!(myColNames %in% colnames(tempTable))]] <- ""
  }
  if(exists("UniprotUpdateTable")){
    UniprotUpdateTable <- rbind(UniprotUpdateTable, tempTable[,myColNames])
  } else {
    UniprotUpdateTable <- tempTable[,myColNames]
  }
}


#
DadosLFQ[,"GeneSymbol_final"] <- sapply(UniprotUpdateTable[match(DadosLFQ[,"UniprotID"], UniprotUpdateTable[,"Entry"]), "Gene Names"],
                                        function(w){
                                          strsplit(w, " ")[[1]][1]
                                        })
DadosLFQ[,"GeneID_final"] <- sapply(UniprotUpdateTable[match(DadosLFQ[,"UniprotID"], UniprotUpdateTable[,"Entry"]), "GeneID"],
                                    function(w){
                                      strsplit(w, " ")[[1]][1]
                                    })

#
AnotationDBI_table <- AnnotationDbi::select(org.Hs.eg.db,
                                            strsplit(paste(DadosLFQ[, "Protein IDs"], collapse=";"), ";")[[1]],
                                            keytype = "UNIPROT",
                                            columns = c("SYMBOL",
                                                        "REFSEQ",
                                                        "ENTREZID") )
#
for(lineNow in which(is.na(DadosLFQ[,"GeneSymbol_final"])) ){
  #
  myTemp <- AnotationDBI_table[match(strsplit(DadosLFQ[lineNow, "Protein IDs"], ";")[[1]],
                                     AnotationDBI_table$UNIPROT),]
  if(dim(myTemp)[1]>0){
    DadosLFQ[lineNow, "GeneSymbol_final"] <- myTemp$SYMBOL[!is.na(myTemp$SYMBOL)][1]
    DadosLFQ[lineNow, "GeneID_final"] <- myTemp$ENTREZID[!is.na(myTemp$ENTREZID)][1]
  }
}

#
for(lineNow in which(is.na(DadosLFQ[,"GeneID_final"])) ){
  #
  myTemp <- AnotationDBI_table[match(strsplit(DadosLFQ[lineNow, "Protein IDs"], ";")[[1]],
                                     AnotationDBI_table$UNIPROT),]
  #
  if(dim(myTemp)[1]>0){
    DadosLFQ[lineNow, "GeneSymbol_final"] <- myTemp$SYMBOL[!is.na(myTemp$SYMBOL)][1]
    DadosLFQ[lineNow, "GeneID_final"] <- myTemp$ENTREZID[!is.na(myTemp$ENTREZID)][1]
  }
}

# Remove ;
for(lineNow in 1:dim(DadosLFQ)[1] ){
  DadosLFQ$GeneID_final[lineNow] <- strsplit(DadosLFQ$GeneID_final[lineNow], ";")[[1]][1]
}
for(lineNow in 1:dim(DadosLFQ)[1] ){
  DadosLFQ$GeneSymbol_final[lineNow] <- strsplit(DadosLFQ$GeneSymbol_final[lineNow], ";")[[1]][1]
}

# Add values to main table of proteins
# Remove empty gene names
DadosLFQ <- DadosLFQ[!(DadosLFQ[,"GeneID_final"] == ""),]
DadosLFQ <- DadosLFQ[!( is.na(DadosLFQ[,"GeneID_final"]) ),]

#
AnotationDBI_table <- AnnotationDbi::select(org.Hs.eg.db,
                                            strsplit(paste(DadosLFQ[, "GeneID_final"], collapse=";"), ";")[[1]],
                                            keytype = "ENTREZID",
                                            columns = c("SYMBOL",
                                                        "REFSEQ",
                                                        "ENTREZID") )


for(lineNow in 1:dim(DadosLFQ)[1] ){
  mySYMBOLNow <- AnotationDBI_table[match(DadosLFQ[lineNow, "GeneID_final"], AnotationDBI_table$ENTREZID), "SYMBOL"]
  
  mySYMBOLNow -> DadosLFQ[lineNow, "GeneSymbol_final"]
  
}


# Remove empty or NA geneSymbol
DadosLFQ <- DadosLFQ[!(DadosLFQ[,"GeneSymbol_final"] == ""),]
DadosLFQ <- DadosLFQ[!( is.na(DadosLFQ[,"GeneSymbol_final"]) ),]

# Resumo final por nomenclatura de proteina
dim(DadosLFQ)

DadosLFQ$Keep <- TRUE

### Duplicate removal sum values from same gene id
for(geneNow in  unique(DadosLFQ[, "GeneSymbol_final"])  ){
  #
  if( sum(geneNow == DadosLFQ[, "GeneSymbol_final"]) > 1 ){
    #
    linesNow <- which((geneNow == DadosLFQ[, "GeneSymbol_final"]))
    #
    DadosLFQ$Keep[linesNow] <- FALSE
    #
    DadosLFQ$Keep[linesNow[(which.min( rowSums(DadosLFQ[linesNow, AmostrasDisponiveis] == 0)  ))]] <- TRUE
  }
}

# Remove lines with zeroes
DadosLFQ <- DadosLFQ[DadosLFQ$Keep, ]

# Resumo final por nomenclatura de proteina
dim(DadosLFQ)

#
colnames(DadosLFQ)[which(colnames(DadosLFQ) == "Fasta headers")] <- "Fasta.headers"


#
write.csv2(
  DadosLFQ[,c("GeneID_final", "GeneSymbol_final", "Protein.name", "Fasta.headers", AmostrasDisponiveis)],
  file=paste0("./", "Results", nameModifier, "/",
              "08_SelectedRawData.csv")
)





##### Imputation section #####
library(data.table)
library(ggplot2)
library(ggfortify)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
# Bioconductor packages
library(biomaRt)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pcaMethods)
library(preprocessCore)
library(randomcoloR)
library(EBSeq)
library(limma)

# Read proteinGroup table
DadosLFQ <- read.csv2(file=paste0("./", "Results", "_00_RawData", "/",
                                  "08_SelectedRawData.csv"),
                      header=TRUE,
                      stringsAsFactors = FALSE,
                      check.names = FALSE,
                      row.names=1)

#
AmostrasDisponiveis <- colnames(DadosLFQ)[ !(colnames(DadosLFQ) %in% c("GeneID_final", "GeneSymbol_final", "Protein.name", "Fasta.headers")) ]

#
nameModifier <- "_01_TodosOsGrupos"

#
dir.create(path = paste0("./", "Results", nameModifier, "/"), showWarnings = FALSE)

# PIPELINE:
# 
# 1) Log-transform and Exchange all zeroes for NA
# 2) Remove genes without at least TWO replicates valid in a single group
# 3) Impute
# 4) Normalize
# 5) DE
#


##### Sample info #####
GroupTable <- read.csv2(file="./Agrupamento_01.csv",
                        header=TRUE,
                        stringsAsFactors = FALSE,
                        check.names = FALSE
)

GroupTable -> GroupTable_factored

GroupTable_factored$Sample <- factor(GroupTable_factored$Sample, levels=unique(GroupTable_factored$Sample))
GroupTable_factored$Lamp <- factor(GroupTable_factored$Lamp, levels=unique(GroupTable_factored$Lamp))
GroupTable_factored$Time <- factor(GroupTable_factored$Time, levels=unique(GroupTable_factored$Time))
GroupTable_factored$LampTime <- factor(GroupTable_factored$LampTime, levels=unique(GroupTable_factored$LampTime))

#  Create all comparisons
GroupTable_summary <- unique(GroupTable[,c("Lamp", "Time", "LampTime")])


ComparisonsTable <- data.frame(matrix(ncol=12,
                                      nrow=0))

colnames(ComparisonsTable) <- c("LampTime_1",
                                "LampTime_2",
                                "LampTime_test_a",
                                "LampTime_test_b",
                                
                                "Time_1",
                                "Time_2",
                                "Time_test_a",
                                "Time_test_b",
                                
                                "Lamp_1",
                                "Lamp_2",
                                "Lamp_test_a",
                                "Lamp_test_b"
                                
)

for(i in 1:dim(GroupTable_summary)[1]){
  for(j in 1:dim(GroupTable_summary)[1]){
    
    if(j>i){
      
      ComparisonsTable[(dim(ComparisonsTable)[1]+1), ] <- c(GroupTable_summary[i,"LampTime"],
                                                            GroupTable_summary[j,"LampTime"],
                                                            paste(GroupTable_summary[i,"LampTime"],
                                                                  GroupTable_summary[j,"LampTime"], sep="-"),
                                                            paste(GroupTable_summary[j,"LampTime"],
                                                                  GroupTable_summary[i,"LampTime"], sep="-"),
                                                            
                                                            GroupTable_summary[i,"Time"],
                                                            GroupTable_summary[j,"Time"],
                                                            
                                                            paste(GroupTable_summary[i,"Time"],
                                                                  GroupTable_summary[j,"Time"], sep="-"),
                                                            paste(GroupTable_summary[j,"Time"],
                                                                  GroupTable_summary[i,"Time"], sep="-"),
                                                            
                                                            GroupTable_summary[i,"Lamp"],
                                                            GroupTable_summary[j,"Lamp"],
                                                            
                                                            paste(GroupTable_summary[i,"Lamp"],
                                                                  GroupTable_summary[j,"Lamp"], sep="-"),
                                                            paste(GroupTable_summary[j,"Lamp"],
                                                                  GroupTable_summary[i,"Lamp"], sep="-")
                                                            
      )
      
    }
    
  }
}

##### Sample processing #####
# Replace zero by NA
for(colNow in AmostrasDisponiveis){
  #
  DadosLFQ[which(DadosLFQ[, colNow] == 0), colNow] <- NA
}

# Make log table
DadosLFQ -> DadosLFQ_log2
log2(DadosLFQ_log2[,AmostrasDisponiveis]) -> DadosLFQ_log2[,AmostrasDisponiveis]


### Check if there is any protein with only 1 sample per group
DadosLFQ_log2[,"Keep"] <- FALSE
for(i in 1:dim(DadosLFQ_log2)[1]){
  print(i)
  NAsNow <- vector()
  for(j in unique(GroupTable_factored$LampTime)){
    NAsNow <- c(NAsNow, sum(is.na(DadosLFQ_log2[i, as.character(GroupTable_factored[which(GroupTable_factored$LampTime == j), "Sample"])  ])))
  }
  if(sum(NAsNow < 2) > 1){
    DadosLFQ_log2[i,"Keep"] <- TRUE
  } else {
  }
}

# Retain only positive
DadosLFQ_log2 <- DadosLFQ_log2[DadosLFQ_log2[, "Keep"],  ]

# Perform imputation
dir.create(path = paste0("./", "Results", nameModifier, "/",
                         "00_ImputeInfo/"), showWarnings = FALSE)
#
# Combine plot pre-Impute
png(file=paste0("./", "Results", nameModifier, "/",
                "00_ImputeInfo/",
                "01_Condition_legend", ".png"),
    width=1000, height=1000, res=200)
display.brewer.pal(n=9, name="Set1")
dev.off()
sink(paste0("./", "Results", nameModifier, "/",
            "00_ImputeInfo/",
            "01_Condition_order", ".txt"))
print(unique(GroupTable$LampTime))
sink()


k <- 1
# Combine plot pre-Impute
png(file=paste0("./", "Results", nameModifier, "/",
                "00_ImputeInfo/",
                "01_Condition_preImpute_All", ".png"),
    width=1000, height=1000, res=200)
unique(GroupTable$LampTime)[1] -> j
tempTable <- DadosLFQ_log2[, GroupTable[which(GroupTable$LampTime == j), "Sample"]  ]

plot(sort(rowMeans(tempTable, na.rm=TRUE)),
     xlab="Protein by abundance index",
     ylab="Protein log-abundance (without zeroes)",
     cex=0.5, pch=16,
     col=add.alpha(col=brewer.pal(n=9, name="Set1"), alpha=0.2)[k],
     ylim=c(8,35),
     xlim=c(0,2200))


for(j in unique(GroupTable$LampTime)[2:length(unique(GroupTable$LampTime))]){
  k <- k+1
  
  tempTable <- DadosLFQ_log2[, GroupTable[which(GroupTable$LampTime == j), "Sample"]  ]
  points(sort(rowMeans(tempTable, na.rm=TRUE)),
         cex=0.5, pch=16,
         col=add.alpha(col=brewer.pal(n=9, name="Set1"), alpha=0.2)[k])
}

abline(v=dim(tempTable)[1]*0.1, lwd=2, col="darkgray")
dev.off()



row.names(DadosLFQ_log2) <- DadosLFQ_log2$GeneSymbol_final
DadosLFQ_log2 -> DadosLFQ_log2_NonImputed


for(j in unique(GroupTable$LampTime)[1:length(unique(GroupTable$LampTime))]){
  # Combine plot pre-Impute
  png(file=paste0("./", "Results", nameModifier, "/",
                  "00_ImputeInfo/",
                  "02_Condition_preImpute_", j, ".png"),
      width=1000, height=1000, res=200)
  
  tempTable <- DadosLFQ_log2[, GroupTable[which(GroupTable$LampTime == j), "Sample"]  ]
  
  plot(sort(rowMeans(tempTable, na.rm=TRUE)),
       xlab="Protein by abundance index",
       ylab="Protein log-abundance (without zeroes)",
       cex=0.5, pch=16,
       col=add.alpha(col=brewer.pal(n=9, name="Set1"), alpha=0.2)[k],
       ylim=c(8,35),
       xlim=c(0,2200))
  abline(v=dim(tempTable)[1]*0.1, lwd=2, col="darkgray")
  dev.off()
  
}

##### imputation per-se #####
# Imputation for low abundant samples by Minimum replacement
for(j in unique(GroupTable$LampTime)){
  
  tempTable <- DadosLFQ_log2[, GroupTable[which(GroupTable$LampTime == j), "Sample"]  ]
  
  PerformMinimumVector <- rowMeans(tempTable, na.rm=TRUE) < quantile(rowMeans(tempTable, na.rm=TRUE), 0.05, na.rm=TRUE)
  
  #
  PerformMinimumVector[is.na(PerformMinimumVector)] <- TRUE
  
  #
  tempTable -> tempTableMinimumImputation
  for(w in 1:dim(tempTableMinimumImputation)[2]){
    tempTableMinimumImputation[which(is.na(tempTableMinimumImputation[, w])), w] <-
      min(tempTableMinimumImputation[, w], na.rm=TRUE)
  }
  
  #
  tempTable[which(rowSums(is.na(tempTable)) != dim(tempTable)[2]), ] -> tempTable_bPCAImputation
  #
  tempTable_bPCAImputation_PCA <- pca(as.matrix(tempTable_bPCAImputation),
                                      method="bpca",
                                      nPcs=2,
                                      eps = 1e-20,
                                      threshold = 1e-10,
                                      maxSteps = 10000)
  tempTable -> tempTable_bPCAImputation
  tempTable_bPCAImputation[rownames(completeObs(tempTable_bPCAImputation_PCA)), ] <- completeObs(tempTable_bPCAImputation_PCA)
  
  ## Perform imputation accordingly the method
  PerformMinimumVector
  tempTable -> tempTable111
  tempTable[PerformMinimumVector, ] <- tempTableMinimumImputation[PerformMinimumVector, ]
  tempTable[!PerformMinimumVector, ] <- tempTable_bPCAImputation[!PerformMinimumVector, ]
  
  tempTable -> DadosLFQ_log2[, GroupTable[which(GroupTable$LampTime == j), "Sample"]  ] 
}

#
k <- 1
# Combine plot POST-Impute
png(file=paste0("./", "Results", nameModifier, "/",
                "00_ImputeInfo/",
                "03_Condition_afterImpute_All", ".png"),
    width=1000, height=1000, res=200)
unique(GroupTable$LampTime)[1] -> j
tempTable <- DadosLFQ_log2[, GroupTable[which(GroupTable$LampTime == j), "Sample"]  ]

plot(sort(rowMeans(tempTable, na.rm=TRUE)),
     xlab="Protein by abundance index",
     ylab="Protein log-abundance (without zeroes)",
     cex=0.5, pch=16,
     col=add.alpha(col=brewer.pal(n=9, name="Set1"), alpha=0.2)[k],
     ylim=c(8,35),
     xlim=c(0,2200))


for(j in unique(GroupTable$LampTime)[2:length(unique(GroupTable$LampTime))]){
  k <- k+1
  
  tempTable <- DadosLFQ_log2[, GroupTable[which(GroupTable$LampTime == j), "Sample"]  ]
  points(sort(rowMeans(tempTable, na.rm=TRUE)),
         cex=0.5, pch=16,
         col=add.alpha(col=brewer.pal(n=9, name="Set1"), alpha=0.2)[k])
}

abline(v=dim(tempTable)[1]*0.1, lwd=2, col="darkgray")
dev.off()




for(j in unique(GroupTable$LampTime)[1:length(unique(GroupTable$LampTime))]){
  # Combine plot pre-Impute
  png(file=paste0("./", "Results", nameModifier, "/",
                  "00_ImputeInfo/",
                  "04_Condition_postImpute_", j, ".png"),
      width=1000, height=1000, res=200)
  
  tempTable <- DadosLFQ_log2[, GroupTable[which(GroupTable$LampTime == j), "Sample"]  ]
  
  plot(sort(rowMeans(tempTable, na.rm=TRUE)),
       xlab="Protein by abundance index",
       ylab="Protein log-abundance (without zeroes)",
       cex=0.5, pch=16,
       col=add.alpha(col=brewer.pal(n=9, name="Set1"), alpha=0.2)[k],
       ylim=c(8,35),
       xlim=c(0,2200))
  abline(v=dim(tempTable)[1]*0.1, lwd=2, col="darkgray")
  dev.off()
  
}

#
# Proccedes with Bullard upper quantile normalization using EBSeq function
# See 
#   Bullard, J. H., Purdom, E., Hansen, K. D., & Dudoit, S. (2010).
#   Evaluation of statistical methods for normalization
#   and differential expression in mRNA-Seq experiments.
#   BMC Bioinformatics, 11, 94. http://doi.org/10.1186/1471-2105-11-94
#


DadosLFQ_log2_UQ <- DadosLFQ_log2

# It will generate an vector for normalization
BullardNormalization <- QuantileNorm(DadosLFQ_log2[, AmostrasDisponiveis],.75)

# Get the normalization values and apply to the list
DadosLFQ_log2_UQ[, AmostrasDisponiveis] <- GetNormalizedMat(DadosLFQ_log2[, AmostrasDisponiveis], BullardNormalization)


# Create random color assingment
colorList <- adjustcolor(randomColor(count = length(AmostrasDisponiveis)), alpha.f = 0.8)


dir.create(path = paste0("./", "Results", nameModifier, "/",
                         "01_NormalizeInfo/"), showWarnings = FALSE)


# Plot before normalization
png(file=paste0("./", "Results", nameModifier, "/",
                "01_NormalizeInfo/",
                "01_DistributionBefore", ".png"),
    width=2000, height=1200, res=300)

plot(density((DadosLFQ_log2[,AmostrasDisponiveis][,1])), lwd=2,
     ylim=c(0,0.25), xlim=c(10,35),
     col=colorList[1])

for(j in 2:length(AmostrasDisponiveis)){
  lines(density((DadosLFQ_log2[,AmostrasDisponiveis][,j])), lwd=2, col=colorList[j])
}
dev.off()
#



# Plot AFTER normalization
png(file=paste0("./", "Results", nameModifier, "/",
                "01_NormalizeInfo/",
                "02_DistributionAfterUQ", ".png"),
    width=2000, height=1200, res=300)

plot(density((DadosLFQ_log2_UQ[,AmostrasDisponiveis][,1])), lwd=2,
     ylim=c(0,0.25), xlim=c(10,35),
     col=colorList[1])

for(j in 2:length(AmostrasDisponiveis)){
  lines(density((DadosLFQ_log2_UQ[,AmostrasDisponiveis][,j])), lwd=2, col=colorList[j])
}
dev.off()
#





for(j in 1:length(AmostrasDisponiveis)){
  # Plot AFTER normalization single
  png(file=paste0("./", "Results", nameModifier, "/",
                  "01_NormalizeInfo/",
                  "04_DistributionBefore_", j, "_", AmostrasDisponiveis[j], ".png"),
      width=2000, height=1200, res=300)
  
  plot(density((DadosLFQ_log2[,AmostrasDisponiveis][,j])), lwd=2,
       ylim=c(0,0.25), xlim=c(10,35),
       col=colorList[j],
       main=AmostrasDisponiveis[j])
  dev.off()
}

for(j in 1:length(AmostrasDisponiveis)){
  # Plot AFTER normalization single
  png(file=paste0("./", "Results", nameModifier, "/",
                  "01_NormalizeInfo/",
                  "05_DistributionAfterUQ_", j, "_", AmostrasDisponiveis[j], ".png"),
      width=2000, height=1200, res=300)
  
  plot(density((DadosLFQ_log2_UQ[,AmostrasDisponiveis][,j])), lwd=2,
       ylim=c(0,0.25), xlim=c(10,35),
       col=colorList[j],
       main=AmostrasDisponiveis[j])
  dev.off()
}



#
DadosLFQ_log2 -> DadosLFQ_log2_scaled
as.data.frame(scale(DadosLFQ_log2[,AmostrasDisponiveis])) -> DadosLFQ_log2_scaled[,AmostrasDisponiveis]




# Plot AFTER normalization
png(file=paste0("./", "Results", nameModifier, "/",
                "01_NormalizeInfo/",
                "03_DistributionAfterScales", ".png"),
    width=2000, height=1200, res=300)

plot(density((DadosLFQ_log2_scaled[,AmostrasDisponiveis][,1])), lwd=2,
     ylim=c(0,0.7),
     xlim=c(-3, 5),
     col=colorList[1])

for(j in 2:length(AmostrasDisponiveis)){
  lines(density((DadosLFQ_log2_scaled[,AmostrasDisponiveis][,j])), lwd=2, col=colorList[j])
}
dev.off()
#


for(j in 1:length(AmostrasDisponiveis)){
  # Plot AFTER normalization single
  png(file=paste0("./", "Results", nameModifier, "/",
                  "01_NormalizeInfo/",
                  "06_DistributionAfterScales_", j, "_", AmostrasDisponiveis[j], ".png"),
      width=2000, height=1200, res=300)
  
  plot(density((DadosLFQ_log2_scaled[,AmostrasDisponiveis][,j])), lwd=2,
       #ylim=c(0,0.25),
       #xlim=c(10,35),
       col=colorList[j],
       main=AmostrasDisponiveis[j])
  dev.off()
}



# Save an object to a file
save(ComparisonsTable,
     GroupTable,
     GroupTable_factored,
     GroupTable_summary,
     AmostrasDisponiveis,
     DadosLFQ_log2_UQ,
     nameModifier,
     file = paste0("./", "Results", nameModifier, "/",
                   "01_dataReady.RData"))

write.csv2(DadosLFQ_log2,
           file=paste0("./", "Results", nameModifier, "/",
                       "01_NormalizeInfo/",
                       "99_DadosLFQ_log2", ".csv"))

write.csv2(DadosLFQ_log2_UQ,
           file=paste0("./", "Results", nameModifier, "/",
                       "01_NormalizeInfo/",
                       "99_DadosLFQ_log2_UQ", ".csv"))



##### LIMMA & DEP analysis ####
# Clean environment
rm(list=ls())

library(limma)
library(RColorBrewer)
library(ggpubr)
library(pheatmap)

library(biomaRt)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(graphite)

library(msigdbr)

### Collect only once
m_df1 <- msigdbr(species = "Homo sapiens", category = "H")
m_df2 <- msigdbr(species = "Homo sapiens", category = "C2")
m_df3 <- msigdbr(species = "Homo sapiens", category = "C3")
m_df4 <- msigdbr(species = "Homo sapiens", category = "C8")
m_df <- rbind(m_df1, m_df2)
m_df <- rbind(m_df, m_df3)
m_df <- rbind(m_df, m_df4)


#
nameModifier <- "_01_TodosOsGrupos"

#
load(file = paste0("./", "Results", nameModifier, "/",
                   "01_dataReady.RData"))


#
nameModifier <- "_01_TodosOsGrupos"

##### limma ####
dir.create(path = paste0("./", "Results", nameModifier, "/",
                         "02_limma/"), showWarnings = FALSE)

# Collect only 48h timepoints
targetsNow <- GroupTable_factored[,c("LampTime"), drop=FALSE]

rownames(targetsNow) <- GroupTable_factored$Sample
colnames(targetsNow)[1] <- "Target"

#
#targetsNow <- targetsNow[grep("_48h", targetsNow$Target), ,drop=FALSE]



#
f <- factor(targetsNow$Target, levels=unique(targetsNow$Target))
design <- model.matrix(~0+f)

colnames(design) <- unique(targetsNow$Target)


# constrast loop
myConstrasts <- vector()
for(a1 in 1:length(unique(GroupTable$LampTime))){
  for(b1 in 1:length(unique(GroupTable$LampTime))){
    #
    if(b1>a1){
      
      #if(grepl("_48h", unique(GroupTable$LampTime)[a1]) &
      #   grepl("_48h", unique(GroupTable$LampTime)[b1])
      #){
        myConstrasts <- c(myConstrasts, paste("", unique(GroupTable$LampTime)[a1], "-", "", unique(GroupTable$LampTime)[b1], sep=""))
      #}
    }
  }
}


#
contrast.matrix <- makeContrasts(contrasts=myConstrasts,
                                 levels=unique(as.character(targetsNow$Target)))



eset <- as.matrix(DadosLFQ_log2_UQ[, AmostrasDisponiveis  ])
#eset <- as.matrix(DadosLFQ_log2_UQ[, AmostrasDisponiveis[grep("_48h", AmostrasDisponiveis)]  ])
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

##### limma results save #####
#
#

dir.create(path = paste0("./", "Results", nameModifier, "/",
                         "02_limma/",
                         "01_Tables/"), showWarnings = FALSE)
list() -> Results_list

SummaryComparison <- data.frame(Comparison=myConstrasts,
                                P.Value=0,
                                adj.P.Val=0)

sink(paste0("./", "Results", nameModifier, "/",
            "02_limma/",
            "01_Tables/", "PValues_summary", ".txt") )

for(ConstrastNow in myConstrasts){
  
  ConstrastNow_i <- match(ConstrastNow, myConstrasts)
  
  #
  ResultTable <- topTable(fit2,
                          coef=ConstrastNow_i,
                          adjust="BH",
                          number=dim(eset)[1])
  
  #
  #
  ResultTable <- cbind(DadosLFQ_log2_UQ[rownames(ResultTable), c("GeneID_final",
                                                                 "GeneSymbol_final",
                                                                 "Protein.name",
                                                                 "Fasta.headers")],
                       ResultTable)
  
  SummaryComparison[ConstrastNow_i, "P.Value"] <- as.numeric(table(ResultTable[,"P.Value"] < 0.05)["TRUE"])
  SummaryComparison[ConstrastNow_i, "adj.P.Val"] <- as.numeric(table(ResultTable[,"adj.P.Val"] < 0.05)["TRUE"])
  

  print(
    paste0(
      "Para a condição ",
      ConstrastNow,
      " temos ",
      table(ResultTable[,"P.Value"] < 0.05)["TRUE"],
      " proteinas com  p-valor e ",
      table(ResultTable[,"adj.P.Val"] < 0.05)["TRUE"],
      " com p-valor ajustado"
    )
  )
  
  write.csv2(ResultTable,
             file=paste0("./", "Results", nameModifier, "/",
                         "02_limma/",
                         "01_Tables/",
                         ConstrastNow, ".csv"))
  
  Results_list[[ (length(Results_list)+1) ]] <- ResultTable
}

sink()

#
names(Results_list) <- myConstrasts

#
write.csv2(SummaryComparison,
           paste0("./", "Results", nameModifier, "/",
                  "02_limma/",
                  "Altered_genes_pvalor", ".csv"))

SummaryComparison$Comparison -> row.names(SummaryComparison)


#### nonAdjPvalue
as.vector(SummaryComparison[,"P.Value"]) -> myBarPlot
names(myBarPlot) <- row.names(SummaryComparison)
##### Save Altered Genes
png(file=paste0("./", "Results", nameModifier, "/",
                "02_limma/",
                "Altered_genes_pvalor", ".png"),
    width=5000, height=(660+(110*length(myBarPlot))), res = 400)
par(mar=c(5,40,2,2),lwd = 2.5)
barplot(myBarPlot,
        main="Significant Genes",
        ylab="",
        xlab="Genes number",
        las=2,
        horiz=TRUE,
        cex.lab=1.7,
        cex.axis=2,
        cex.main=1.5,
        axes=FALSE,
        col=brewer.pal(8,"Set1")[1], 
        border="white",
        cex.names=1.5,
        space=0.001
)
#abline(v=1:floor(max(myBarPlot, na.rm=TRUE)  ),
#       col="white", lwd=4)
axis(side = 1, lwd = 1.5, cex.axis=1.4)
dev.off()



#### AdjPvalue
as.vector(SummaryComparison[,"adj.P.Val"]) -> myBarPlot
names(myBarPlot) <- row.names(SummaryComparison)
##### Save Altered Genes
png(file=paste0("./", "Results", nameModifier, "/",
                "02_limma/",
                "Altered_genes_adjpvalor", ".png"),
    width=5000, height=(660+(110*length(myBarPlot))), res = 400)
par(mar=c(5,40,2,2),lwd = 2.5)
barplot(myBarPlot,
        main="Significant Genes",
        ylab="",
        xlab="Genes number",
        las=2,
        horiz=TRUE,
        cex.lab=1.7,
        cex.axis=2,
        cex.main=1.5,
        axes=FALSE,
        col=brewer.pal(8,"Set1")[1], 
        border="white",
        cex.names=1.5,
        space=0.001
)
#abline(v=1:floor(max(myBarPlot, na.rm=TRUE)  ),
#       col="white", lwd=4)
axis(side = 1, lwd = 1.5, cex.axis=1.4)
dev.off()

##
# 


###### Comparison analysis ######
GroupTable_factored$Sample -> row.names(GroupTable_factored)

dir.create(path = paste0("./", "Results", nameModifier, "/",
                         "02_limma/",
                         "02_Comp/"), showWarnings = FALSE)

###
DEGenes <- vector()


### names(Results_list)[1:length(names(Results_list))]
for(ComparisonNow in c("L0nm_48h-L254nm_48h",
                       "L222nm_48h-L254nm_48h",
                       "L0nm_48h-L222nm_48h")){
  
  
  dir.create(path = paste0("./", "Results", nameModifier, "/",
                           "02_limma/",
                           "02_Comp/",
                           ComparisonNow, "/"), showWarnings = FALSE)
  
  
  dir.create(path = paste0("./", "Results", nameModifier, "/",
                           "02_limma/",
                           "02_Comp/",
                           ComparisonNow, "/",
                           "0_Proteins", "/"), showWarnings = FALSE)
  dir.create(path = paste0("./", "Results", nameModifier, "/",
                           "02_limma/",
                           "02_Comp/",
                           ComparisonNow, "/",
                           "0_GEA", "/"), showWarnings = FALSE)
  
  
  #
  ComparisonNow_i <- match(ComparisonNow, names(Results_list))
  
  #
  tempTable <- Results_list[[ComparisonNow]]
  
  #
  tempTable[,"DE"] <- FALSE
  tempTable[which(tempTable[, "P.Value"] < 0.05),"DE"] <- TRUE
  
  #
  DEGenes <- c(DEGenes, rownames(tempTable[tempTable$DE,]))
  
  #
  tempTable$DEP <- "no"
  #tempTable$adj.P.Val
  tempTable$DEP[which(tempTable$P.Value < 0.05 &
                        tempTable$logFC > 0)] <- "up"
  tempTable$DEP[which(tempTable$P.Value < 0.05 &
                        tempTable$logFC < 0)] <- "down"
  
  
  
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("down", "up", "no")
  
  #
  png(paste0("./", "Results", nameModifier, "/",
             "02_limma/",
             "02_Comp/",
             ComparisonNow, "/Volcano_pvalue.png"),
         width=1400,
         height=1000,
         res=300)
  
  print(ggplot(data=tempTable, aes(x=logFC,
                                   y=-log10(P.Value),
                                   col=DEP
                                   #,
                                   #label=delabel
  )) + 
    geom_point() + 
    theme_minimal() +
    scale_colour_manual(values = mycolors) +
    geom_vline(xintercept=c(0), col="black") +
    geom_hline(yintercept=-log10(0.05), col="red"))
  
  dev.off()
  
  
  
  
  #
  tempTable$DEP <- "no"
  #tempTable$adj.P.Val
  tempTable$DEP[which(tempTable$adj.P.Val < 0.05 &
                        tempTable$logFC > 0)] <- "up"
  tempTable$DEP[which(tempTable$adj.P.Val < 0.05 &
                        tempTable$logFC < 0)] <- "down"
  
  
  
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("down", "up", "no")
  
  #
  png(paste0("./", "Results", nameModifier, "/",
             "02_limma/",
             "02_Comp/",
             ComparisonNow, "/Volcano_adjpvalue.png"),
      width=1400,
      height=1000,
      res=300)
  
  
  print(
    ggplot(data=tempTable, aes(x=logFC,
                               y=-log10(adj.P.Val),
                               col=DEP
                               #,
                               #label=delabel
    )) + 
      geom_point() + 
      theme_minimal() +
      scale_colour_manual(values = mycolors) +
      geom_vline(xintercept=c(0), col="black") +
      geom_hline(yintercept=-log10(0.05), col="red")
  )
  
  
  dev.off()
  
  #
  print(paste0("Plotting proteins from ", ComparisonNow))
  
  
  print(paste0("Proteins to be plotted: ", length(rownames(tempTable[tempTable$DE,]))) )
  ## Plot each gene comparison
  for(protNow in rownames(tempTable[tempTable$DE,]) ){
    
    barplotTemp <- data.frame(Group=GroupTable_factored$LampTime,
                              value=as.numeric(DadosLFQ_log2_UQ[protNow, as.character(GroupTable_factored$Sample)]),
                              Compared="No")
    "Yes" -> barplotTemp$Compared[which(as.character(GroupTable_factored$LampTime) %in% strsplit(ComparisonNow, "-")[[1]])]
    
    # Plot ggpubr
    # par(mar=c(5.1,20.1,4.1,2.1))
    if("A"=="B"){
      p <- ggbarplot(barplotTemp,
                     fill = "Compared",
                     x = "Group",
                     y = "value",
                     title = protNow,
                     ylab = "log2(LFQ-UQ)",
                     add = c("mean_se"),
                     size=0.4,
                     add.params = list(size = 0.4,
                                       alpha = 1,
                                       color="black")) +
        rotate_x_text(90) +
        geom_jitter(color = "black",
                    width = 0.2,
                    size = 0.9,
                    alpha = 0.4,
                    shape=16)
      png(paste0("./", "Results", nameModifier, "/",
                 "02_limma/",
                 "02_Comp/",
                 ComparisonNow, "/",
                 "0_Proteins", "/",
                 protNow, ".png"),
          width=(260+(110*length(unique(as.character(barplotTemp$Group))))),
          height=1200, res=400)
      print(ggpar(p,
                  palette = "npg"
      ))
      dev.off()
    }

    write.csv2(barplotTemp,
               file=paste0("./", "Results", nameModifier, "/",
                           "02_limma/",
                           "02_Comp/",
                           ComparisonNow, "/",
                           "0_Proteins", "/",
                           protNow, ".csv"))
    
  }
  
  #
  
  print(paste0("Plotting heatmap from ", ComparisonNow))
  ### Plot heatmap
  HeatmapTable <- DadosLFQ_log2_UQ[rownames(tempTable[tempTable$DE,]),
                                   as.character(GroupTable_factored$Sample)]
  
  
  # Overall comparisons
  my_palette <- c(colorRampPalette(c("green","green","green","black"))(20),
                  colorRampPalette(c("black","red","red","red"))(20)
  )
  
  #
  mycolors1 <- c("#e7298a", "#66a61e", "#e6ab02")
  #
  mycolors2 <- c("#1b9e77", "#d95f02", "#7570b3")
  
  
  #
  names(mycolors1) <- c( unique(GroupTable_factored$Lamp) )
  #
  names(mycolors2) <- c( unique(GroupTable_factored$Time) )
  
  
  #
  mycolors <- list(Lamp = mycolors1,
                   Time = mycolors2)
  
  
  #
  heatmap2 <- pheatmap(mat = as.matrix(HeatmapTable),
                       border_color = "grey60",
                       color = my_palette,
                       na_col = "white",
                       cluster_rows = TRUE,
                       cluster_cols = TRUE,
                       scale="row", # c("row", "column", "none")
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "correlation",
                       clustering_method = "complete",
                       annotation_col = GroupTable_factored[,c("Lamp", "Time")],
                       annotation_colors = mycolors
  )
  
  ###
  png(file=paste0("./", "Results", nameModifier, "/",
                  "02_limma/",
                  "02_Comp/",
                  "Heatmap_",
                  ComparisonNow, ".png"),
      width = 1500,
      height=(390+ (24*dim(HeatmapTable)[1]) ),
      res=200)
  print(heatmap2)
  dev.off()
  
  
  #
  dev.off()

  
  
  ##
  
  ### Heatmap only conditions now
  #
  heatmap2 <- pheatmap(mat = as.matrix(HeatmapTable[,rownames(GroupTable_factored)[GroupTable_factored$LampTime %in% unique(c(unlist( strsplit(ComparisonNow,  "-")), "L0nm_48h"))]]),
                       border_color = "grey60",
                       color = my_palette,
                       na_col = "white",
                       cluster_rows = TRUE,
                       cluster_cols = TRUE,
                       scale="row", # c("row", "column", "none")
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "correlation",
                       clustering_method = "complete",
                       annotation_col = GroupTable_factored[,c("Lamp", "Time")],
                       annotation_colors = mycolors
  )

  ###
  png(file=paste0("./", "Results", nameModifier, "/",
                  "02_limma/",
                  "02_Comp/",
                  "SimpleHeatmap_",
                  ComparisonNow, ".png"),
      width = 1500,
      height=(800+ (8*dim(HeatmapTable)[1]) ),
      res=200)
  print(heatmap2)
  dev.off()
  
  
  #
  dev.off()
  
  
  
  
  
  ### Plot heatmaponly padj
  HeatmapTable <- DadosLFQ_log2_UQ[rownames(tempTable[ which(tempTable$adj.P.Val < 0.05) ,]),
                                   as.character(GroupTable_factored$Sample) ]
  
  
  if(dim(HeatmapTable)[1]>2){
    ### Heatmap only conditions now
    #
    heatmap2 <- pheatmap(mat = as.matrix(HeatmapTable[,rownames(GroupTable_factored)[GroupTable_factored$LampTime %in% unique(c(unlist( strsplit(ComparisonNow,  "-")), "L0nm_48h"))]]),
                         border_color = "grey60",
                         color = my_palette,
                         na_col = "white",
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         scale="row", # c("row", "column", "none")
                         clustering_distance_rows = "correlation",
                         clustering_distance_cols = "correlation",
                         clustering_method = "complete",
                         annotation_col = GroupTable_factored[,c("Lamp", "Time")],
                         annotation_colors = mycolors
    )
    
    ###
    png(file=paste0("./", "Results", nameModifier, "/",
                    "02_limma/",
                    "02_Comp/",
                    "SimpleHeatmap_padjusted_",
                    ComparisonNow, ".png"),
        width = 1500,
        height=(800+ (8*dim(HeatmapTable)[1]) ),
        res=200)
    print(heatmap2)
    dev.off()
    
    
    #
    dev.off()
  }
  
}














###### Overall heatmap #####

#
unique(DEGenes) -> DEGenes


### Plot heatmap
HeatmapTable <- DadosLFQ_log2_UQ[DEGenes,
                                 as.character(GroupTable_factored$Sample)]


# Overall comparisons
my_palette <- c(colorRampPalette(c("green","green","green","black"))(20),
                colorRampPalette(c("black","red","red","red"))(20)
)

#
mycolors1 <- c("#e7298a", "#66a61e", "#e6ab02")
#
mycolors2 <- c("#1b9e77", "#d95f02", "#7570b3")


#
names(mycolors1) <- c( unique(GroupTable_factored$Lamp) )
#
names(mycolors2) <- c( unique(GroupTable_factored$Time) )

#
mycolors <- list(Lamp = mycolors1,
                 Time = mycolors2)


#
heatmap2 <- pheatmap(mat = as.matrix(HeatmapTable),
                     border_color = "grey60",
                     color = my_palette,
                     na_col = "white",
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     scale="row", # c("row", "column", "none")
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "correlation",
                     clustering_method = "complete",
                     annotation_col = GroupTable_factored[,c("Lamp", "Time")],
                     annotation_colors = mycolors
)

###
png(file=paste0("./", "Results", nameModifier, "/",
                "02_limma/",
                "02_Comp/",
                "Heatmap_",
                "0-All_DE_Genes", ".png"),
    width = 1500,
    height=(390+ (24*dim(HeatmapTable)[1]) ),
    res=200)
print(heatmap2)
dev.off()


#
dev.off()

#
# Save an object to a file
save(ComparisonsTable,
     GroupTable,
     GroupTable_factored,
     GroupTable_summary,
     AmostrasDisponiveis,
     DadosLFQ_log2_UQ,
     nameModifier,
     file = paste0("./", "Results", nameModifier, "/",
                   "02_dataReady.RData"))




##### Pathway analysis #####

##### Load packages #####
rm(list=ls())



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

##### Collect MSigDB #####

#

### Collect only once
m_df <- rbind( msigdbr(species = "Homo sapiens", category = "H"),
               msigdbr(species = "Homo sapiens", category = "C1") )

m_df <- rbind( m_df,
               msigdbr(species = "Homo sapiens", category = "C2") )

m_df <- rbind( m_df,
               msigdbr(species = "Homo sapiens", category = "C3") )

m_df <- rbind( m_df,
               msigdbr(species = "Homo sapiens", category = "C4") )

m_df <- rbind( m_df,
               msigdbr(species = "Homo sapiens", category = "C5") )

m_df <- rbind( m_df,
               msigdbr(species = "Homo sapiens", category = "C6") )

m_df <- rbind( m_df,
               msigdbr(species = "Homo sapiens", category = "C7") )

m_df <- rbind( m_df,
               msigdbr(species = "Homo sapiens", category = "C8") )

data.frame(m_df) -> m_df


##### Collect Data and Pathways #####
# A
dir.create(path=paste0("./Results_pathways/"))

# A
AbundanceData <- read.csv2(paste0("./Results_01_TodosOsGrupos/01_NormalizeInfo/99_DadosLFQ_log2_UQ.csv"),
                           row.names=1)

#
Groups <- read.csv2(paste0("./Agrupamento_01.csv"),
                    row.names=1)

m_df$Direction <- ""

SelectedPathways <- m_df[, c("gs_name",
                             "Direction",
                             "gs_name",
                             "gs_id",
                             "gs_description")]

colnames(SelectedPathways) <- c("Group.Name",
                                "Direction",
                                "Standard_name",
                                "Systematic_name",
                                "Source")

unique(SelectedPathways) -> SelectedPathways

#
SelectedPathways$Direction <- "Up"

#
SelectedPathways <- SelectedPathways[SelectedPathways$Systematic_name != "",]

#
SelectedPathways$Group.Name <- str_to_title(gsub("_", " ", SelectedPathways$Group.Name))

#
SelectedPathways$Group.Name <- gsub(" Dn$", "", SelectedPathways$Group.Name)
SelectedPathways$Group.Name <- gsub(" Up$", "", SelectedPathways$Group.Name)


# Label accordingly
SelectedPathways$PresentProteins <- 0

#
if("A" == "B"){
  SelectedPathways$PresentProteins[1:1000] <- unname(sapply(
    1:1000, #dim(SelectedPathways)[1]
    function(pathwayNow){length(which(m_df[ which(m_df$gs_id %in% SelectedPathways$Systematic_name[pathwayNow]), "human_entrez_gene"] %in% AbundanceData$GeneID_final))}))
  
  #
  SelectedPathways$PresentProteins[1001:5000] <- unname(sapply(
    1001:5000, #dim(SelectedPathways)[1]
    function(pathwayNow){length(which(m_df[ which(m_df$gs_id %in% SelectedPathways$Systematic_name[pathwayNow]), "human_entrez_gene"] %in% AbundanceData$GeneID_final))}))
  
  #
  SelectedPathways$PresentProteins[5001:10000] <- unname(sapply(
    5001:10000, #dim(SelectedPathways)[1]
    function(pathwayNow){length(which(m_df[ which(m_df$gs_id %in% SelectedPathways$Systematic_name[pathwayNow]), "human_entrez_gene"] %in% AbundanceData$GeneID_final))}))
  
  #
  SelectedPathways$PresentProteins[10001:20000] <- unname(sapply(
    10001:20000, #dim(SelectedPathways)[1]
    function(pathwayNow){length(which(m_df[ which(m_df$gs_id %in% SelectedPathways$Systematic_name[pathwayNow]), "human_entrez_gene"] %in% AbundanceData$GeneID_final))}))
  
  #
  SelectedPathways$PresentProteins[20001:25000] <- unname(sapply(
    20001:25000, #dim(SelectedPathways)[1]
    function(pathwayNow){length(which(m_df[ which(m_df$gs_id %in% SelectedPathways$Systematic_name[pathwayNow]), "human_entrez_gene"] %in% AbundanceData$GeneID_final))}))
  
  #
  SelectedPathways$PresentProteins[25001:dim(SelectedPathways)[1]] <- unname(sapply(
    25001:dim(SelectedPathways)[1], #dim(SelectedPathways)[1]
    function(pathwayNow){length(which(m_df[ which(m_df$gs_id %in% SelectedPathways$Systematic_name[pathwayNow]), "human_entrez_gene"] %in% AbundanceData$GeneID_final))}))
  
  
  #
  write.xlsx(SelectedPathways,
             "./Results_pathways/SelectedPathways.xlsx")
  
}


#
SelectedPathways <- read.xlsx("./Results_pathways/SelectedPathways.xlsx")


# Add keep column
SelectedPathways$Keep <- FALSE

dim(SelectedPathways)

#
k <- 0
for(pathwayNow in unique(SelectedPathways$Group.Name)  ){
  k <- k+1
  print(k)
  
  #
  if( sum(SelectedPathways[which(SelectedPathways$Group.Name == pathwayNow), "PresentProteins"]) > 4 ){
    #
    SelectedPathways[which(SelectedPathways$Group.Name == pathwayNow), "Keep"] <- TRUE
  }
}


#
SelectedPathways <- SelectedPathways[SelectedPathways$Keep,]



#####
for(pathwayNow in 1:dim(SelectedPathways)[1]){
  print(pathwayNow)
  
  # SelectedPathways$PresentProteins[pathwayNow] <- length(which(m_df[ which(m_df$gs_id %in% SelectedPathways$Systematic_name[pathwayNow]), "human_entrez_gene"] %in% AbundanceData$GeneID_final))
  
  if(SelectedPathways$Group.Name[pathwayNow] != str_to_title(gsub("_", " ", SelectedPathways$Standard_name[pathwayNow])) ){
    
    if( grepl("_DN", SelectedPathways$Standard_name[pathwayNow]) ){
      SelectedPathways$Direction[pathwayNow] <- "Down"
    }
    if( grepl("_UP", SelectedPathways$Standard_name[pathwayNow]) ){
      SelectedPathways$Direction[pathwayNow] <- "Up"
    } 
  }
  
}

#
dim(SelectedPathways)

#
write.xlsx(SelectedPathways,
           file="./Results_pathways/FilteredPathways.xlsx")




##### Loop to score calculation #####
rankData <- rankGenes(expreMatrix=AbundanceData[, rownames(Groups) ],
                      tiesMethod = "min",
                      stableGenes = NULL)


#grep("Controle_", colnames(AbundanceData))
#grep("Lampada_", colnames(AbundanceData))

Groups -> ScoreData

Groups -> GroupTable_factored

GroupTable_factored$Lamp <- as.factor(GroupTable_factored$Lamp)
GroupTable_factored$LampTime <- as.factor(GroupTable_factored$LampTime)
GroupTable_factored$Time    <- as.factor(GroupTable_factored$Time)

#
#SelectedPathways <- SelectedPathways[!(SelectedPathways$Group.Name %in% c("Hallmarka Repair",
#                                                                          "Hallmarka Repair",
#                                                                          "Acevedo Methylated In Liver Cancer",
#                                                                          "Boyault Liver Cancer Subclass G5",
#                                                                          "Choi Atl Chronic Vs Acute")),]


#### Scoring Loop
w <- 0
if(exists("FinalGenePerPathTable")){
  remove(FinalGenePerPathTable)
}

for(pathwayNow in unique(SelectedPathways$Group.Name)[1:length(unique(SelectedPathways$Group.Name))]  ){
  
  #
  w <- w+1
  #
  print(paste0(w, " ", pathwayNow))
  
  #####
  # DANG_MYC_TARGETS
  UpGenes <- vector()
  DownGenes <- vector()
  
  k <- 0
  for(directNow in SelectedPathways[which(SelectedPathways$Group.Name == pathwayNow), "Direction"]){
    k <- k+1
    #
    if(directNow == "Up"){
      #
      #print("Up")
      UpGenes <- c(UpGenes, row.names(AbundanceData)[which(AbundanceData$GeneID_final %in% m_df[which(m_df$gs_id %in% SelectedPathways[which(SelectedPathways$Group.Name == pathwayNow)[k], "Systematic_name"]), "human_entrez_gene"])])
    }
    if(directNow == "Down"){
      #print("Down")
      DownGenes <- c(DownGenes, row.names(AbundanceData)[which(AbundanceData$GeneID_final %in% m_df[which(m_df$gs_id %in% SelectedPathways[which(SelectedPathways$Group.Name == pathwayNow)[k], "Systematic_name"]), "human_entrez_gene"])])
    }
  }
  
  
  #
  unique(UpGenes) -> UpGenes
  unique(DownGenes) -> DownGenes
  
  #
  #
  if(exists("TempTable")){
    remove(TempTable)
  } else {
    #
  }
  
  # Gene table  
  TempTable <- m_df[which(m_df$gs_id %in% SelectedPathways[which(SelectedPathways$Group.Name == pathwayNow), "Systematic_name"]),    ]
  TempTable[which(TempTable$human_gene_symbol %in% UpGenes), "Direction"] <- "Up"
  TempTable[which(TempTable$human_gene_symbol %in% DownGenes), "Direction"] <- "Down"
  
  
  
  # 
  if(exists("scoredf")){
    remove(scoredf)
  }
  
  # Calculate score
  if( ( length(UpGenes) > 0 & length(DownGenes) > 0 ) ){
    
    #
    
    #
    scoredf <- simpleScore(rankData,
                           upSet = GeneSet(geneIds = UpGenes,
                                           geneIdType = SymbolIdentifier(),
                                           setName = pathwayNow
                           ),
                           downSet = GeneSet(geneIds = DownGenes,
                                             geneIdType = SymbolIdentifier(),
                                             setName = pathwayNow
                           ) )
  } else {
    if(length(UpGenes) > 0){
      #
      scoredf <- simpleScore(rankData,
                             upSet = GeneSet(geneIds = UpGenes,
                                             geneIdType = SymbolIdentifier(),
                                             setName = pathwayNow
                             ) )
    }
    
    if(length(DownGenes) > 0){
      #
      scoredf <- simpleScore(rankData,
                             upSet = GeneSet(geneIds = DownGenes,
                                             geneIdType = SymbolIdentifier(),
                                             setName = pathwayNow
                             ) )
    }
  }
  
  # Add score to the table
  #SelectedPathways[which(SelectedPathways$Group.Name == pathwayNow), "Standard_name"]
  ScoreData[, gsub(" ", "_", pathwayNow) ] <- scoredf[rownames(ScoreData), "TotalScore"]
  
  #
  if(exists("FinalGenePerPathTable")){
    #
    FinalGenePerPathTable <- rbind(FinalGenePerPathTable, TempTable)
  } else {
    #
    FinalGenePerPathTable <- TempTable
  }
  
  #
  
}


dim(FinalGenePerPathTable)

# Classify genes
FinalGenePerPathTable[which(FinalGenePerPathTable$Direction == ""), "Direction"] <- "NotInProteomics"
#
write.csv2(FinalGenePerPathTable,
           file="./Results_pathways/PathwaysWithGenes.csv")


write.xlsx(ScoreData,
           file="./Results_pathways/ScoreData.xlsx")





# read it back
FinalGenePerPathTable <- read.csv2(file="./Results_pathways/PathwaysWithGenes.csv",
                                   row.names = 1, stringsAsFactors = FALSE)

#
ScoreData <- read.xlsx("./Results_pathways/ScoreData.xlsx")



#ScoreData[1:5, 1:5]

##### t-test loop #####


# Remove columns with only zeroes
ScoreData <- ScoreData[, c(1:3, (which( colSums(ScoreData[,4:dim(ScoreData)[2]]) > 0 )+3)  )  ]



#
PathwayEnrichmentGSEA <- data.frame(row.names=colnames(ScoreData)[ !(colnames(ScoreData) %in% c("Lamp", "Time", "LampTime"))  ])

#
PathwayEnrichmentGSEA[,"FoldChange"] <- ""
PathwayEnrichmentGSEA[,"pvalue"] <- ""


PathwayEnrichmentGSEA[,"FoldChange_222_toCtl"] <- ""
PathwayEnrichmentGSEA[,"pvalue_222_toCtl"] <- ""

PathwayEnrichmentGSEA[,"FoldChange_254_toCtl"] <- ""
PathwayEnrichmentGSEA[,"pvalue_254_toCtl"] <- ""


dir.create( paste0("./Results_pathways/", "Paths/") )

#
for( PathwayNow in colnames(ScoreData)[ !(colnames(ScoreData) %in% c("Lamp", "Time", "LampTime"))  ] ){
  
  #
  TempTable <- ScoreData[,  c("Lamp", "Time", "LampTime", PathwayNow)    ]
  
  #
  colnames(TempTable)[4] <- "Value"
  
  
  # Select data 222nm vs 254nm
  TempTable2 <- TempTable[ c(TempTable$Time == "48h" &
                               TempTable$Lamp != "0nm" ) , ]
  
  TempTable2$LampTime <- factor(TempTable2$LampTime, levels=c("L222nm_48h", "L254nm_48h"))
  
  
  if(exists("ttestnow")){
    remove(ttestnow)
  }
  
  
  #
  tryCatch(
    expr = {
      # Your code...
      # goes here...
      t.test(formula=Value~LampTime,
             data=TempTable2,
             paired=FALSE,
             var.equal=FALSE) -> ttestnow
    },
    error = function(e){ 
      # (Optional)
      # Do this if an error is caught...
      print("A")
    },
    warning = function(w){
      # (Optional)
      # Do this if an warning is caught...
      #print("B")
    },
    finally = {
      # (Optional)
      # Do this at the end before quitting the tryCatch structure...
      #print("C")
    }
  )
  
  
  if(exists("ttestnow") ){
    #
    PathwayEnrichmentGSEA[PathwayNow, "FoldChange"] <- c(ttestnow$estimate[1]/ttestnow$estimate[2])
    PathwayEnrichmentGSEA[PathwayNow, "pvalue"] <- c(ttestnow$p.value)
    
    
    # Plot if
    if(  ((!is.nan(ttestnow$p.value)) & (ttestnow$p.value < 0.05) & exists("ttestnow"))   ){
      
      # Create mean table
      df_mean <- ScoreData[,c("Lamp", "Time", "LampTime", PathwayNow)] %>% 
        group_by(Time, Lamp) %>% 
        summarize(average = mean(get(PathwayNow))) %>%
        ungroup()
      
      if(nchar(PathwayNow)>40){
        #
        nameNow <- paste0(strsplit(PathwayNow, "")[[1]][1:40], collapse="")
      } else {
        nameNow <- PathwayNow
      }
      
      png(file=paste0("./Results_pathways/", "Paths/", nameNow, ".png"),
          width=1200,
          height=800,
          res=300)
      print(
        # Plot
        ggplot(mapping = aes(x = Time, y = get(PathwayNow), fill = Lamp),
               data = ScoreData[,c("Lamp", "Time", "LampTime", PathwayNow)]) + 
          geom_boxplot(alpha=0.6) +
          ylab(paste0(gsub("_", " ", PathwayNow), " score")) +
          theme_classic() + 
          geom_jitter(width=0.3,
                      alpha=0.3) +
          xlab("Lamp") + 
          facet_wrap(~ Lamp,
                     ncol = 3) +
          geom_point(mapping = aes(x = Time, y = average),
                     data = df_mean,
                     alpha=0) +
          geom_line(mapping = aes(x = Time, y = average, group=Lamp,
                                  color= Lamp),
                    data = df_mean
          ) +
          scale_fill_manual( values=c("#969696", "#525252", "#000000") ) +
          scale_color_manual( values=c("#000000", "#000000", "#000000") ) +
          theme(axis.text.x = element_text(angle = 45,
                                           hjust = 1,
                                           color = "black"),
                axis.text.y = element_text(angle = 45,
                                           hjust = 1,
                                           color = "black")
          )
      )
      dev.off()
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Select data 222nm vs control
  TempTable2 <- TempTable[ c(TempTable$Time == "48h" &
                               TempTable$Lamp != "254nm" ) , ]
  
  TempTable2$LampTime <- factor(TempTable2$LampTime, levels=c("L222nm_48h", "L0nm_48h"))
  
  
  if(exists("ttestnow")){
    remove(ttestnow)
  }
  
  
  #
  tryCatch(
    expr = {
      # Your code...
      # goes here...
      t.test(formula=Value~LampTime,
             data=TempTable2,
             paired=FALSE,
             var.equal=FALSE) -> ttestnow
    },
    error = function(e){ 
      # (Optional)
      # Do this if an error is caught...
      print("A")
    },
    warning = function(w){
      # (Optional)
      # Do this if an warning is caught...
      #print("B")
    },
    finally = {
      # (Optional)
      # Do this at the end before quitting the tryCatch structure...
      #print("C")
    }
  )
  
  
  if(exists("ttestnow") ){
    #
    PathwayEnrichmentGSEA[PathwayNow, "FoldChange_222_toCtl"] <- c(ttestnow$estimate[1]/ttestnow$estimate[2])
    PathwayEnrichmentGSEA[PathwayNow, "pvalue_222_toCtl"] <- c(ttestnow$p.value)
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Select data 254nm vs control
  TempTable2 <- TempTable[ c(TempTable$Time == "48h" &
                               TempTable$Lamp != "222nm" ) , ]
  
  TempTable2$LampTime <- factor(TempTable2$LampTime, levels=c("L254nm_48h", "L0nm_48h"))
  
  
  if(exists("ttestnow")){
    remove(ttestnow)
  }
  
  
  #
  tryCatch(
    expr = {
      # Your code...
      # goes here...
      t.test(formula=Value~LampTime,
             data=TempTable2,
             paired=FALSE,
             var.equal=FALSE) -> ttestnow
    },
    error = function(e){ 
      # (Optional)
      # Do this if an error is caught...
      print("A")
    },
    warning = function(w){
      # (Optional)
      # Do this if an warning is caught...
      #print("B")
    },
    finally = {
      # (Optional)
      # Do this at the end before quitting the tryCatch structure...
      #print("C")
    }
  )
  
  
  if(exists("ttestnow") ){
    #
    PathwayEnrichmentGSEA[PathwayNow, "FoldChange_254_toCtl"] <- c(ttestnow$estimate[1]/ttestnow$estimate[2])
    PathwayEnrichmentGSEA[PathwayNow, "pvalue_254_toCtl"] <- c(ttestnow$p.value)
    
  }
  
  
  
  
}





### 


#
PathwayEnrichmentGSEA$FDR <- p.adjust(p=PathwayEnrichmentGSEA$pvalue, method = "fdr")

#
PathwayEnrichmentGSEA$pvalue_signif <- "no"
PathwayEnrichmentGSEA$pvalue_signif[PathwayEnrichmentGSEA$pvalue < 0.05] <- "yes"
PathwayEnrichmentGSEA$fdr_signif <- "no"
PathwayEnrichmentGSEA$fdr_signif[PathwayEnrichmentGSEA$FDR < 0.05] <- "yes"

table(PathwayEnrichmentGSEA$pvalue_signif)

table(PathwayEnrichmentGSEA$fdr_signif)

#
PathwayEnrichmentGSEA$FoldChange <- as.numeric(PathwayEnrichmentGSEA$FoldChange)
PathwayEnrichmentGSEA$pvalue <- as.numeric(PathwayEnrichmentGSEA$pvalue)

PathwayEnrichmentGSEA$Use <- PathwayEnrichmentGSEA$pvalue_signif


#
PathwayEnrichmentGSEA$FoldChange_222_toCtl <- as.numeric(PathwayEnrichmentGSEA$FoldChange_222_toCtl)
PathwayEnrichmentGSEA$pvalue_222_toCtl <- as.numeric(PathwayEnrichmentGSEA$pvalue_222_toCtl)

PathwayEnrichmentGSEA$FoldChange_254_toCtl <- as.numeric(PathwayEnrichmentGSEA$FoldChange_254_toCtl)
PathwayEnrichmentGSEA$pvalue_254_toCtl <- as.numeric(PathwayEnrichmentGSEA$pvalue_254_toCtl)

#
write.xlsx(PathwayEnrichmentGSEA,
           paste0("./Results_pathways/", "PathwayEnrichmentGSEA.xlsx"),
           rowNames=TRUE)




# Make a copy and read it back
# Add "yes" to the use column.
PathwayEnrichmentGSEA <- read.xlsx(paste0("./Results_pathways/", "PathwayEnrichmentGSEA_Manual.xlsx"),
                                   rowNames=TRUE)


##


#
FinalGenePerPathTable$Group.Name <- str_to_title(gsub("_", " ", FinalGenePerPathTable$gs_name))
FinalGenePerPathTable$Group.Name <- gsub(" Dn$", "", FinalGenePerPathTable$Group.Name)
FinalGenePerPathTable$Group.Name <- gsub(" Up$", "", FinalGenePerPathTable$Group.Name)
FinalGenePerPathTable$Group.Name <- gsub(" ", "_", FinalGenePerPathTable$Group.Name)

write.xlsx(x=FinalGenePerPathTable[which(FinalGenePerPathTable$Group.Name %in% rownames(PathwayEnrichmentGSEA)[which(PathwayEnrichmentGSEA$Use == "yes")]),],
           file="Pathways_with_genes_significant_selected.xlsx")

write.xlsx(x=FinalGenePerPathTable[which(FinalGenePerPathTable$Group.Name %in% rownames(PathwayEnrichmentGSEA)[which(PathwayEnrichmentGSEA$pvalue_signif == "yes")]),],
           file="Pathways_with_genes_significant.xlsx")

###


ScoreData[
  with(ScoreData, order(Lamp, Time, LampTime)),
] -> ScoreData

ScoreData[,c("Lamp", "Time", "LampTime")] <- NULL

ScoreData <- ScoreData[,which(colSums(ScoreData) != 0) ]


# Add rowname
# Add value to the score
ScoreData <- cbind(
  GroupTable_factored[
    with(GroupTable_factored, order(Lamp, Time, LampTime)),
  ], ScoreData
)


##### Heatmap plotting ######
#
mycolors1 <- c("#e7298a", "#66a61e", "#e6ab02")
#
mycolors2 <- c("#1b9e77", "#d95f02", "#7570b3")

#
names(mycolors1) <- c( unique(GroupTable_factored$Lamp) )
#
names(mycolors2) <- c( unique(GroupTable_factored$Time) )



# Overall comparisons
my_palette <- c(colorRampPalette(c("green","green","green","black"))(20),
                colorRampPalette(c("black","red","red","red"))(20)
)

#
mycolors <- list(Lamp = mycolors1,
                 Time = mycolors2)





####### All ##########

HeatmapNow <- ScoreData[,row.names(PathwayEnrichmentGSEA)[PathwayEnrichmentGSEA$pvalue_signif == "yes"]]


HeatmapNow[,c("Lamp", "Time", "LampTime")] <- NULL

### Todos os timepoints
heatmap2 <- pheatmap(mat = as.matrix(t(HeatmapNow)),
                     border_color = "grey60",
                     color = my_palette,
                     na_col = "white",
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     scale="row", # c("row", "column", "none")
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "correlation",
                     clustering_method = "complete",
                     annotation_col = GroupTable_factored[,c("Lamp", "Time")],
                     annotation_colors = mycolors
)



###
png(file=paste0("./", "Results_pathways", "/",
                "0-AllTimepoints", ".png"),
    width = 4500,
    height=(390+ (24*dim(ScoreData)[1]) ),
    res=200)
print(heatmap2)
dev.off()


#
dev.off()

### Only 0h
heatmap2 <- pheatmap(mat = as.matrix(t(HeatmapNow[grep("_0h", rownames(HeatmapNow)), ])),
                     border_color = "grey60",
                     color = my_palette,
                     na_col = "white",
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     scale="row", # c("row", "column", "none")
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "correlation",
                     clustering_method = "complete",
                     annotation_col = GroupTable_factored[,c("Lamp", "Time")],
                     annotation_colors = mycolors
)

###
png(file=paste0("./", "Results_pathways", "/",
                "1-Timepoint0h", ".png"),
    width = 2500,
    height=(390+ (38*dim(ScoreData)[1]) ),
    res=200)
print(heatmap2)
dev.off()


#
dev.off()







### Only 24h
heatmap2 <- pheatmap(mat = as.matrix(t(HeatmapNow[grep("_24h", rownames(HeatmapNow)), ])),
                     border_color = "grey60",
                     color = my_palette,
                     na_col = "white",
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     scale="row", # c("row", "column", "none")
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "correlation",
                     clustering_method = "complete",
                     annotation_col = GroupTable_factored[,c("Lamp", "Time")],
                     annotation_colors = mycolors
)

###
png(file=paste0("./", "Results_pathways", "/",
                "1-Timepoint24h", ".png"),
    width = 2500,
    height=(390+ (38*dim(ScoreData)[1]) ),
    res=200)
print(heatmap2)
dev.off()


#
dev.off()








### Only 48h
heatmap2 <- pheatmap(mat = as.matrix(t(HeatmapNow[grep("_48h", rownames(HeatmapNow)), ])),
                     border_color = "grey60",
                     color = my_palette,
                     na_col = "white",
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     scale="row", # c("row", "column", "none")
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "correlation",
                     clustering_method = "complete",
                     annotation_col = GroupTable_factored[,c("Lamp", "Time")],
                     annotation_colors = mycolors
)

###
png(file=paste0("./", "Results_pathways", "/",
                "1-Timepoint48h", ".png"),
    width = 2500,
    height=(390+ (38*dim(ScoreData)[1]) ),
    res=200)
print(heatmap2)
dev.off()


#
dev.off()




##### Only selected #####

HeatmapNow <- ScoreData[,row.names(PathwayEnrichmentGSEA)[PathwayEnrichmentGSEA$Use != "no"]]


### Todos os timepoints
heatmap2 <- pheatmap(mat = as.matrix(t(HeatmapNow)),
                     border_color = "grey60",
                     color = my_palette,
                     na_col = "white",
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     scale="row", # c("row", "column", "none")
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "correlation",
                     clustering_method = "complete",
                     annotation_col = GroupTable_factored[,c("Lamp", "Time")],
                     annotation_colors = mycolors
)



###
png(file=paste0("./", "Results_pathways", "/",
                "2-AllTimepoints", ".png"),
    width = 4500,
    height=(390+ (24*dim(ScoreData)[1]) ),
    res=200)
print(heatmap2)
dev.off()


#
dev.off()












### Apenas 0h
heatmap2 <- pheatmap(mat = as.matrix(t(HeatmapNow[grep("_0h", rownames(HeatmapNow)), ])),
                     border_color = "grey60",
                     color = my_palette,
                     na_col = "white",
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     scale="row", # c("row", "column", "none")
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "correlation",
                     clustering_method = "complete",
                     annotation_col = GroupTable_factored[,c("Lamp", "Time")],
                     annotation_colors = mycolors
)

###
png(file=paste0("./", "Results_pathways", "/",
                "2-Timepoint0h", ".png"),
    width = 2500,
    height=(390+ (38*dim(ScoreData)[1]) ),
    res=200)
print(heatmap2)
dev.off()


#
dev.off()







### Apenas 24h
heatmap2 <- pheatmap(mat = as.matrix(t(HeatmapNow[grep("_24h", rownames(HeatmapNow)), ])),
                     border_color = "grey60",
                     color = my_palette,
                     na_col = "white",
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     scale="row", # c("row", "column", "none")
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "correlation",
                     clustering_method = "complete",
                     annotation_col = GroupTable_factored[,c("Lamp", "Time")],
                     annotation_colors = mycolors
)

###
png(file=paste0("./", "Results_pathways", "/",
                "2-Timepoint24h", ".png"),
    width = 2500,
    height=(390+ (38*dim(ScoreData)[1]) ),
    res=200)
print(heatmap2)
dev.off()


#
dev.off()





# Overall comparisons
my_palette <- c(colorRampPalette(c("green","green","green","black"))(10),
                colorRampPalette(c("black", "red", "red","red"))(10)
)


### Apenas 48h
heatmap2 <- pheatmap(mat = as.matrix(t(HeatmapNow[grep("_48h", rownames(HeatmapNow)), ])),
                     border_color = "grey60",
                     color = my_palette,
                     na_col = "white",
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     scale="row", # c("row", "column", "none")
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method = "complete",
                     annotation_col = GroupTable_factored[,c("Lamp", "Time")],
                     annotation_colors = mycolors
)

###
png(file=paste0("./", "Results_pathways", "/",
                "2-Timepoint48h", ".png"),
    width = 2000,
    height=(390+ (58*dim(ScoreData)[1]) ),
    res=200)
print(heatmap2)
dev.off()


#
dev.off()



##### Line plotting #####

# Add value to the score
# ScoreData <- cbind(GroupTable_factored[,c("Lamp", "Time", "LampTime")], ScoreData)



#
Tukey_test_conditions <- data.frame(row.names = row.names(PathwayEnrichmentGSEA)[PathwayEnrichmentGSEA$pvalue_signif != "no"] )

#
Tukey_test_conditions[,c("L254nm_48h-L222nm_48h",
                         "L222nm_48h-L0nm_48h",
                         "L254nm_48h-L0nm_48h")] <- 1



#
dir.create(paste0("./Results_pathways/Paths2/"))

# 
for(pathwayNow in row.names(PathwayEnrichmentGSEA)[which(PathwayEnrichmentGSEA$pvalue_signif != "no" & PathwayEnrichmentGSEA$Use == "yes")] ){
  
  # Create mean table
  df_mean <- ScoreData[,c("Lamp", "Time", "LampTime", pathwayNow)] %>% 
    group_by(Time, Lamp) %>% 
    summarize(average = mean(get(pathwayNow))) %>%
    ungroup()
  
  
  if(nchar(pathwayNow)>40){
    #
    nameNow <- paste0(strsplit(pathwayNow, "")[[1]][1:40], collapse="")
  } else {
    nameNow <- pathwayNow
  }
  
  
  #
  fm1 <- aov(get(pathwayNow) ~ LampTime,
             data = ScoreData[,c("Lamp",
                                 "Time",
                                 "LampTime",
                                 pathwayNow)])
  
  #
  Tukey_test_conditions[pathwayNow,c("L254nm_48h-L222nm_48h",
                                     "L222nm_48h-L0nm_48h",
                                     "L254nm_48h-L0nm_48h") ] <-  TukeyHSD(fm1)$LampTime[c("L254nm_48h-L222nm_48h",
                                                                                           "L222nm_48h-L0nm_48h",
                                                                                           "L254nm_48h-L0nm_48h"), "p adj"]
  
  
  png(file=paste0("./", "Results_pathways", "/Paths2", "/4-",
                  nameNow, ".png"),
      width=1200,
      height=800,
      res=300)
  
  print(
    
    # Plot
    ggplot(mapping = aes(x = Time, y = get(pathwayNow), fill = Lamp),
           data = ScoreData[,c("Lamp", "Time", "LampTime", pathwayNow)]) + 
      geom_boxplot(alpha=0.6) +
      ylab(paste0(gsub("_", " ", pathwayNow), " score")) +
      theme_classic() + 
      geom_jitter(width=0.3,
                  alpha=0.3) +
      xlab("Lamp") + 
      
      facet_wrap(~ Lamp,
                 ncol = 3) +
      
      geom_point(mapping = aes(x = Time, y = average),
                 data = df_mean,
                 alpha=0) +
      geom_line(mapping = aes(x = Time, y = average, group=Lamp,
                              color= Lamp),
                data = df_mean
      )+
      
      scale_fill_manual( values=c("#969696", "#525252", "#000000") ) +
      scale_color_manual( values=c("#000000", "#000000", "#000000") ) +
      
      theme(axis.text.x = element_text(angle = 45,
                                       hjust = 1,
                                       color = "black"),
            axis.text.y = element_text(angle = 45,
                                       hjust = 1,
                                       color = "black")
      )
    
  )
  
  
  
  dev.off()
  
  
  
}



#
write.xlsx(Tukey_test_conditions,
           rowNames=TRUE,
           file="./Results_pathways/Paths2/TukeyPValues.xlsx")



#### Heatmap selected pathways ####

rm(list=ls())

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
Groups <- read.csv2(paste0("./Agrupamento_01.csv"),
                    row.names=1)

Groups -> GroupTable_factored

GroupTable_factored$Lamp <- as.factor(GroupTable_factored$Lamp)
GroupTable_factored$LampTime <- as.factor(GroupTable_factored$LampTime)
GroupTable_factored$Time    <- as.factor(GroupTable_factored$Time)


# read it back
FinalGenePerPathTable <- read.csv2(file="./Results_pathways/PathwaysWithGenes.csv",
                                   row.names = 1, stringsAsFactors = FALSE)

#
ScoreData <- read.xlsx("./Results_pathways/ScoreData.xlsx")



ScoreData[
  with(ScoreData, order(Lamp, Time, LampTime)),
] -> ScoreData

ScoreData[,c("Lamp", "Time", "LampTime")] <- NULL

ScoreData <- ScoreData[,which(colSums(ScoreData) != 0) ]


# Add rowname
# Add value to the score
ScoreData <- cbind(
  GroupTable_factored[
    with(GroupTable_factored, order(Lamp, Time, LampTime)),
  ], ScoreData
)




#
PathwayEnrichmentGSEA <- read.xlsx(paste0("./Results_pathways/", "PathwayEnrichmentGSEA_Manual.xlsx"),
                                   rowNames=TRUE)



#ScoreData[1:5, 1:5]

##### Heatmap plotting ######
#
mycolors1 <- c("#e7298a", "#66a61e", "#e6ab02")
#
mycolors2 <- c("#1b9e77", "#d95f02", "#7570b3")

#
names(mycolors1) <- c( unique(GroupTable_factored$Lamp) )
#
names(mycolors2) <- c( unique(GroupTable_factored$Time) )
#
mycolors <- list(Lamp = mycolors1,
                 Time = mycolors2)




# Overall comparisons
my_palette <- c(colorRampPalette(c("green","green","green","black"))(10),
                colorRampPalette(c("black", "red", "red","red"))(10)
)




#
row.names(ScoreData) <- row.names(GroupTable_factored)

#
HeatmapNow <- ScoreData[,row.names(PathwayEnrichmentGSEA)[PathwayEnrichmentGSEA$Use != "no"]]

colnames(HeatmapNow) <- gsub("_", " ", colnames(HeatmapNow))

### Apenas 48h
heatmap2 <- pheatmap(mat = as.matrix(t(HeatmapNow[grep("_48h", ScoreData$LampTime   ), ])),
                     border_color = "grey60",
                     color = my_palette,
                     na_col = "white",
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     scale="row", # c("row", "column", "none")
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method = "complete",
                     annotation_col = GroupTable_factored[,c("Lamp", "Time")],
                     annotation_colors = mycolors
)

###
png(file=paste0("./Results_pathways/",
                "2-Timepoint48h", ".png"),
    width = 2200,
    height=(390+ (58*dim(ScoreData)[1]) ),
    res=230)
print(heatmap2)
dev.off()

dev.off()



####



##### Line plotting #####
#
dir.create(path = paste0("./", "Results_pathways/Final_timecourses/"), showWarnings = FALSE)

#
Tukey_test_conditions <- data.frame(row.names = row.names(PathwayEnrichmentGSEA)[PathwayEnrichmentGSEA$Use != "no"] )


#
Tukey_test_conditions[,c("L222nm_48h-L222nm_0h",
                         "L254nm_48h-L254nm_0h",
                         
                         "L254nm_48h-L222nm_48h",
                         "L222nm_48h-L0nm_48h",
                         "L254nm_48h-L0nm_48h",
                         
                         "ttest_L222nm_48h-L222nm_0h",
                         "ttest_L254nm_48h-L254nm_0h",
                         
                         "ttest_L254nm_48h-L222nm_48h",
                         "ttest_L222nm_48h-L0nm_48h",
                         "ttest_L254nm_48h-L0nm_48h"
)] <- 1

#
for(pathwayNow in row.names(PathwayEnrichmentGSEA)[PathwayEnrichmentGSEA$Use != "no"] ){
  
  # Create mean table
  df_mean <- ScoreData[,c("Lamp", "Time", "LampTime", pathwayNow)] %>% 
    group_by(Time, Lamp) %>% 
    summarize(average = mean(get(pathwayNow))) %>%
    ungroup()
  
  if(nchar(pathwayNow)>40){
    #
    nameNow <- paste0(strsplit(pathwayNow, "")[[1]][1:40], collapse="")
  } else {
    nameNow <- pathwayNow
  }
  
  
  #
  fm1 <- aov(get(pathwayNow) ~ LampTime,
             data = ScoreData[,c("Lamp",
                                 "Time",
                                 "LampTime",
                                 pathwayNow)])
  
  #
  Tukey_test_conditions[pathwayNow,c("L222nm_48h-L222nm_0h",
                                     "L254nm_48h-L254nm_0h",
                                     
                                     "L254nm_48h-L222nm_48h",
                                     "L222nm_48h-L0nm_48h",
                                     "L254nm_48h-L0nm_48h") ] <-  TukeyHSD(fm1)$LampTime[c("L222nm_48h-L222nm_0h",
                                                                                           "L254nm_48h-L254nm_0h",
                                                                                           
                                                                                           "L254nm_48h-L222nm_48h",
                                                                                           "L222nm_48h-L0nm_48h",
                                                                                           "L254nm_48h-L0nm_48h"), "p adj"]
  
  #
  Tukey_test_conditions[pathwayNow, "ttest_L222nm_48h-L222nm_0h"] <- t.test(x = ScoreData[which(ScoreData$LampTime == "L222nm_0h"), pathwayNow ],
                                                                            y = ScoreData[which(ScoreData$LampTime == "L222nm_48h"), pathwayNow ],
                                                                            alternative = c("two.sided"),
                                                                            mu = 0,
                                                                            paired = FALSE,
                                                                            var.equal = FALSE,
                                                                            conf.level = 0.95)$p.value
  
  #
  Tukey_test_conditions[pathwayNow, "ttest_L254nm_48h-L254nm_0h"] <- t.test(x = ScoreData[which(ScoreData$LampTime == "L254nm_0h"), pathwayNow ],
                                                                            y = ScoreData[which(ScoreData$LampTime == "L254nm_48h"), pathwayNow ],
                                                                            alternative = c("two.sided"),
                                                                            mu = 0,
                                                                            paired = FALSE,
                                                                            var.equal = FALSE,
                                                                            conf.level = 0.95)$p.value
  
  
  
  
  
  Tukey_test_conditions[pathwayNow, "ttest_L254nm_48h-L222nm_48h"] <- t.test(x = ScoreData[which(ScoreData$LampTime == "L222nm_48h"), pathwayNow ],
                                                                             y = ScoreData[which(ScoreData$LampTime == "L254nm_48h"), pathwayNow ],
                                                                             alternative = c("two.sided"),
                                                                             mu = 0,
                                                                             paired = FALSE,
                                                                             var.equal = FALSE,
                                                                             conf.level = 0.95)$p.value
  Tukey_test_conditions[pathwayNow, "ttest_L222nm_48h-L0nm_48h"] <- t.test(x = ScoreData[which(ScoreData$LampTime == "L0nm_48h"), pathwayNow ],
                                                                           y = ScoreData[which(ScoreData$LampTime == "L222nm_48h"), pathwayNow ],
                                                                           alternative = c("two.sided"),
                                                                           mu = 0,
                                                                           paired = FALSE,
                                                                           var.equal = FALSE,
                                                                           conf.level = 0.95)$p.value
  Tukey_test_conditions[pathwayNow, "ttest_L254nm_48h-L0nm_48h"] <- t.test(x = ScoreData[which(ScoreData$LampTime == "L0nm_48h"), pathwayNow ],
                                                                           y = ScoreData[which(ScoreData$LampTime == "L254nm_48h"), pathwayNow ],
                                                                           alternative = c("two.sided"),
                                                                           mu = 0,
                                                                           paired = FALSE,
                                                                           var.equal = FALSE,
                                                                           conf.level = 0.95)$p.value
  ######
  #
  
  png(file=paste0("./", "./Results_pathways/Final_timecourses/", "", "4-",
                  nameNow, ".png"),
      width=1100,
      height=650,
      res=300)
  
  print(
    
    # Plot
    ggplot(mapping = aes(x = Time, y = get(pathwayNow), fill = Lamp),
           data = ScoreData[,c("Lamp", "Time", "LampTime", pathwayNow)]) + 
      geom_boxplot(alpha=0.6) +
      ylab(paste0(gsub("_", " ", pathwayNow), " score")) +
      theme_classic() + 
      geom_jitter(width=0.3,
                  size=0.8,
                  alpha=0.3) +
      xlab("Lamp") + 
      
      facet_wrap(~ Lamp,
                 ncol = 3) +
      
      geom_point(mapping = aes(x = Time, y = average),
                 data = df_mean,
                 alpha=0) +
      geom_line(mapping = aes(x = Time, y = average, group=Lamp,
                              color= Lamp),
                data = df_mean
      ) +
      
      scale_fill_manual( values=c("#FFFFFF", "#00441B", "#A50F15"    ) ) +
      
      scale_color_manual( values=c("#000000", "#000000", "#000000") ) +
      
      theme(axis.text.x = element_text(angle = 45,
                                       hjust = 1,
                                       color = "black"),
            axis.text.y = element_text(angle = 45,
                                       hjust = 1,
                                       color = "black")
      ) #+
    #geom_signif(
    #   comparisons = list(c("0h", "48h") ),
    #  map_signif_level = TRUE, textsize = 6
    # )
  )
  
  dev.off()
  
}


#
write.xlsx(Tukey_test_conditions,
           rowNames=TRUE,
           file="./Results_pathways/Final_timecourses/TukeyPValues.xlsx")



