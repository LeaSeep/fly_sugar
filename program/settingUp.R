# Load initial data ----
library(SummarizedExperiment)
load("../data/20230818_DESeq_fit_Object (1).Rdata")

MaleSampleData <- as.data.frame(colData(dds)[,-ncol(colData(dds))])
MaleSampleData$Day <- as.factor(gsub("[^0-9.-]", "", MaleSampleData$Group))
write.csv(MaleSampleData, "../data/SampleTableMale.csv")
maleCM_fromDDs <- assay(dds)

load("../data/20230818_DESeq_fit_Object.Rdata")

FeMaleSampleData <- as.data.frame(colData(dds)[,-ncol(colData(dds))])
FeMaleSampleData$Day <- as.factor(gsub("[^0-9.-]", "", FeMaleSampleData$Group))
write.csv(FeMaleSampleData, "../data/SampleTableFemale.csv")
femaleCM_fromDDs <- assay(dds)


# prep SummarizedExperiment objects ----
allSampleData <- rbind(MaleSampleData,FeMaleSampleData)
allSampleData$Day <- as.factor(gsub("[^0-9.-]", "", allSampleData$Group))

write.csv(allSampleData,"SampleData_prepped.csv")
# M44 and M8 added based on Drosophila_RNA_smaples_SBa
allSampleData <- read.csv("../data/SampleData_prepped.csv",row.names = 1)

## load CountMatrix ----

rawCounts <- read.csv("../data/salmon.merged.gene_counts.tsv",sep="\t")
rowAnno <- rawCounts[,1:2]
rownames(rowAnno) <- rowAnno$gene_id

rownames(rawCounts) <- rawCounts$gene_id
rawCounts <- rawCounts[,-c(1:2)]

colnames(rawCounts) <- gsub("_.*","",gsub("^AAC2CNJHV_","",colnames(rawCounts)))
colnames(rawCounts) <- toupper(colnames(rawCounts))

rawCounts <- rawCounts[,!grepl("S",colnames(rawCounts))]
dim(rawCounts)
dim(allSampleData)

rawCounts <- rawCounts[,rownames(allSampleData)]

allSampleData$Diet <- factor(allSampleData$Diet, levels = c("CD","HFD","HSD"))
allSampleData$Gender <- factor(allSampleData$Gender, levels = c("Male","Female"))
allSampleData$Day <- factor(allSampleData$Day, levels = c("0","10","50"))
allSampleData$Condition <- as.factor(allSampleData$Condition)

se_fly_all <- SummarizedExperiment(assays = as.matrix(round(rawCounts)),   ## If salon used should use quant.sf.gz files
                                    colData = allSampleData,
                                    rowData = as.data.frame(rowAnno))


## here we load kallisto data ----
se_fly_all <- dds_kallisto
se_fly_all$Gender <- as.factor(se_fly_all$Gender)
se_fly_all$Condition <- as.factor(se_fly_all$Condition)
se_fly_all$Day <- factor(se_fly_all$Day,levels=c("0","10","50"))
se_fly_all$Diet <- as.factor(se_fly_all$Diet)

## remove Outliers identifed via PCA ----
se_fly_all <- se_fly_all[,!colnames(se_fly_all)%in%c("M8","F12")]

# create sex se
se_fly_male <- se_fly_all[,se_fly_all$Gender=="Male"]
se_fly_female <- se_fly_all[,se_fly_all$Gender=="Female"]

# combine (lets use common ground for a testing)


# Do inital filtering ----

library("DESeq2")

ddsSE <- DESeqDataSet(se_fly_all, design = ~ Condition)
ddsSE$Day <- relevel(ddsSE$Day, ref = "0")
ddsSE$Diet <- relevel(ddsSE$Diet, ref = "CD")

design(ddsSE)
ddsSE_preFilter <- ddsSE
ddsSE <- preprocessing(ddsSE_preFilter,10,
                        protCodingOnly=T,
                        removeConstRows=T,
                        filterPerSample=T)
de_seq_all <- DESeq(ddsSE) 
resultsNames(de_seq_all)
summary(de_seq_all)
res <- results(de_seq_all)

# Do male only ----
ddsSE_male <- DESeqDataSet(se_fly_male, design = ~ Condition)
ddsSE_male$Day <- relevel(ddsSE_male$Day, ref = "0")
ddsSE_male$Diet <- relevel(ddsSE_male$Diet, ref = "CD")

design(ddsSE_male)
ddsSE_preFilter_male <- ddsSE_male
ddsSE_male <- preprocessing(ddsSE_preFilter_male,10,
                       protCodingOnly=T,
                       removeConstRows=T,
                       filterPerSample=T)
de_seq_all_vst_male <- vst(ddsSE_male,blind=T)
library(ggplot2)
PCA_Data<-doPCA(
  de_seq_all_vst_male,
  #xPC="PC2",
  #yPC="PC3",
  colorTheme = colorTheme,
  shapeVar = "Day", # one of colnames in colData(dds)
  colorVar = "Diet" # one of colnames in colData(dds)
)
View(PCA_Data$data)

doPCA(
  de_seq_all_vst_male,
  #xPC="PC2",
  #yPC="PC3",
  colorTheme = colorTheme,
  shapeVar = "Diet", # one of colnames in colData(dds)
  colorVar = "Day" # one of colnames in colData(dds)
)

# Do female only ----
ddsSE_female <- DESeqDataSet(se_fly_female, design = ~ Condition)
ddsSE_female$Day <- relevel(ddsSE_female$Day, ref = "0")
ddsSE_female$Diet <- relevel(ddsSE_female$Diet, ref = "CD")

design(ddsSE_female)
ddsSE_preFilter_male <- ddsSE_female
ddsSE_female <- preprocessing(ddsSE_preFilter_male,10,
                            protCodingOnly=T,
                            removeConstRows=T,
                            filterPerSample=T)
de_seq_all_vst_female <- vst(ddsSE_female,blind=T)
library(ggplot2)
PCA_Data<-doPCA(
  de_seq_all_vst_female,
  xPC="PC1",
  yPC="PC2",
  colorTheme = colorTheme,
  shapeVar = "Day", # one of colnames in colData(dds)
  colorVar = "Diet" # one of colnames in colData(dds)
)
View(PCA_Data$data)
doPCA(
  de_seq_all_vst_female,
  xPC="PC1",
  yPC="PC2",
  colorTheme = colorTheme,
  shapeVar = "Diet", # one of colnames in colData(dds)
  colorVar = "Day" # one of colnames in colData(dds)
)

# PCA all ----
de_seq_all_vst <- vst(de_seq_all,blind=T)
de_seq_all_vst <- vst(ddsSE_split,blind=T)

de_seq_all_vst$Gender = factor(de_seq_all_vst$Gender ,
                            levels = c("Female","Male"),
                            ordered = T)
#reorder coulumns for plotting purposes
idx_2 <- which(grepl("M8|F12",colnames(de_seq_all_vst)))
idx_1 <- which(!grepl("M8|F12",colnames(de_seq_all_vst)))
de_seq_all_vst <- de_seq_all_vst[,c(idx_1,idx_2)]
colData(de_seq_all)
colorTheme <- c("#797979","#F52F92","#4581BE","#e28743","#76b5c5","#dde5b6","#7cb518","#31572c")
names(colorTheme) <- c("CD","HSD","HFD","Male","Female","0","10","50")

library(ggplot2)
doPCA(
  de_seq_all_vst,
 # xPC="PC2",
  #yPC="PC3",
  colorTheme = colorTheme,
  shapeVar = "Diet", # one of colnames in colData(dds)
  colorVar = "Gender" # one of colnames in colData(dds)
)
#View(PCA_Data$data)
doPCA(
  de_seq_all_vst,
  xPC="PC2",
  yPC="PC3",
  colorTheme = colorTheme,
  shapeVar = "Gender", # one of colnames in colData(dds)
  colorVar = "Day" # one of colnames in colData(dds)
)

doPCA(
  de_seq_all_vst,
  xPC="PC2",
  yPC="PC3",
  colorTheme = colorTheme,
  shapeVar = "Day", # one of colnames in colData(dds)
  colorVar = "Diet" # one of colnames in colData(dds)
)


## DE - stats Overlap

# dif male vs female ----
resultsListAll <- list()
namesOnly <- list()
counter=0
for(i in c("0","10","50")){
  for(j in c("CD","HSD","HFD")){
    counter = counter+1
     term = paste0(j,"_Day_",i)
     dif = paste0("sex_",term)
     compA = paste0(term,"_Male_Progeny")
     compB = paste0(term,"_Female_Progeny")
     res_tmp <- results(de_seq_all, contrast = c("Condition",compA,compB))
     resultsListAll[[dif]] <- res_tmp
     
     namesOnly[[paste0(dif,"_UP")]] <- rownames(subset(
       results(de_seq_all, contrast = c("Condition",compA,compB)),
       padj < 0.05 & log2FoldChange > 1))
     namesOnly[[paste0(dif,"_DOWN")]] <- rownames(subset(
       results(de_seq_all, contrast = c("Condition",compA,compB)),
       padj < 0.05 & log2FoldChange < -1))
     print(paste0("Done: ",counter,"/9"))
  }
}

library(UpSetR)
# split it up a bit
names(namesOnly)
UpSetR::upset(fromList(namesOnly), 
              sets = names(namesOnly)[grepl("UP",names(namesOnly))],
              order.by = "freq", 
              #group.by = "sets", 
              cutoff = 10,text.scale = 2)

UpSetR::upset(fromList(namesOnly), 
              sets = names(namesOnly)[grepl("DOWN",names(namesOnly))],
              order.by = "freq", 
              #group.by = "sets", 
              cutoff = 10,text.scale = 2)


# dif within Gender Diet diff ----
library(DESeq2)
resultsListAll_withinGender <- list()
namesOnly_withinGender <- list()

counter = 0
for(i in c("0","10","50")){
  for(j in c("Female","Male")){
    counter = counter+1
    term = paste0("_Day_",i,"_",j,"_Progeny")
    
    difA = paste0("pAdjOnly_Diet_HSD_CD",term)
    difB = paste0("pAdjOnly_Diet_HFD_CD",term)
    difC = paste0("pAdjOnly_Diet_HSD_HFD",term)
    
    compA=paste0("CD",term)
    compB=paste0("HFD",term)
    compC=paste0("HSD",term)
    
    res_tmp <- results(de_seq_all, contrast = c("Condition",compC,compA))
    resultsListAll_withinGender[[difA]] <- res_tmp
    namesOnly_withinGender[[paste0(difA,"_UP")]] <- rownames(subset(
      results(de_seq_all, contrast = c("Condition",compC,compA)),
      padj < 0.05 & log2FoldChange > 0)) #1
    namesOnly_withinGender[[paste0(difA,"_DOWN")]] <- rownames(subset(
      results(de_seq_all, contrast = c("Condition",compC,compA)),
      padj < 0.05 & log2FoldChange < 0))
    
    res_tmp <- results(de_seq_all, contrast = c("Condition",compB,compA))
    resultsListAll_withinGender[[difB]] <- res_tmp
    namesOnly_withinGender[[paste0(difB,"_UP")]] <- rownames(subset(
      results(de_seq_all, contrast = c("Condition",compB,compA)),
      padj < 0.05 & log2FoldChange > 0))
    namesOnly_withinGender[[paste0(difB,"_DOWN")]] <- rownames(subset(
      results(de_seq_all, contrast = c("Condition",compB,compA)),
      padj < 0.05 & log2FoldChange < 0))
    
    res_tmp <- results(de_seq_all, contrast = c("Condition",compB,compC))
    resultsListAll_withinGender[[difC]] <- res_tmp
    namesOnly_withinGender[[paste0(difC,"_UP")]] <- rownames(subset(
      results(de_seq_all, contrast = c("Condition",compB,compC)),
      padj < 0.05 & log2FoldChange > 0))
    namesOnly_withinGender[[paste0(difC,"_DOWN")]] <- rownames(subset(
      results(de_seq_all, contrast = c("Condition",compB,compC)),
      padj < 0.05 & log2FoldChange < 0))
    
    
    
    
    print(paste0("Done: ",counter,"/9"))
  }
}

# split it up a bit
names(namesOnly)
names(namesOnly_withinGender)
UpSetR::upset(fromList(namesOnly), 
              sets = names(namesOnly)[grepl("Diet.*UP",names(namesOnly))],
              order.by = "freq", 
              #group.by = "sets", 
              cutoff = 10,text.scale = 2)

UpSetR::upset(fromList(namesOnly), 
              sets = names(namesOnly)[grepl("Diet.*Female_Progeny_UP",names(namesOnly))],
              order.by = "freq", 
              #group.by = "sets", 
              cutoff = 10,text.scale = 2)

UpSetR::upset(fromList(namesOnly), 
              sets = names(namesOnly)[grepl("Diet.*Male_Progeny_UP",names(namesOnly))],
              order.by = "freq", 
              #group.by = "sets", 
              cutoff = 10,text.scale = 2)

UpSetR::upset(fromList(namesOnly), 
              sets = names(namesOnly)[grepl("Diet.*DOWN",names(namesOnly))],
              order.by = "freq", 
              #group.by = "sets", 
              cutoff = 10,text.scale = 2)

UpSetR::upset(fromList(namesOnly), 
              sets = names(namesOnly)[grepl("Diet.*Female_Progeny_DOWN",names(namesOnly))],
              order.by = "freq", 
              #group.by = "sets", 
              cutoff = 10,text.scale = 2)

UpSetR::upset(fromList(namesOnly), 
              sets = names(namesOnly)[grepl("Diet.*Male_Progeny_DOWN",names(namesOnly))],
              order.by = "freq", 
              #group.by = "sets", 
              cutoff = 10,text.scale = 2)



# dif within diet and gender Day diff ---
counter = 0
resultsListAll <- list()
namesOnly <- list()

for(i in c("CD","HFD","HSD")){
  for(j in c("Male","Female")){
    counter = counter+1
    term = paste0(j,"_",i)
    
    #dif = paste0("sex_",term)
    
    difA = paste0(term,"_Day_10_0_")
    difB = paste0(term,"_Day_50_0_")
    difC = paste0(term,"_Day_50_10_")
    
    compA=paste0(i,"_Day_0_",j,"_Progeny")
    compB=paste0(i,"_Day_10_",j,"_Progeny")
    compC=paste0(i,"_Day_50_",j,"_Progeny")
    
    res_tmp <- results(de_seq_all, contrast = c("Condition",compB,compA))
    resultsListAll[[difA]] <- res_tmp
    namesOnly[[paste0(difA,"_UP")]] <- rownames(subset(
      results(de_seq_all, contrast = c("Condition",compB,compA)),
      padj < 0.05 & log2FoldChange > 1))
    namesOnly[[paste0(difA,"_DOWN")]] <- rownames(subset(
      results(de_seq_all, contrast = c("Condition",compB,compA)),
      padj < 0.05 & log2FoldChange < -1))
    namesOnly[[paste0(difA,"_DEG")]] <- rownames(subset(
      results(de_seq_all, contrast = c("Condition",compB,compA)),
      padj < 0.05 & abs(log2FoldChange) > 1))
    
    res_tmp <- results(de_seq_all, contrast = c("Condition",compC,compA))
    resultsListAll[[difB]] <- res_tmp
    namesOnly[[paste0(difB,"_UP")]] <- rownames(subset(
      results(de_seq_all, contrast = c("Condition",compC,compA)),
      padj < 0.05 & log2FoldChange > 1))
    namesOnly[[paste0(difB,"_DOWN")]] <- rownames(subset(
      results(de_seq_all, contrast = c("Condition",compC,compA)),
      padj < 0.05 & log2FoldChange < -1))
    namesOnly[[paste0(difB,"_DEG")]] <- rownames(subset(
      results(de_seq_all, contrast = c("Condition",compC,compA)),
      padj < 0.05 & abs(log2FoldChange) > 1))
    
    res_tmp <- results(de_seq_all, contrast = c("Condition",compC,compB))
    resultsListAll[[difC]] <- res_tmp
    namesOnly[[paste0(difC,"_UP")]] <- rownames(subset(
      results(de_seq_all, contrast = c("Condition",compC,compB)),
      padj < 0.05 & log2FoldChange > 1))
    namesOnly[[paste0(difC,"_DOWN")]] <- rownames(subset(
      results(de_seq_all, contrast = c("Condition",compC,compB)),
      padj < 0.05 & log2FoldChange < -1))
    namesOnly[[paste0(difC,"_DEG")]] <- rownames(subset(
      results(de_seq_all, contrast = c("Condition",compC,compB)),
      padj < 0.05 & abs(log2FoldChange) > 1))

    print(paste0("Done: ",counter,"/6"))
  }
}

names(namesOnly)


UpSetR::upset(fromList(namesOnly), 
              sets = names(namesOnly)[grepl("^Male.*DEG",names(namesOnly))],
              order.by = "degree", 
              nintersects = 1,
              #group.by = "sets", 
              cutoff = 10,
              text.scale = 1,
              )
intersectAll <- Reduce(intersect,
       namesOnly[names(namesOnly)[grepl("^Male.*DEG",names(namesOnly))]])
rowData(de_seq_all)[intersectAll,"SYMBOL"]


UpSetR::upset(fromList(namesOnly), 
              sets = names(namesOnly)[grepl("^Female.*DEG",names(namesOnly))],
              order.by = "degree", 
              nintersects = 1,
              #group.by = "sets", 
              cutoff = 10,
              text.scale = 1,
)

# Now idea to get all DEG within Gender,
# within Diet between TIme points

OverlapAll_male <- Reduce(unique,
                       namesOnly[names(namesOnly)[grepl("^Male.*DEG",names(namesOnly))]])

length(OverlapAll_male)

OverlapAll_female <- Reduce(unique,
                          namesOnly[names(namesOnly)[grepl("^Female.*DEG",names(namesOnly))]])

length(OverlapAll_female)



# Do heatmap on intresting set ----
ofInterest <- intersectAll

subsetData <- de_seq_all_vst[rownames(de_seq_all_vst) %in% ofInterest,]
subsetData_df <- as.data.frame(t(assay(subsetData)))
subsetData_df <- cbind(subsetData_df,
                       colData(de_seq_all_vst)
                       [rownames(subsetData_df),c("Diet","Gender","Day","Merged")])
#colnames(subsetData_df)[ncol(subsetData_df)] <- c("Diet","Gender","Day","Merged")
# mean_per_condition = aggregate(x = subsetData_df[,-ncol(subsetData_df)],
#                                by = list(as.factor(subsetData_df[,"Merged"])),
#                                FUN = mean)
# 
# 
# rownames(mean_per_condition) <- mean_per_condition$Group.1
# mean_per_condition$Group.1 <- NULL
# colnames(mean_per_condition) <- rowData(CoCena_Input_HC)[rownames(rowData(CoCena_Input_HC)) %in% colnames(mean_per_condition),"SYMBOL"]
# z_values <- mean_per_condition
z_values <- subsetData_df
#colnames(z_values) <- colnames(mean_per_condition)

#get Colors
nBins_each <- 5
myColor <- c(colorRampPalette(c("#0077b6", "#d9effa"))(nBins_each),
             "white",
             colorRampPalette(c("#ffe8e8","#D62828"))(nBins_each)
)

annoCol <- list(group = colorTheme)
col_annp_df <- data.frame(rownames(z_values),z_values[,c("Diet","Gender","Day","Merged")])
rownames(col_annp_df) = col_annp_df[,1]
col_annp_df[,1]<-NULL
str(col_annp_df)
col_annp_df$Merged <- NULL

annoCol <- list(Diet = colorTheme[c("CD","HFD","HSD")],
                Gender = colorTheme[c("Male","Female")],
                Day = colorTheme[c("0","10","50")])

library(pheatmap)
toPlot <- t(z_values[,!colnames(z_values)%in%c("Diet","Gender","Day","Merged")])
heatmap <- pheatmap(toPlot,
                    color = myColor,
                    scale="row",
                    fontsize_number = 15,
                    #cellheight = 10,
                    #cellwidth = 10,
                    annotation_colors = annoCol,
                    annotation_col  = col_annp_df,
                    angle_col = 90,
                    cluster_rows = T
)

maleOnly <- toPlot[,grepl("M",colnames(toPlot))]
heatmap <- pheatmap(maleOnly,
                    color = myColor,
                    scale="row",
                    fontsize_number = 15,
                    #cellheight = 10,
                    #cellwidth = 10,
                    annotation_colors = annoCol,
                    annotation_col  = col_annp_df,
                    angle_col = 90,
                    cluster_rows = T
)


## Work with DEG Lists ----
names(namesOnly_withinGender)

# get all UP in HSD over all Days in Female
library(UpSetR)
clipr::write_clip(c(namesOnly_withinGender[[selectedSets[5]]],namesOnly_withinGender[[selectedSets[6]]]))

selectedSets <- names(namesOnly_withinGender)[grepl(".*Day_0_.*Female.",names(namesOnly_withinGender))]
print(selectedSets)

UpSetR::upset(fromList(namesOnly_withinGender), 
              sets = selectedSets,
              nsets = length(selectedSets),order.by = "degree"
)

selectedSets <- names(namesOnly_withinGender)[grepl(".*Day_10_.*Female",names(namesOnly_withinGender))]
print(selectedSets)

UpSetR::upset(fromList(namesOnly_withinGender), 
              sets = selectedSets,
              nsets = length(selectedSets),order.by = "degree"
)

selectedSets <- names(namesOnly_withinGender)[grepl("^Diet.*Day_50_.*Female.*DOWN",names(namesOnly_withinGender))]
print(selectedSets)

UpSetR::upset(fromList(namesOnly_withinGender), 
              sets = selectedSets,
              nsets = length(selectedSets),order.by = "degree"
)


# get all of the sets 
# List of intersections and their entries
library(dplyr)
getSets <- function(selectedList){
  df2 <- data.frame(gene=unique(unlist(selectedList)))
  df1 <- lapply(selectedList,function(x){
    data.frame(gene = x)
  }) %>% 
    bind_rows(.id = "path")
  
  df_int <- lapply(df2$gene,function(x){
    # pull the name of the intersections
    intersection <- df1 %>% 
      dplyr::filter(gene==x) %>% 
      arrange(path) %>% 
      pull("path") %>% 
      paste0(collapse = "|")
    
    # build the dataframe
    data.frame(gene = x,int = intersection)
  }) %>% 
    bind_rows()
}
DEG_list <- list()
list2Excel <- list()
for(Day in c("0","10","50")){
  selectedSets_Day <- names(namesOnly_withinGender)[grepl(paste0("^Diet.*Day_",Day),names(namesOnly_withinGender))]
  for(sex in c("Female","Male")){
    selectedSets_Day_Gender <- names(namesOnly_withinGender)[grepl(paste0("^Diet.*Day_",Day,".*",sex),names(namesOnly_withinGender))]
    new_names <- unique(gsub("_..?.?.$","",selectedSets_Day_Gender))
    
    for(i in 1:length(new_names)){
      DEG_list[[new_names[i]]] <- c(unlist(namesOnly_withinGender[grepl(paste0("^",new_names[i]),names(namesOnly_withinGender))],use.names = F))
      
    }
    
    new_names <- unique(gsub("_..?.?.$","",selectedSets_Day_Gender))
    tmp <- getSets(DEG_list[new_names])
    summary <- tmp %>% 
      group_by(int) %>% 
      summarise(n=n()) %>% 
      arrange(desc(n))
    
    summary$SetName <- paste0("Set_",LETTERS[1:nrow(summary)],"_",summary$n)
    splittedForExcel <- split(tmp$gene, tmp$int)
    # rename for Excel Purposes
    for(i in 1:length(splittedForExcel)){
      names(splittedForExcel)[i] <- c(summary[summary$int%in%names(splittedForExcel)[i],"SetName"])
    }
    
    list2Excel[[paste0("Overview_",Day,"_",sex,"_DEG")]] <- as.data.frame(summary)
    
    list2Excel[paste0("singleSets_",Day,"_",sex,"_DEG_",names(splittedForExcel))] <- splittedForExcel
    
    
    for(direction in c("UP","DOWN")){
      selectedSets_Day_Gender_dir <- names(namesOnly_withinGender)[grepl(paste0("^Diet.*Day_",Day,".*",sex,".*",direction),names(namesOnly_withinGender))]
    
      tmp <- getSets(namesOnly_withinGender[selectedSets_Day_Gender_dir])
      summary <- tmp %>% 
        group_by(int) %>% 
        summarise(n=n()) %>% 
        arrange(desc(n))
      
      summary$SetName <- paste0("Set_",LETTERS[1:nrow(summary)],"_",summary$n)
      splittedForExcel <- split(tmp$gene, tmp$int)
      # rename for Excel Purposes
      for(i in 1:length(splittedForExcel)){
        names(splittedForExcel)[i] <- c(summary[summary$int%in%names(splittedForExcel)[i],"SetName"])
      }
      # construct excel
      # overview above (to resemble UpsetR)
      # for each set a sheet
      # have in their UP , DOWN and DEG
      
      list2Excel[[paste0("Overview_",Day,"_",sex,"_",direction)]] <- summary
      list2Excel[paste0("singleSets_",Day,"_",sex,"_",direction,"_",names(splittedForExcel))] <- splittedForExcel
    }
  }
  
  new_names <- unique(gsub("_..?.?.$","",selectedSets_Day))
  tmp <- getSets(DEG_list[new_names])
  summary <- tmp %>% 
    group_by(int) %>% 
    summarise(n=n()) %>% 
    arrange(desc(n))
  
  summary$SetName <- paste0("Set_",LETTERS[1:nrow(summary)],"_",summary$n)
  splittedForExcel <- split(tmp$gene, tmp$int)
  # rename for Excel Purposes
  for(i in 1:length(splittedForExcel)){
    names(splittedForExcel)[i] <- c(summary[summary$int%in%names(splittedForExcel)[i],"SetName"])
  }

  list2Excel[[paste0("Overview_",Day,"_DEG")]] <- summary
  list2Excel[paste0("singleSets_",Day,"_DEG","_",names(splittedForExcel))] <- splittedForExcel
}


#for convenience add gene names to IDs
rowAnno <- as.data.frame(rowData(de_seq_all))

names(resultsListAll)
names(resultsListAll_withinGender)

list2Excel_setsOnly <- list()

for(i in names(resultsListAll_withinGender)){
  DEG_tmp <- subset(resultsListAll_withinGender[[i]],padj < 0.05)
  up_tmp <- subset(DEG_tmp, log2FoldChange > 0)
  down_tmp <- subset(DEG_tmp, log2FoldChange < 0)
  setName <- gsub("pAdjOnly_","",i)
  
  list2Excel_setsOnly[[paste0(setName,"_DEG")]] <- data.frame(
    GeneID = rownames(DEG_tmp),
    GeneName=rowAnno[rownames(DEG_tmp),c("SYMBOL")],
    DEG_tmp
  )
                                                            
}


for(sex in c("Female","Male")){
  tmp_list <- list2Excel_setsOnly[
    sort(names(list2Excel_setsOnly)[grepl(sex,names(list2Excel_setsOnly))])]
  names(tmp_list) <- gsub("Diet_|_Progeny","",names(tmp_list))
  openxlsx::write.xlsx(tmp_list,
                       paste0("./DEG_lists_withinDay_DietDiff_",sex,".xlsx"))
}




for(i in names(list2Excel)[grepl("singleSets",names(list2Excel))]){
  list2Excel[[i]] <- data.frame(ID=list2Excel[[i]],
                                GeneName=rowAnno[list2Excel[[i]],c("SYMBOL")],
                                pAdj=rowAnno[list2Excel[[i]],c("SYMBOL")])
}

# For seyhmus no datasets overlap but sets all


# For excel split the entire list into / per Day and all Overview_*day*_DEG as special lists

library(openxlsx)

for(Day in c("0","10","50")){
  tmp_list <- list2Excel[sort(names(list2Excel)[grepl(paste0("_",Day,"_DEG"),names(list2Excel))])]
  names(tmp_list) <- gsub("singleSets_","",names(tmp_list))
  openxlsx::write.xlsx(tmp_list,
                       paste0("../results/Day",Day,"_DEG_bothSex.xlsx"))
  for(sex in c("Female","Male")){
    tmp_list <- list2Excel[sort(names(list2Excel)[grepl(paste0("_",Day,"_",sex),names(list2Excel))])]
    names(tmp_list) <- gsub("singleSets_","",names(tmp_list))
    openxlsx::write.xlsx(tmp_list,
                         paste0("../results/Day",Day,"_",sex,".xlsx"))
  }
}



source("./utils/doORA.R")
  
library(msigdbr)
library(clusterProfiler)
library(org.Dm.eg.db)
  
union_all_genes_in_network <- rownames(rowAnno)


universe_entrez <- clusterProfiler::bitr(union_all_genes_in_network,
                                         fromType = "FLYBASE",
                                         toType = "ENTREZID",
                                         OrgDb = org.Dm.eg.db)$ENTREZID


ORA_cluster_results <- list()
setwd("../results/EnrichmentResults")
setsOf_Interes <- names(list2Excel)[grepl("singleSets_",names(list2Excel))]
for(i in setsOf_Interes){ 
  gene_set <- list2Excel[[i]]
  # assumes always ENSEMBL IDs
  geneSetChoice_tranlsated <- clusterProfiler::bitr(gene_set$ID,
                                                    fromType="FLYBASE",
                                                    toType="ENTREZID",
                                                    OrgDb=org.Dm.eg.db)$ENTREZID
  fileName = i
  ORA_cluster_results[[i]] <- doOra(
    geneSetChoice_tranlsated, # ENSEBML
    type=c("GO","KEGG","HALLMARK"),
    levelGOTerms=6,
    universe_entrez,
    filename = fileName # will be ORA_[filename][type].png
  )
  #ORA_cluster_results[[i]][["GeneSet"]] <- clusterGenes[[i]]
}

setwd("../../")

## Save for html file of data overview
saveRDS(ORA_cluster_results,"Transcriptomics_ORA_CoCena_results_KC.rds")


# do Diet vs CD comparisons only ----
list2Excel_setsOnly

# get all DEG from Day0 HFD vs CD, Day 10 , and Day 50

term <- "Diet_HFD_CD.*Male" 
term <- "Diet_HSD_CD.*Female"

DataOfInterest<-names(list2Excel_setsOnly)[grepl(term,names(list2Excel_setsOnly))]

library(ggvenn)
# Make venn with all DEG from HFD vs CD 

subsetList2plot <- lapply(list2Excel_setsOnly[DataOfInterest],function(x){x[,"GeneName"]})
names(subsetList2plot) <- gsub(".*(Day_[0-9]+).*", "\\1", names(subsetList2plot))
colorTheme_wDays <- colorTheme
names(colorTheme_wDays)[6:8] <- c("Day_0","Day_10","Day_50")

test_colors <- unlist(colorTheme_wDays[names(subsetList2plot)],use.names = F)
names(test_colors) <- NULL

ggvenn(subsetList2plot,
       fill_color = test_colors) +
  ggtitle(term)

# get genes that are DE in at least two conditions

genesInAtLeast2 <- unlist(subsetList2plot)[duplicated(unlist(subsetList2plot))]
genesInAtLeast2_ID <- rownames(rowAnno[rowAnno$SYMBOL %in% genesInAtLeast2,])

# split in Up and Down
library(rlist)
tmp_list <- rlist::list.rbind(list2Excel_setsOnly[DataOfInterest])
tmd_df <- tmp_list[tmp_list$GeneID %in% genesInAtLeast2_ID,]
tmd_df$LFC_cat <- "UPandDOWN"
for(i in unique(tmd_df$GeneID)){
  idx <- which(tmd_df$GeneID == i)
  if(all(tmd_df[idx,"log2FoldChange"] > 0)){
    tmd_df[idx,"LFC_cat"] <- "UP"
  }
  if(all(tmd_df[idx,"log2FoldChange"] < 0)){
    tmd_df[idx,"LFC_cat"] <- "DOWN"
  }
}

table(tmd_df$GeneID,tmd_df$LFC_cat)

uniqueSets_UP <- unique(tmd_df[tmd_df$LFC_cat == "UP","GeneID"])
uniqueSets_DOWN <- unique(tmd_df[tmd_df$LFC_cat == "DOWN","GeneID"])
uniqueSets_UPandDOWN <- unique(tmd_df[tmd_df$LFC_cat == "UPandDOWN","GeneID"])

print(paste0("UP: ",length(uniqueSets_UP)))
print(paste0("DOWN: ",length(uniqueSets_DOWN)))
print(paste0("UPandDOWN: ",length(uniqueSets_UPandDOWN)))

listToEnrich <- list(unqiueSets_UP = uniqueSets_UP,
                     uniqueSets_DOWN = uniqueSets_DOWN,
                     uniqueSets_UPandDOWN = uniqueSets_UPandDOWN)

ORA_cluster_results_singleSets <- list()

union_all_genes_in_network <- rownames(rowAnno)


universe_entrez <- clusterProfiler::bitr(union_all_genes_in_network,
                                         fromType = "FLYBASE",
                                         toType = "ENTREZID",
                                         OrgDb = org.Dm.eg.db)$ENTREZID
length(universe_entrez)
#setwd("../results/EnrichmentResults_pics/oneComparisonOnly/")
for(i in names(listToEnrich)){ 
  gene_set <- listToEnrich[[i]]
  # assumes always ENSEMBL IDs
  geneSetChoice_tranlsated <- clusterProfiler::bitr(gene_set,
                                                    fromType="FLYBASE",
                                                    toType="ENTREZID",
                                                    OrgDb=org.Dm.eg.db)$ENTREZID
  fileName = paste0(term,"_",i)
  ORA_cluster_results_singleSets[[i]] <- doOra(
    geneSetChoice_tranlsated, # ENSEBML
    type=c("GO","KEGG","HALLMARK"),
    levelGOTerms=6,
    universe_entrez,
    filename = fileName # will be ORA_[filename][type].png
  )
  #ORA_cluster_results[[i]][["GeneSet"]] <- clusterGenes[[i]]
}


# Do Vis of 3 gens ----
geneOfinterest <- c("Rel","Myd88","upd3")
geneOfinterest_ID <-rownames(rowData(de_seq_all)[rowData(de_seq_all)$SYMBOL %in% geneOfinterest,])

prepData <- cbind(t(assay(de_seq_all_vst)[geneOfinterest_ID,]),
                  colData(de_seq_all_vst)[colnames(assay(de_seq_all_vst)[geneOfinterest_ID,]),])

colnames(prepData)[1:3] <- c("Rel","Myd88","upd3")
prepData <- as.data.frame(prepData)
prepData$DietDay <- paste0(prepData$Diet,"_",prepData$Day)

colorTheme <- c(colorTheme,c("#dde5b6", "#7cb518", "#31572c"))
names(colorTheme)[9:11] <- c("0","10","50")
library(ggpubr)
ggplot(prepData, 
       aes(x=DietDay, 
           y=upd3,
           color=Day)) + 
  scale_color_manual(values=colorTheme[c("0","10","50")])+
  geom_point(color="black")+
  geom_boxplot(alpha=0.7)+
  theme_bw()+
  facet_grid(~Gender)


ggplot(prepData, 
       aes(x = Diet, 
           y = Myd88,
           color = Day)) + 
  scale_color_manual(values=colorTheme[c("0","10","50")])+
  geom_point(color="black")+
  geom_boxplot(alpha=0.7)+
  theme_bw()+
  facet_grid(~Day+Gender)


resultsListAll_withinGender$pAdjOnly_Diet_HFD_CD_Day_50_Male_Progeny[geneOfinterest_ID,]


# Change DESeq formula ----
library(DESeq2)
ddsSE_split <- DESeqDataSet(se_fly_all, design = ~ Diet + Day + Diet:Day)


ddsSE_split$Day <- relevel(ddsSE_split$Day, ref = "0")
ddsSE_split$Diet <- relevel(ddsSE_split$Diet, ref = "CD")
ddsSE_split$Gender <- relevel(ddsSE_split$Gender, ref = "Male")

ddsSE_split <- preprocessing(ddsSE_split,10,
                       protCodingOnly=T,
                       removeConstRows=T,
                       filterPerSample=T)


ddsSE_split <- DESeq(ddsSE_split) 
resultsNames(ddsSE_split)
summary(results(ddsSE_split))

# The key point to remember about designs with interaction terms is that, 
# unlike for a design ~genotype + condition, 
# where the condition effect represents the overall effect controlling for differences 
# due to genotype, 
# by adding genotype:condition, 
# the main condition effect only represents the effect of condition for the reference level.
# The interaction terms genotypeII.conditionB and genotypeIII.conditionB 
# give the difference between the condition effect for a given genotype 
# and the condition effect for the reference genotype.

mainEffect_HFD <- results(ddsSE_split, contrast=c("Diet","HFD","CD"))
summary(mainEffect_HFD)
mainEffect_HSD <- results(ddsSE_split, contrast=c("Diet","HSD","CD"))
summary(mainEffect_HSD)

# the condition effect for Diet HFD
# this is the main effect *plus* the interaction term (Diet:Day)
# (the extra condition effect in HFD compared to CD).
main_interaction_HFD_10 <- results(ddsSE_split, contrast=list( c("Diet_HFD_vs_CD","Gender_Female_vs_Male","DietHFD.Day10") ))
summary(main_interaction_HFD_10)
main_interaction_HFD_50 <- results(ddsSE_split, contrast=list( c("Diet_HFD_vs_CD","Gender_Female_vs_Male","DietHFD.Day50") ))
summary(main_interaction_HFD_50)

main_interaction_HSD_10 <- results(ddsSE_split, contrast=list( c("Diet_HSD_vs_CD","Gender_Female_vs_Male","DietHSD.Day10") ))
summary(main_interaction_HSD_10)
main_interaction_HSD_50 <- results(ddsSE_split, contrast=list( c("Diet_HSD_vs_CD","Gender_Female_vs_Male","DietHSD.Day50") ))
summary(main_interaction_HSD_50)


main_interaction_HFD_10 <- results(ddsSE_split, contrast=list( c("Diet_HFD_vs_CD","DietHFD.Day10") ))
summary(main_interaction_HFD_10)
main_interaction_HFD_50 <- results(ddsSE_split, contrast=list( c("Diet_HFD_vs_CD","DietHFD.Day50") ))
summary(main_interaction_HFD_50)

main_interaction_HSD_10 <- results(ddsSE_split, contrast=list( c("Diet_HSD_vs_CD","DietHSD.Day10") ))
summary(main_interaction_HSD_10)
main_interaction_HSD_50 <- results(ddsSE_split, contrast=list( c("Diet_HSD_vs_CD","DietHSD.Day50") ))
summary(main_interaction_HSD_50)
library(ggVennDiagram)
ggVennDiagram(list(HFD_Day_50_UP =rownames(subset(main_interaction_HFD_50,padj<0.1 & log2FoldChange >0)),
            HFD_Day_50_DOWN = rownames(subset(main_interaction_HFD_50,padj<0.1 & log2FoldChange <0)),
            HSD_Day_50_UP = rownames(subset(main_interaction_HSD_50,padj<0.1 & log2FoldChange >0)),
            HSD_Day_50_DOWN = rownames(subset(main_interaction_HSD_50,padj<0.1 & log2FoldChange <0))),
       label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)

# getting this 839 genes
sharedUpregulated <- intersect(rownames(subset(main_interaction_HFD_50,padj<0.1 & log2FoldChange >0)),
         rownames(subset(main_interaction_HSD_50,padj<0.1 & log2FoldChange >0)))

HFD_unique <- setdiff(rownames(subset(main_interaction_HFD_50,padj<0.1 & log2FoldChange >0)),
                      sharedUpregulated)

HSD_unique <- setdiff(rownames(subset(main_interaction_HSD_50,padj<0.1 & log2FoldChange >0)),
                      sharedUpregulated)

# assumes always ENSEMBL IDs
geneSetChoice_tranlsated <- clusterProfiler::bitr(HSD_unique,
                                                  fromType="FLYBASE",
                                                  toType="ENTREZID",
                                                  OrgDb=org.Dm.eg.db)$ENTREZID


ORA_cluster_results_interaction<- doOra(
  geneSetChoice_tranlsated, # ENSEBML
  type=c("GO","KEGG","HALLMARK"),
  levelGOTerms=2,
  universe_entrez,
  filename = "UniqueUpregulated_day50_HSD" # will be ORA_[filename][type].png
)


# the interaction term for Diet effect in day 0 vs day 10.
# this tests if the Day effect is different in HFD compared to CD
Day10_effect_onDiet_HFD <- results(ddsSE_split, name="DietHFD.Day10")
Day50_effect_onDiet_HFD <- results(ddsSE_split, name="DietHFD.Day50")
Day10_effect_onDiet_HSD <- results(ddsSE_split, name="DietHSD.Day10")
Day50_effect_onDiet_HSD <- results(ddsSE_split, name="DietHSD.Day50")

summary(Day10_effect_onDiet_HFD)

list_interactionEffects_down <- list(Day10_effect_onDiet_HFD =rownames(subset(Day10_effect_onDiet_HFD,padj<0.1 & log2FoldChange <0)),
                                 Day50_effect_onDiet_HFD = rownames(subset(Day50_effect_onDiet_HFD,padj<0.1 & log2FoldChange <0)),
                                 Day10_effect_onDiet_HSD = rownames(subset(Day10_effect_onDiet_HSD,padj<0.1 & log2FoldChange <0)),
                                 Day50_effect_onDiet_HSD = rownames(subset(Day50_effect_onDiet_HSD,padj<0.1 & log2FoldChange <0)))
list_interactionEffects_up <- list(Day10_effect_onDiet_HFD =rownames(subset(Day10_effect_onDiet_HFD,padj<0.1 & log2FoldChange >0)),
                                     Day50_effect_onDiet_HFD = rownames(subset(Day50_effect_onDiet_HFD,padj<0.1 & log2FoldChange >0)),
                                     Day10_effect_onDiet_HSD = rownames(subset(Day10_effect_onDiet_HSD,padj<0.1 & log2FoldChange >0)),
                                     Day50_effect_onDiet_HSD = rownames(subset(Day50_effect_onDiet_HSD,padj<0.1 & log2FoldChange >0)))


ggVennDiagram(list_interactionEffects_down,
              label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)


HFD_unique_UP <- setdiff(unlist(list_interactionEffects_up["Day50_effect_onDiet_HFD"]),
                         unlist(list_interactionEffects_up[names(list_interactionEffects_up)!="Day50_effect_onDiet_HFD"])
)
length(HFD_unique_UP)

HSD_unique_UP <- setdiff(unlist(list_interactionEffects_up["Day50_effect_onDiet_HSD"]),
                         unlist(list_interactionEffects_up[names(list_interactionEffects_up)!="Day50_effect_onDiet_HSD"])
)

shared_interaction <- Reduce(intersect,list_interactionEffects_up)
tmp <- clusterProfiler::bitr(shared_interaction,
                             fromType="FLYBASE",
                             toType="SYMBOL",
                             OrgDb="org.Dm.eg.db")$SYMBOL
length(shared_interaction)


# assumes always ENSEMBL IDs
geneSetChoice_tranlsated <- clusterProfiler::bitr(shared_interaction,
                                                  fromType="FLYBASE",
                                                  toType="ENTREZID",
                                                  OrgDb=org.Dm.eg.db)$ENTREZID


ORA_cluster_results_interaction<- doOra(
  geneSetChoice_tranlsated, # ENSEBML
  type=c("GO","KEGG","HALLMARK"),
  levelGOTerms=2,
  universe_entrez,
  filename = "UniqueUpregulated_interaction_alldays_shared" # will be ORA_[filename][type].png
)


# test 50 vs 10
HSD_50vs10 <- results(ddsSE_split, contrast = list("DietHSD.Day50", "DietHSD.Day10"))
HSD_10vs0 <- results(ddsSE_split, contrast = list("DietHSD.Day10"))
HSD_50vs0 <- results(ddsSE_split, contrast = list("DietHSD.Day50"))

HSD_interactionOverTime_UP <- list(
  HSD_50vs10_UP = rownames(subset(HSD_50vs10,padj<0.1 & log2FoldChange >0)),
 # HSD_50vs10_DOWN = rownames(subset(HSD_50vs10,padj<0.1,log2FoldChange <0)),
  HSD_10vs0_UP = rownames(subset(HSD_10vs0,padj<0.1&log2FoldChange >0)),
  #HSD_10vs0_DOWN = rownames(subset(HSD_10vs0,padj<0.1,log2FoldChange <0)),
  HSD_50vs0_UP = rownames(subset(HSD_50vs0,padj<0.1&log2FoldChange >0))
  #HSD_50vs0_DOWN = rownames(subset(HSD_50vs0,padj<0.1,log2FoldChange <0))
)

HSD_interactionOverTime_DOWN <- list(
  #HSD_50vs10_UP = rownames(subset(HSD_50vs10,padj<0.1,log2FoldChange >0)),
  HSD_50vs10_DOWN = rownames(subset(HSD_50vs10,padj<0.1&log2FoldChange <0)),
  #HSD_10vs0_UP = rownames(subset(HSD_10vs0,padj<0.1,log2FoldChange >0)),
  HSD_10vs0_DOWN = rownames(subset(HSD_10vs0,padj<0.1&log2FoldChange <0)),
  #HSD_50vs0_UP = rownames(subset(HSD_50vs0,padj<0.1,log2FoldChange >0)),
  HSD_50vs0_DOWN = rownames(subset(HSD_50vs0,padj<0.1&log2FoldChange <0))
)

library(ggVennDiagram)
library(ggplot2)
ggVennDiagram(HSD_interactionOverTime_UP,
              label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)

ggVennDiagram(HSD_interactionOverTime_DOWN,
              label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)



# test 50 vs 10
HFD_50vs10 <- results(ddsSE_split, contrast = list("DietHFD.Day50", "DietHFD.Day10"))
HFD_10vs0 <- results(ddsSE_split, contrast = list("DietHFD.Day10"))
HFD_50vs0 <- results(ddsSE_split, contrast = list("DietHFD.Day50"))

HFD_interactionOverTime_UP <- list(
  HFD_50vs10_UP = rownames(subset(HFD_50vs10,padj<0.1 & log2FoldChange >0)),
  # HSD_50vs10_DOWN = rownames(subset(HSD_50vs10,padj<0.1,log2FoldChange <0)),
  HFD_10vs0_UP = rownames(subset(HFD_10vs0,padj<0.1&log2FoldChange >0)),
  #HSD_10vs0_DOWN = rownames(subset(HSD_10vs0,padj<0.1,log2FoldChange <0)),
  HFD_50vs0_UP = rownames(subset(HFD_50vs0,padj<0.1&log2FoldChange >0))
  #HSD_50vs0_DOWN = rownames(subset(HSD_50vs0,padj<0.1,log2FoldChange <0))
)

HFD_interactionOverTime_DOWN <- list(
  #HSD_50vs10_UP = rownames(subset(HSD_50vs10,padj<0.1,log2FoldChange >0)),
  HFD_50vs10_DOWN = rownames(subset(HFD_50vs10,padj<0.1&log2FoldChange <0)),
  #HSD_10vs0_UP = rownames(subset(HSD_10vs0,padj<0.1,log2FoldChange >0)),
  HFD_10vs0_DOWN = rownames(subset(HFD_10vs0,padj<0.1&log2FoldChange <0)),
  #HSD_50vs0_UP = rownames(subset(HSD_50vs0,padj<0.1,log2FoldChange >0)),
  HFD_50vs0_DOWN = rownames(subset(HFD_50vs0,padj<0.1&log2FoldChange <0))
)

library(ggVennDiagram)
library(ggplot2)
ggVennDiagram(HFD_interactionOverTime_UP,
              label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)

ggVennDiagram(HFD_interactionOverTime_DOWN,
              label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)

# this tests if the Day effect is different in Diet HFD vs HSD.
Day10_effects_vsDiets <- results(ddsSE_split, contrast=list("DietHFD.Day10", "DietHSD.Day10"))
summary(Day10_effects_vsDiets)

resultsNames(ddsSE_split)
Gender_Female_vs_Male <- results(ddsSE_split, contrast=list("Gender_Female_vs_Male"))

Diet_HFD_vs_CD <- results(ddsSE_split, contrast=list("Diet_HFD_vs_CD"))
summary(Diet_HFD_vs_CD)
# CHeck expression from Seyhmus provided sets ----
library(readxl)
path <- "../data/geneLists_Seyhmus/Immune_related_genes_from_GO_Biological_Process_2018.xls"
Go_list_provided <- list()
for(i in readxl::excel_sheets(path)){
  Go_list_provided[[  gsub(" ","_",i)]] <- 
    read_excel(path,sheet=i)
}
lapply(Go_list_provided,nrow)

path <- "../data/geneLists_Seyhmus/Longevity_regulating_pathway_from_KEG20192.xls"
Kegg_list_provided <- list()
for(i in readxl::excel_sheets(path)){
  Kegg_list_provided[[  gsub(" ","_",i)]] <- 
    read_excel(path,sheet=i)
}

path <- "../data/geneLists_Seyhmus/neurodegenerative diseases from human_diseases_FlyBase_2017.xls"
Flybase_list_provided <- list()
for(i in readxl::excel_sheets(path)){
  Flybase_list_provided[[  gsub(" ","_",i)]] <- 
    read_excel(path,sheet=i)
}

path <- "../data/geneLists_Seyhmus/Selected_GENES_for_heatmap_drosophila.xlsx"
Selected_list_provided <- list()
for(i in readxl::excel_sheets(path)){
  Selected_list_provided[[  gsub(" ","_",i)]] <- 
    read_excel(path,sheet=i)
}

# plot the given genes in a pheatmap to check upon clustering
dfToPlot <- as.matrix(assay(de_seq_all_vst)) #maybe we want vst data?
table(apply(dfToPlot,1,sd)==0)
rownames(dfToPlot) <- rowData(de_seq_all_vst)$SYMBOL
translatedSubset <- unlist(Go_list_provided$GO_0045089,use.names = F)
# colnames(Go_list_provided$GO_0045089)
translatedSubset <- unlist(Flybase_list_provided$neurodegenerative_diseases_from,use.names = F)
translatedSubset <- unlist(Selected_list_provided$Is_r_genes,use.names = F)
translatedSubset <- unlist(Kegg_list_provided$Longevity_regulating_pathway_fr,use.names = F)


#check which are still present
keep <- translatedSubset[translatedSubset %in% rownames(dfToPlot)]
colData(de_seq_all_vst)$Group <- as.factor(colData(de_seq_all_vst)$Group)

# separate by gender the heatmaps
samples2keep <- colData(de_seq_all_vst)$Gender %in% "Male" 
samples2keep <- colData(de_seq_all_vst)$Gender %in% "Female" 
#samples2keep <- colData(de_seq_all_vst)$Gender %in% "Female" & colData(de_seq_all_vst)$Day %in% 50
#names(annoCol[[3]]) <- gsub("Day","",names(annoCol[[3]]))

#annoCol_tmp <- annoCol
#annoCol_tmp$Day <- annoCol_tmp$Day[-4]

names(which(apply(dfToPlot[keep,samples2keep],1,sd)==0))

pheatmap::pheatmap(dfToPlot[keep[!apply(dfToPlot[keep,samples2keep],1,sd)==0],samples2keep],
                   color = myColor,
                   scale = "row",
                   fontsize_number = 15,
                   cluster_rows = T,
                   cluster_cols = T,
                   show_rownames = T,
                   show_colnames = F,
                   annotation_colors = annoCol_tmp,
                   annotation_col  = col_annp_df)


# change DESeq formula separate for male & female ----
ddsSE_preFilter_male <- ddsSE_male

design(ddsSE_male) <- ~ Diet + Day + Diet:Day
design(ddsSE_female) <- ~ Diet + Day + Diet:Day

## male ----
ddsSE_male$Day <- relevel(ddsSE_male$Day, ref = "0")
ddsSE_male$Diet <- relevel(ddsSE_male$Diet, ref = "CD")

ddsSE_male <- preprocessing(ddsSE_male,10,
                             protCodingOnly=T,
                             removeConstRows=T,
                             filterPerSample=T)


DE_SE_male <- DESeq(ddsSE_male) 
resultsNames(DE_SE_male)
summary(results(DE_SE_male))
table(ddsSE_male$Day,ddsSE_male$Diet)

# Main effect of Diet to reference level CD
mainEffect_HFD_male <- results(DE_SE_male, contrast=c("Diet","HFD","CD"))
summary(mainEffect_HFD_male)
mainEffect_HSD_male <- results(DE_SE_male, contrast=c("Diet","HSD","CD"))
summary(mainEffect_HSD_male)
res_list <- list(HFD_male = mainEffect_HFD_male,HSD_male=mainEffect_HSD_male)


mainEffectVis(res_list,color=annoCol)

main_interaction_HFD_10_male <- results(DE_SE_male, contrast=list( c("Diet_HFD_vs_CD","DietHFD.Day10") ))
summary(main_interaction_HFD_10_male)
main_interaction_HFD_50_male <- results(DE_SE_male, contrast=list( c("Diet_HFD_vs_CD","DietHFD.Day50") ))
summary(main_interaction_HFD_50_male)

main_interaction_HSD_10_male <- results(DE_SE_male, contrast=list( c("Diet_HSD_vs_CD","DietHSD.Day10") ))
summary(main_interaction_HSD_10_male)
main_interaction_HSD_50_male <- results(DE_SE_male, contrast=list( c("Diet_HSD_vs_CD","DietHSD.Day50") ))
summary(main_interaction_HSD_50_male)

res_list <- list(male_Day10_HFD = main_interaction_HFD_10_male,
                 male_Day50_HFD = main_interaction_HFD_50_male)

#names(annoCol$Day)<-paste0("Day",names(annoCol$Day))

mainEffectVis(res_list,color=annoCol[-1])

library(ggVennDiagram)
ggVennDiagram(list(main_interaction_HFD_10_male_UP =rownames(subset(main_interaction_HFD_10_male,padj<0.1 & log2FoldChange >0)),
                   main_interaction_HFD_10_male_DOWN =rownames(subset(main_interaction_HFD_10_male,padj<0.1 & log2FoldChange <0)),
                   main_interaction_HFD_50_male_UP = rownames(subset(main_interaction_HFD_50_male,padj<0.1 & log2FoldChange >0)),
                   main_interaction_HFD_50_male_DOWN = rownames(subset(main_interaction_HFD_50_male,padj<0.1 & log2FoldChange <0))
                   ),
              label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)

### hsd ----


main_interaction_HSD_10_male <- results(DE_SE_male, contrast=list( c("Diet_HSD_vs_CD","DietHSD.Day10") ))
summary(main_interaction_HSD_10_male)
main_interaction_HSD_50_male <- results(DE_SE_male, contrast=list( c("Diet_HSD_vs_CD","DietHSD.Day50") ))
summary(main_interaction_HSD_50_male)

res_list <- list(male_Day10_HSD = main_interaction_HSD_10_male,
                 male_Day50_HSD = main_interaction_HSD_50_male)

#names(annoCol$Day)<-paste0("Day",names(annoCol$Day))

mainEffectVis(res_list,color=annoCol[-1])

ggVennDiagram(list(main_interaction_HSD_10_male_UP =rownames(subset(main_interaction_HSD_10_male,padj<0.1 & log2FoldChange >0)),
                   main_interaction_HSD_10_male_DOWN =rownames(subset(main_interaction_HSD_10_male,padj<0.1 & log2FoldChange <0)),
                   main_interaction_HSD_50_male_UP = rownames(subset(main_interaction_HSD_50_male,padj<0.1 & log2FoldChange >0)),
                   main_interaction_HSD_50_male_DOWN = rownames(subset(main_interaction_HSD_50_male,padj<0.1 & log2FoldChange <0))
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)


# the interaction term for Diet effect in day 0 vs day 10.
# this tests if the Day effect is different in HFD compared to CD
Day10_effect_onDiet_HFD <- results(DE_SE_male, name="DietHFD.Day10")
Day50_effect_onDiet_HFD <- results(DE_SE_male, name="DietHFD.Day50")
Day10_effect_onDiet_HSD <- results(DE_SE_male, name="DietHSD.Day10")
Day50_effect_onDiet_HSD <- results(DE_SE_male, name="DietHSD.Day50")


res_list <- list(Day10_effect_onDiet_HFD = Day10_effect_onDiet_HFD,
                 Day50_effect_onDiet_HFD = Day50_effect_onDiet_HFD)

mainEffectVis(res_list,color=annoCol[-1])

ggVennDiagram(list(Day10_effect_onDiet_HFD_UP =rownames(subset(Day10_effect_onDiet_HFD,padj<0.1 & log2FoldChange >0)),
                   Day10_effect_onDiet_HFD_DOWN =rownames(subset(Day10_effect_onDiet_HFD,padj<0.1 & log2FoldChange <0)),
                   Day50_effect_onDiet_HFD_UP = rownames(subset(Day50_effect_onDiet_HFD,padj<0.1 & log2FoldChange >0)),
                   Day50_effect_onDiet_HFD_DOWN = rownames(subset(Day50_effect_onDiet_HFD,padj<0.1 & log2FoldChange <0))
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)

res_list <- list(
                 Day10_effect_onDiet_HSD = Day10_effect_onDiet_HSD,
                 Day50_effect_onDiet_HSD = Day50_effect_onDiet_HSD)

mainEffectVis(res_list,color=annoCol[-1])
ggVennDiagram(list(Day10_effect_onDiet_HSD_UP =rownames(subset(Day10_effect_onDiet_HSD,padj<0.1 & log2FoldChange >0)),
                   Day10_effect_onDiet_HSD_DOWN =rownames(subset(Day10_effect_onDiet_HSD,padj<0.1 & log2FoldChange <0)),
                   Day50_effect_onDiet_HSD_UP = rownames(subset(Day50_effect_onDiet_HSD,padj<0.1 & log2FoldChange >0)),
                   Day50_effect_onDiet_HSD_DOWN = rownames(subset(Day50_effect_onDiet_HSD,padj<0.1 & log2FoldChange <0))
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)


ggVennDiagram(list(mainEffect_HFD =rownames(subset(mainEffect_HFD_male,padj<0.1 )),
                   Day10_effect_onDiet_HFD=rownames(subset(Day10_effect_onDiet_HFD,padj<0.1)),
                   Day50_effect_onDiet_HFD = rownames(subset(Day50_effect_onDiet_HFD,padj<0.1))
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)


## Venn accross Diets ----
ggVennDiagram(list(mainEffect_HFD_UP =rownames(subset(mainEffect_HFD_male,padj<0.1 & log2FoldChange >0)),
                   mainEffect_HFD_DOWN =rownames(subset(mainEffect_HFD_male,padj<0.1 & log2FoldChange <0)),
                   mainEffect_HSD_UP = rownames(subset(mainEffect_HSD_male,padj<0.1 & log2FoldChange >0)),
                   mainEffect_HSD_DOWN = rownames(subset(mainEffect_HSD_male,padj<0.1 & log2FoldChange <0))
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)

ggVennDiagram(list(mainEffect_HFD_UP =rownames(subset(main_interaction_HFD_10_male,padj<0.1 & log2FoldChange >0)),
                   mainEffect_HFD_DOWN =rownames(subset(main_interaction_HFD_10_male,padj<0.1 & log2FoldChange <0)),
                   mainEffect_HSD_UP = rownames(subset(main_interaction_HSD_10_male,padj<0.1 & log2FoldChange >0)),
                   mainEffect_HSD_DOWN = rownames(subset(main_interaction_HSD_10_male,padj<0.1 & log2FoldChange <0))
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)

# get the shared UP and the shared down to ORA + unique HSD
sharedUpregulated <- intersect(rownames(subset(mainEffect_HSD_male,padj<0.1 & log2FoldChange >0)),
                               rownames(subset(mainEffect_HSD_male,padj<0.1 & log2FoldChange <0)))
length(sharedUpregulated)
HSD_unique <- setdiff(rownames(subset(mainEffect_HSD_male,padj<0.1 & log2FoldChange >0)),
                      c(rownames(subset(Day10_effect_onDiet_HSD,padj<0.1)),
                        rownames(subset(Day50_effect_onDiet_HSD,padj<0.1))
                        ))

length(HSD_unique)
sharedUpregulated <- HSD_unique
library(org.Dm.eg.db)
# take union of males & female identified genes
universe_entrez_all<- clusterProfiler::bitr(rownames(ddsSE_split),
                                         fromType = "FLYBASE",
                                         toType = "ENTREZID",
                                         OrgDb = org.Dm.eg.db)$ENTREZID

geneSetChoice_tranlsated <- clusterProfiler::bitr(sharedUpregulated,
                                                  fromType="FLYBASE",
                                                  toType="ENTREZID",
                                                  OrgDb=org.Dm.eg.db)$ENTREZID

length(geneSetChoice_tranlsated)
library(msigdbr)
ORA_cluster_results_interaction <- doOra(
  geneSetChoice_tranlsated, # ENSEBML
  type=c("GO","KEGG","HALLMARK"),
  levelGOTerms=8,
  universe_entrez_all,
  GoFilter=T,
  filename = "uniqHSD_UPregulated_mainEffect_noInteraction_male_" # will be ORA_[filename][type].png
)



ggVennDiagram(list(main_interaction_HFD_50_UP = rownames(subset(main_interaction_HFD_50_male,padj<0.1 & log2FoldChange >0)),
                   main_interaction_HFD_50_DOWN = rownames(subset(main_interaction_HFD_50_male,padj<0.1 & log2FoldChange <0)),
                    main_interaction_HSD_50_UP =rownames(subset(main_interaction_HSD_50_male,padj<0.1 & log2FoldChange >0)),
                   main_interaction_HSD_50_DOWN =rownames(subset(main_interaction_HSD_50_male,padj<0.1 & log2FoldChange <0))
              
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)


sharedUpregulated <- intersect(rownames(subset(main_interaction_HFD_50_male,padj<0.1 & log2FoldChange >0)),
                               rownames(subset(main_interaction_HSD_50_male,padj<0.1 & log2FoldChange >0)))
length(sharedUpregulated)
HSD_unique <- setdiff(rownames(subset(main_interaction_HSD_50_male,padj<0.1 & log2FoldChange >0)),
                      rownames(subset(main_interaction_HFD_50_male,padj<0.1 & log2FoldChange >0)))

HSD_unique <- setdiff(HSD_unique,intersect(rownames(subset(main_interaction_HSD_50_male,padj<0.1 & log2FoldChange >0)),
                                           rownames(subset(main_interaction_HFD_50_male,padj<0.1 & log2FoldChange <0))))
length(HSD_unique)
sharedUpregulated <- HSD_unique
library(org.Dm.eg.db)
# take union of males & female identified genes
universe_entrez_all <- clusterProfiler::bitr(rownames(ddsSE_split),
                                            fromType = "FLYBASE",
                                            toType = "ENTREZID",
                                            OrgDb = org.Dm.eg.db)$ENTREZID

geneSetChoice_tranlsated <- clusterProfiler::bitr(sharedUpregulated,
                                                  fromType="FLYBASE",
                                                  toType="ENTREZID",
                                                  OrgDb=org.Dm.eg.db)$ENTREZID

length(geneSetChoice_tranlsated)
library(msigdbr)
ORA_cluster_results_interaction <- doOra(
  geneSetChoice_tranlsated, # ENSEBML
  type=c("GO","KEGG","HALLMARK"),
  levelGOTerms=8,
  universe_entrez_all,
  GoFilter=T,
  filename = "uniqHSD_UPregulated_totalEffect_male" # will be ORA_[filename][type].png
)




## female ----
ddsSE_female$Day <- relevel(ddsSE_female$Day, ref = "0")
ddsSE_female$Diet <- relevel(ddsSE_female$Diet, ref = "CD")


DE_SE_female <- DESeq(ddsSE_female) 
resultsNames(DE_SE_female)
summary(results(DE_SE_female))
table(DE_SE_female$Day,DE_SE_female$Diet)


# Main effect of Diet to reference level CD
mainEffect_HFD_female <- results(DE_SE_female, contrast=c("Diet","HFD","CD"))
summary(mainEffect_HFD_female)
mainEffect_HSD_female <- results(DE_SE_female, contrast=c("Diet","HSD","CD"))
summary(mainEffect_HSD_female)
res_list <- list(HFD_female = mainEffect_HFD_female,HSD_female=mainEffect_HSD_female)


mainEffectVis(res_list,color=annoCol)

main_interaction_HFD_10_female <- results(DE_SE_female, contrast=list( c("Diet_HFD_vs_CD","DietHFD.Day10") ))
summary(main_interaction_HFD_10_female)
main_interaction_HFD_50_female <- results(DE_SE_female, contrast=list( c("Diet_HFD_vs_CD","DietHFD.Day50") ))
summary(main_interaction_HFD_50_female)

main_interaction_HSD_10_female <- results(DE_SE_female, contrast=list( c("Diet_HSD_vs_CD","DietHSD.Day10") ))
summary(main_interaction_HSD_10_female)
main_interaction_HSD_50_female <- results(DE_SE_female, contrast=list( c("Diet_HSD_vs_CD","DietHSD.Day50") ))
summary(main_interaction_HSD_50_female)

res_list <- list(female_Day10_HFD = main_interaction_HFD_10_female,
                 female_Day50_HFD = main_interaction_HFD_50_female)

#names(annoCol$Day)<-paste0("Day",names(annoCol$Day))

mainEffectVis(res_list,color=annoCol[-1])

library(ggVennDiagram)
ggVennDiagram(list(main_interaction_HFD_10_female_UP =rownames(subset(main_interaction_HFD_10_female,padj<0.1 & log2FoldChange >0)),
                   main_interaction_HFD_10_female_DOWN =rownames(subset(main_interaction_HFD_10_female,padj<0.1 & log2FoldChange <0)),
                   main_interaction_HFD_50_female_UP = rownames(subset(main_interaction_HFD_50_female,padj<0.1 & log2FoldChange >0)),
                   main_interaction_HFD_50_female_DOWN = rownames(subset(main_interaction_HFD_50_female,padj<0.1 & log2FoldChange <0))
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)

### hsd ----


main_interaction_HSD_10_female <- results(DE_SE_female, contrast=list( c("Diet_HSD_vs_CD","DietHSD.Day10") ))
summary(main_interaction_HSD_10_female)
main_interaction_HSD_50_female <- results(DE_SE_female, contrast=list( c("Diet_HSD_vs_CD","DietHSD.Day50") ))
summary(main_interaction_HSD_50_female)

res_list <- list(female_Day10_HSD = main_interaction_HSD_10_female,
                 female_Day50_HSD = main_interaction_HSD_50_female)

#names(annoCol$Day)<-paste0("Day",names(annoCol$Day))

mainEffectVis(res_list,color=annoCol[-1])

ggVennDiagram(list(main_interaction_HSD_10_female_UP =rownames(subset(main_interaction_HSD_10_female,padj<0.1 & log2FoldChange >0)),
                   main_interaction_HSD_10_female_DOWN =rownames(subset(main_interaction_HSD_10_female,padj<0.1 & log2FoldChange <0)),
                   main_interaction_HSD_50_female_UP = rownames(subset(main_interaction_HSD_50_female,padj<0.1 & log2FoldChange >0)),
                   main_interaction_HSD_50_female_DOWN = rownames(subset(main_interaction_HSD_50_female,padj<0.1 & log2FoldChange <0))
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)


# the interaction term for Diet effect in day 0 vs day 10.
# this tests if the Day effect is different in HFD compared to CD
Day10_effect_onDiet_HFD <- results(DE_SE_female, name="DietHFD.Day10")
Day50_effect_onDiet_HFD <- results(DE_SE_female, name="DietHFD.Day50")
Day10_effect_onDiet_HSD <- results(DE_SE_female, name="DietHSD.Day10")
Day50_effect_onDiet_HSD <- results(DE_SE_female, name="DietHSD.Day50")


res_list <- list(Day10_effect_onDiet_HFD = Day10_effect_onDiet_HFD,
                 Day50_effect_onDiet_HFD = Day50_effect_onDiet_HFD)

mainEffectVis(res_list,color=annoCol[-1])

ggVennDiagram(list(Day10_effect_onDiet_HFD_UP =rownames(subset(Day10_effect_onDiet_HFD,padj<0.1 & log2FoldChange >0)),
                   Day10_effect_onDiet_HFD_DOWN =rownames(subset(Day10_effect_onDiet_HFD,padj<0.1 & log2FoldChange <0)),
                   Day50_effect_onDiet_HFD_UP = rownames(subset(Day50_effect_onDiet_HFD,padj<0.1 & log2FoldChange >0)),
                   Day50_effect_onDiet_HFD_DOWN = rownames(subset(Day50_effect_onDiet_HFD,padj<0.1 & log2FoldChange <0))
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)

res_list <- list(
  Day10_effect_onDiet_HSD = Day10_effect_onDiet_HSD,
  Day50_effect_onDiet_HSD = Day50_effect_onDiet_HSD)

mainEffectVis(res_list,color=annoCol[-1])
ggVennDiagram(list(Day10_effect_onDiet_HSD_UP =rownames(subset(Day10_effect_onDiet_HSD,padj<0.1 & log2FoldChange >0)),
                   Day10_effect_onDiet_HSD_DOWN =rownames(subset(Day10_effect_onDiet_HSD,padj<0.1 & log2FoldChange <0)),
                   Day50_effect_onDiet_HSD_UP = rownames(subset(Day50_effect_onDiet_HSD,padj<0.1 & log2FoldChange >0)),
                   Day50_effect_onDiet_HSD_DOWN = rownames(subset(Day50_effect_onDiet_HSD,padj<0.1 & log2FoldChange <0))
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)





## Venn accross Diets ----
ggVennDiagram(list(mainEffect_HFD_UP =rownames(subset(mainEffect_HFD_female,padj<0.1 & log2FoldChange >0)),
                   mainEffect_HFD_DOWN =rownames(subset(mainEffect_HFD_female,padj<0.1 & log2FoldChange <0)),
                   mainEffect_HSD_UP = rownames(subset(mainEffect_HSD_female,padj<0.1 & log2FoldChange >0)),
                   mainEffect_HSD_DOWN = rownames(subset(mainEffect_HSD_female,padj<0.1 & log2FoldChange <0))
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)


# get the shared UP and the shared down to ORA + unique HSD
sharedUpregulated <- intersect(rownames(subset(mainEffect_HSD_female,padj<0.1 & log2FoldChange >0)),
                               rownames(subset(mainEffect_HSD_female,padj<0.1 & log2FoldChange <0)))
length(sharedUpregulated)
HSD_unique <- setdiff(rownames(subset(Day50_effect_onDiet_HSD,padj<0.1 & log2FoldChange >0)),
                      c(rownames(subset(Day10_effect_onDiet_HSD,padj<0.1 & log2FoldChange >0)),
                        rownames(subset(Day50_effect_onDiet_HSD,padj<0.1 & log2FoldChange <0)),
                        rownames(subset(Day10_effect_onDiet_HSD,padj<0.1 & log2FoldChange <0))))

length(HSD_unique)
sharedUpregulated <- HSD_unique
library(org.Dm.eg.db)
# take union of males & female identified genes
universe_entrez_all<- clusterProfiler::bitr(rownames(ddsSE_split),
                                            fromType = "FLYBASE",
                                            toType = "ENTREZID",
                                            OrgDb = org.Dm.eg.db)$ENTREZID

geneSetChoice_tranlsated <- clusterProfiler::bitr(sharedUpregulated,
                                                  fromType="FLYBASE",
                                                  toType="ENTREZID",
                                                  OrgDb=org.Dm.eg.db)$ENTREZID

length(geneSetChoice_tranlsated)
library(msigdbr)
ORA_cluster_results_interaction <- doOra(
  geneSetChoice_tranlsated, # ENSEBML
  type=c("GO","KEGG","HALLMARK"),
  levelGOTerms=8,
  universe_entrez_all,
  GoFilter=T,
  filename = "uniqHSD_UPregulated_interaction_female" # will be ORA_[filename][type].png
)



ggVennDiagram(list(main_interaction_HFD_50_UP = rownames(subset(main_interaction_HFD_50_female,padj<0.1 & log2FoldChange >0)),
                   main_interaction_HFD_50_DOWN = rownames(subset(main_interaction_HFD_50_female,padj<0.1 & log2FoldChange <0)),
                   main_interaction_HSD_50_UP =rownames(subset(main_interaction_HSD_50_female,padj<0.1 & log2FoldChange >0)),
                   main_interaction_HSD_50_DOWN =rownames(subset(main_interaction_HSD_50_female,padj<0.1 & log2FoldChange <0))
                   
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)


ggVennDiagram(list(main_interaction_HFD_10_UP = rownames(subset(main_interaction_HFD_10_female,padj<0.1 & log2FoldChange >0)),
                   main_interaction_HFD_10_DOWN = rownames(subset(main_interaction_HFD_10_female,padj<0.1 & log2FoldChange <0)),
                   main_interaction_HSD_10_UP =rownames(subset(main_interaction_HSD_10_female,padj<0.1 & log2FoldChange >0)),
                   main_interaction_HSD_10_DOWN =rownames(subset(main_interaction_HSD_10_female,padj<0.1 & log2FoldChange <0))
                   
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)



sharedUpregulated <- intersect(rownames(subset(main_interaction_HFD_50_female,padj<0.1 & log2FoldChange <0)),
                               rownames(subset(main_interaction_HSD_50_female,padj<0.1 & log2FoldChange <0)))
length(sharedUpregulated)
HSD_unique <- setdiff(rownames(subset(main_interaction_HSD_50_male,padj<0.1 & log2FoldChange >0)),
                      rownames(subset(main_interaction_HFD_50_male,padj<0.1 & log2FoldChange >0)))

HSD_unique <- setdiff(HSD_unique,intersect(rownames(subset(main_interaction_HSD_50_male,padj<0.1 & log2FoldChange >0)),
                                           rownames(subset(main_interaction_HFD_50_male,padj<0.1 & log2FoldChange <0))))
length(HSD_unique)
sharedUpregulated <- HSD_unique
library(org.Dm.eg.db)
# take union of males & female identified genes
universe_entrez_all <- clusterProfiler::bitr(rownames(ddsSE_split),
                                             fromType = "FLYBASE",
                                             toType = "ENTREZID",
                                             OrgDb = org.Dm.eg.db)$ENTREZID

geneSetChoice_tranlsated <- clusterProfiler::bitr(sharedUpregulated,
                                                  fromType="FLYBASE",
                                                  toType="ENTREZID",
                                                  OrgDb=org.Dm.eg.db)$ENTREZID

length(geneSetChoice_tranlsated)
library(msigdbr)
ORA_cluster_results_interaction <- doOra(
  geneSetChoice_tranlsated, # ENSEBML
  type=c("GO","KEGG","HALLMARK"),
  levelGOTerms=8,
  universe_entrez_all,
  GoFilter=T,
  filename = "shared_DOWNregulated_totalEffect_female" # will be ORA_[filename][type].png
)


# check if there are any genes having a main diet effect but no interaction effect
HSD_unique <- setdiff(rownames(subset(main_interaction_HSD_50_male,padj<0.1)),
                      rownames(subset(main_interaction_HFD_50_male,padj<0.1 & log2FoldChange >0)))


list(Day10_effect_onDiet_HFD_UP =rownames(subset(Day10_effect_onDiet_HFD,padj<0.1 & log2FoldChange >0)),
     Day10_effect_onDiet_HFD_DOWN =rownames(subset(Day10_effect_onDiet_HFD,padj<0.1 & log2FoldChange <0)),
     Day50_effect_onDiet_HFD_UP = rownames(subset(Day50_effect_onDiet_HFD,padj<0.1 & log2FoldChange >0)),
     Day50_effect_onDiet_HFD_DOWN = rownames(subset(Day50_effect_onDiet_HFD,padj<0.1 & log2FoldChange <0))
)

ggVennDiagram(list(mainEffect_HSD =rownames(subset(mainEffect_HSD_female,padj<0.1 )),
                   Day10_effect_onDiet_HSD=rownames(subset(Day10_effect_onDiet_HSD,padj<0.1)),
                   Day50_effect_onDiet_HSD = rownames(subset(Day50_effect_onDiet_HSD,padj<0.1))
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)


HSD_unique <- setdiff(rownames(subset(mainEffect_HSD_female,padj<0.1 & log2FoldChange >0)),
                      c(rownames(subset(Day10_effect_onDiet_HSD,padj<0.1)),
                        rownames(subset(Day50_effect_onDiet_HSD,padj<0.1))
                      ))

geneSetChoice_tranlsated <- clusterProfiler::bitr(HSD_unique,
                                                  fromType="FLYBASE",
                                                  toType="SYMBOL",
                                                  
                                                  OrgDb=org.Dm.eg.db)



## get total effects ----
# get the total effects and do venn diagram with stuff that also as main Day effect



main_interaction_HSD_50_10_female <- results(DE_SE_female, contrast=list( c("Diet_HSD_vs_CD","DietHSD.Day10","DietHSD.Day50") ))
summary(main_interaction_HSD_50_10_female)

mainEffect_Day_50_10_female <- results(DE_SE_female, contrast=c("Day","50","10"))
mainEffect_Day_50_0_female <- results(DE_SE_female, contrast=c("Day","50","0"))
mainEffect_Day_10_0_female <- results(DE_SE_female, contrast=c("Day","10","0"))

summary(mainEffect_Day_10_0_female)

summary(main_interaction_HSD_10_female)


ggVennDiagram(list(main_interaction_HSD_10_female =rownames(subset(main_interaction_HSD_10_female,padj<0.1 )),
                   main_interaction_HSD_50_female=rownames(subset(main_interaction_HSD_50_female,padj<0.1)),
                   main_interaction_HSD_50_10_female=rownames(subset(main_interaction_HSD_50_10_female,padj<0.1)),
                   mainEffect_HSD_female = rownames(subset(mainEffect_HSD_female,padj<0.1)),
                   mainEffect_Days_all = c(rownames(subset(mainEffect_Day_50_10_female,padj<0.1)),
                                           rownames(subset(mainEffect_Day_50_0_female,padj<0.1)),
                                           rownames(subset(mainEffect_Day_10_0_female,padj<0.1)))
                   ),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)

# get those unique numbers
unique_main_interaction_HSD_10_female<-
  setdiff(rownames(subset(main_interaction_HSD_10_female,padj<0.1 )),
        c(rownames(subset(main_interaction_HSD_50_female,padj<0.1 )),
          rownames(subset(main_interaction_HSD_50_10_female,padj<0.1)),
          rownames(subset(mainEffect_HSD_female,padj<0.1)),
          c(rownames(subset(mainEffect_Day_50_10_female,padj<0.1)),
            rownames(subset(mainEffect_Day_50_0_female,padj<0.1)),
            rownames(subset(mainEffect_Day_10_0_female,padj<0.1)))
        ))

length(unique_main_interaction_HSD_10_female)

unique_main_interaction_HSD_50_female<-
  setdiff(rownames(subset(main_interaction_HSD_50_female,padj<0.1 )),
          c(rownames(subset(main_interaction_HSD_10_female,padj<0.1 )),
            rownames(subset(main_interaction_HSD_50_10_female,padj<0.1)),
            rownames(subset(mainEffect_HSD_female,padj<0.1)),
            c(rownames(subset(mainEffect_Day_50_10_female,padj<0.1)),
              rownames(subset(mainEffect_Day_50_0_female,padj<0.1)),
              rownames(subset(mainEffect_Day_10_0_female,padj<0.1)))
          ))

length(unique_main_interaction_HSD_50_female)

unique_main_interaction_HSD_50_10_female <-
  setdiff(rownames(subset(main_interaction_HSD_50_10_female,padj<0.1 )),
          c(rownames(subset(main_interaction_HSD_10_female,padj<0.1 )),
            rownames(subset(main_interaction_HSD_50_female,padj<0.1)),
            rownames(subset(mainEffect_HSD_female,padj<0.1)),
            c(rownames(subset(mainEffect_Day_50_10_female,padj<0.1)),
              rownames(subset(mainEffect_Day_50_0_female,padj<0.1)),
              rownames(subset(mainEffect_Day_10_0_female,padj<0.1)))
          ))

length(unique_main_interaction_HSD_50_10_female)

# subset to the to evaluate up & down

uniq_main_interaction_HSD_10_female_subset <- main_interaction_HSD_10_female[unique_main_interaction_HSD_10_female,]
uniq_main_interaction_HSD_50_female_subset <- main_interaction_HSD_50_female[unique_main_interaction_HSD_50_female,]
uniq_main_interaction_HSD_50_10_female_subset <- main_interaction_HSD_50_10_female[unique_main_interaction_HSD_50_10_female,]

UP_and_down <- list(uniq_main_interaction_HSD_10_female_subset_UP = rownames(subset(uniq_main_interaction_HSD_10_female_subset,padj<0.1 & log2FoldChange>0 )),
                   uniq_main_interaction_HSD_10_female_subset_DOWN = rownames(subset(uniq_main_interaction_HSD_10_female_subset,padj<0.1 & log2FoldChange<0 )),
                   uniq_main_interaction_HSD_50_female_subset_UP = rownames(subset(uniq_main_interaction_HSD_50_female_subset,padj<0.1 & log2FoldChange>0 )),
                   uniq_main_interaction_HSD_50_female_subset_DOWN = rownames(subset(uniq_main_interaction_HSD_50_female_subset,padj<0.1 & log2FoldChange<0 )),
                   uniq_main_interaction_HSD_50_10_female_subset_UP = rownames(subset(uniq_main_interaction_HSD_50_10_female_subset,padj<0.1 & log2FoldChange>0 )),
                   uniq_main_interaction_HSD_50_10_female_subset_DOWN = rownames(subset(uniq_main_interaction_HSD_50_10_female_subset,padj<0.1 & log2FoldChange<0 )))

mainEffectVis(list(uniq_main_interaction_Day10_female_subset=uniq_main_interaction_HSD_10_female_subset,
                   uniq_main_interaction_main_DayBetween_female_subset=uniq_main_interaction_HSD_50_10_female_subset,
                   uniq_main_interaction_Day50_female_subset=uniq_main_interaction_HSD_50_female_subset),color=annoCol)
                   


library(org.Dm.eg.db)
# take union of males & female identified genes
universe_entrez_all <- clusterProfiler::bitr(rownames(DE_SE_female),
                                             fromType = "FLYBASE",
                                             toType = "ENTREZID",
                                             OrgDb = org.Dm.eg.db)$ENTREZID

geneSetChoice_tranlsated <- clusterProfiler::bitr(UP_and_down$uniq_main_interaction_HSD_50_10_female_subset_DOWN,
                                                  fromType="FLYBASE",
                                                  toType="ENTREZID",
                                                  OrgDb=org.Dm.eg.db)$ENTREZID

length(geneSetChoice_tranlsated)
library(msigdbr)
ORA_cluster_results_interaction <- doOra(
  geneSetChoice_tranlsated, # ENSEBML
  type=c("GO","KEGG","HALLMARK"),
  levelGOTerms=6,
  universe_entrez_all,
  GoFilter=T,
  filename = "uniq_DOWNregulated_totalEffectDay50_10_noDay_female" # will be ORA_[filename][type].png
)
View(ORA_cluster_results_interaction$GO)


## total effects across time ----
### males ----
main_interaction_HSD_50_10_male <- results(DE_SE_male, contrast=list( c("Diet_HSD_vs_CD","DietHSD.Day10","DietHSD.Day50") ))

ggVennDiagram(list(main_interaction_HSD_10_male =rownames(subset(main_interaction_HSD_10_male,padj<0.1 )),
                   main_interaction_HSD_50_male=rownames(subset(main_interaction_HSD_50_male,padj<0.1)),
                   main_interaction_HSD_50_10_female=rownames(subset(main_interaction_HSD_50_10_male,padj<0.1))
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)

main_interaction_HFD_50_10_male <- results(DE_SE_male, contrast=list( c("Diet_HFD_vs_CD","DietHFD.Day10","DietHFD.Day50") ))

ggVennDiagram(list(main_interaction_HFD_10_male =rownames(subset(main_interaction_HFD_10_male,padj<0.1 )),
                   main_interaction_HFD_50_male=rownames(subset(main_interaction_HFD_50_male,padj<0.1)),
                   main_interaction_HFD_50_10_male=rownames(subset(main_interaction_HFD_50_10_male,padj<0.1))
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)


ggVennDiagram(list(main_interaction_HFD_10_male_UP =rownames(subset(main_interaction_HFD_10_male,padj<0.1&log2FoldChange <0)),
                   main_interaction_HFD_50_male_UP=rownames(subset(main_interaction_HFD_50_male,padj<0.1& log2FoldChange <0)),
                   main_interaction_HFD_50_10_male_UP=rownames(subset(main_interaction_HFD_50_10_male,padj<0.1&log2FoldChange <0))
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)


union <- unique(unlist(list(main_interaction_HFD_10_male_UP =rownames(subset(main_interaction_HFD_10_male,padj<0.1&log2FoldChange <0)),
                            main_interaction_HFD_50_male_UP=rownames(subset(main_interaction_HFD_50_male,padj<0.1& log2FoldChange <0)),
                            main_interaction_HFD_50_10_male_UP=rownames(subset(main_interaction_HFD_50_10_male,padj<0.1&log2FoldChange <0)))))


geneSetChoice_tranlsated <- clusterProfiler::bitr(union,
                                                  fromType="FLYBASE",
                                                  toType="ENTREZID",
                                                  OrgDb=org.Dm.eg.db)$ENTREZID

length(geneSetChoice_tranlsated)
library(msigdbr)
ORA_cluster_results_interaction <- doOra(
  geneSetChoice_tranlsated, # ENSEBML
  type=c("GO","KEGG","HALLMARK"),
  levelGOTerms=6,
  universe_entrez_all,
  GoFilter=T,
  filename = "union_DOWNregulated_HFD_totalEffectDay_allDays_male" # will be ORA_[filename][type].png
)
View(ORA_cluster_results_interaction$GO)

ORA_cluster_results_interaction$HALLMARK[grepl("GLYCOLYSI",ORA_cluster_results_interaction$HALLMARK$Description),c("Description","GeneRatio","qvalue","geneNames")]
ORA_cluster_results_interaction$KEGG[grepl("GLYCOLYSI",ORA_cluster_results_interaction$KEGG$Description),c("Description","GeneRatio","qvalue","geneNames")]
ORA_cluster_results_interaction$HALLMARK[grepl("INFLAMMATORY",ORA_cluster_results_interaction$HALLMARK$Description),c("Description","GeneRatio","qvalue","geneNames")]
ORA_cluster_results_interaction$KEGG[grepl("OXIDATIVE",ORA_cluster_results_interaction$KEGG$Description),c("Description","GeneRatio","qvalue","geneNames")]
ORA_cluster_results_interaction$GO[grepl("immun",ORA_cluster_results_interaction$GO$Description),c("Description","GeneRatio","qvalue","geneNames")]


# Mothers ----
## here we load kallisto data ----
se_fly_mothers <- dds_kallisto_mothers
se_fly_mothers$Gender <- "female"
se_fly_mothers$Condition <- gsub("s","S",se_fly_mothers$Condition)
se_fly_mothers$Condition <- as.factor(se_fly_mothers$Condition)
se_fly_mothers$Merged <- se_fly_mothers$Condition

## Do inital filtering ----

library("DESeq2")

ddsSE_mothers <- DESeqDataSet(se_fly_mothers, design = ~ Condition)
ddsSE_mothers$Diet <- relevel(ddsSE_mothers$Condition, ref = "CD")

design(ddsSE_mothers)
ddsSE_preFilter_mothers <- ddsSE_mothers
ddsSE_mothers <- preprocessing(
  ddsSE_preFilter_mothers,
  10,
  protCodingOnly=T,
  removeConstRows=T,
  filterPerSample=T)

de_seq_mothers <- DESeq(ddsSE_mothers) 

de_seq_all_vst_mothers <- vst(de_seq_mothers, blind=T)

library(ggplot2)
PCA_Data <- doPCA(
  de_seq_all_vst_mothers,
  xPC = "PC1",
  yPC = "PC2",
  colorTheme = colorTheme,
  shapeVar = NULL, # one of colnames in colData(dds)
  colorVar = "Condition" # one of colnames in colData(dds)
)

PCA_Data


ggplot(PCA_Data$data, 
       aes(x = PC1, 
           y = PC2, 
           color = Conc.pg_Per_uL.)) +
  scale_color_gradient(low = "blue", high = "red") +
  geom_point(size = 3) +
  coord_fixed()+
  theme_classic()+
  theme(text=element_text(size = 21),aspect.ratio = 1)

ggplot(PCA_Data$data, 
       aes(x = PC1, 
           y = PC2, 
           color = Rin)) +
  scale_color_gradient(low = "blue", high = "red") +
  geom_point(size = 3) +
  coord_fixed()+
  theme_classic()+
  theme(text=element_text(size = 21),aspect.ratio = 1)


summary(results(de_seq_mothers,contrast = list("Condition_HFD_vs_CD")))

summary(results(de_seq_mothers,contrast = list("Condition_HSD_vs_CD")))
