# This main script can be sourced to run the entire analysis pipeline
# Input:
# - DDS_kallisto_mothers.rds
# - DDS_kallisto_allrds

# Output:
# - PCA-plots
# - Venn-diagrams based on DEG
# - Enrichment analysis
# - heatmaps for genes of interest

# Set up ----
library(DESeq2)
library(svglite)
library(ggplot2)
library(ggVennDiagram)
library(org.Dm.eg.db)
library(msigdbr)
library(UpSetR)
library(readxl)

setwd("program")

source("utils/doPCA.R")
source("utils/doOra.R")
source("utils/preprocessing_dds.R")
source("utils/vis_routine.R")

# changed to due to color mix up
colorTheme <- c("#797979","#4581BE","#F52F92","#e28743","#76b5c5","#dde5b6","#7cb518","#31572c")
names(colorTheme) <- c("CD","HSD","HFD","Male","Female","0","10","50")

# Mothers ----
## Preparation data ----
dds_kallisto_mothers <- readRDS("../data/DDS_kallisto_mothers.rds")
se_fly_mothers <- dds_kallisto_mothers

write.table(file ="../data/countMatrix_mothers.csv", assay(dds_kallisto_mothers),sep = ",")
se_fly_mothers$Gender <- "female"
se_fly_mothers$Condition <- gsub("s","S",se_fly_mothers$Condition)
se_fly_mothers$Condition <- as.factor(se_fly_mothers$Condition)
se_fly_mothers$Merged <- se_fly_mothers$Condition

# Do inital filtering 

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

## DESeq2 ----
de_seq_mothers <- DESeq(ddsSE_mothers) 

## PCA ----
de_seq_all_vst_mothers <- vst(de_seq_mothers, blind=T)

PCA_Data <- doPCA(
  de_seq_all_vst_mothers,
  xPC = "PC1",
  yPC = "PC2",
  colorTheme = colorTheme,
  shapeVar = NULL, # one of colnames in colData(dds)
  colorVar = "Condition" # one of colnames in colData(dds)
)+
  theme_bw(base_size = 8)+
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        plot.margin = margin(1, 1, 1, 1),   
        aspect.ratio = 1)

svglite::svglite("../results/PCA_mothers.svg",width=2.7,height=2.7)
plot(PCA_Data)
dev.off()

summary(results(de_seq_mothers,contrast = list("Condition_HFD_vs_CD")))
summary(results(de_seq_mothers,contrast = list("Condition_HSD_vs_CD")))

# do Venn diagrams for up and down reg

mainEffect_HFD_mothers <- results(de_seq_mothers, contrast=list( c("Condition_HFD_vs_CD") ))
summary(mainEffect_HFD_mothers)
UP_HFD <- rownames(mainEffect_HFD_mothers)[mainEffect_HFD_mothers$padj<0.1 & mainEffect_HFD_mothers$log2FoldChange>0]
UP_HFD <- UP_HFD[!is.na(UP_HFD)]
DOWN_HFD <- rownames(mainEffect_HFD_mothers)[mainEffect_HFD_mothers$padj<0.1 & mainEffect_HFD_mothers$log2FoldChange<0]
DOWN_HFD <- DOWN_HFD[!is.na(DOWN_HFD)]

mainEffect_HSD_mothers <- results(de_seq_mothers, contrast=list( c("Condition_HSD_vs_CD") ))
summary(mainEffect_HSD_mothers)
UP_HSD <- rownames(mainEffect_HSD_mothers)[mainEffect_HSD_mothers$padj<0.1 & mainEffect_HSD_mothers$log2FoldChange>0]
UP_HSD <- UP_HSD[!is.na(UP_HSD)]
DOWN_HSD <- rownames(mainEffect_HSD_mothers)[mainEffect_HSD_mothers$padj<0.1 & mainEffect_HSD_mothers$log2FoldChange<0]
DOWN_HSD <- DOWN_HSD[!is.na(DOWN_HSD)]


svglite::svglite("../results/VennDiagram_mothers_UP.svg",width = 10, height = 10)
ggVennDiagram(list(UP_HFD = UP_HFD,
                   UP_HSD = UP_HSD
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)
dev.off()

intersect(UP_HFD,UP_HSD)

svglite::svglite("../results/VennDiagram_mothers_DOWN.svg",width = 10, height = 10)
ggVennDiagram(list(DOWN_HFD = DOWN_HFD,
                   DOWN_HSD = DOWN_HSD
),
label = "count")+ 
  scale_x_continuous(expand = expansion(mult = .2))+ 
  scale_fill_distiller(palette = "Reds", direction = 1)
dev.off()

allSets <- list(UP_HFD = UP_HFD,
                UP_HSD = UP_HSD,
                DOWN_HFD = DOWN_HFD,
                DOWN_HSD = DOWN_HSD)

universe_entrez_all_mothers <- clusterProfiler::bitr(rownames(de_seq_mothers),
                                                     fromType = "FLYBASE",
                                                     toType = "ENTREZID",
                                                     OrgDb = org.Dm.eg.db)$ENTREZID


for(i in names(allSets)){
  geneSetChoice_tranlsated <- clusterProfiler::bitr(allSets[[i]],
                                                    fromType="FLYBASE",
                                                    toType="ENTREZID",
                                                    OrgDb=org.Dm.eg.db)$ENTREZID
  
  
  length(geneSetChoice_tranlsated)
  ORA_cluster_results_interaction <- doOra(
    geneSetChoice_tranlsated, # ENSEBML
    type=c("GO","KEGG","HALLMARK"),
    levelGOTerms=6,
    universe_entrez_all_mothers,
    GoFilter=T,
    filename = paste0(i,"_mothers") # will be ORA_[filename][type].png
  )
  
  
  if(!("KEGG" %in% names(ORA_cluster_results_interaction))){
    ORA_cluster_results_interaction$KEGG <- c("Nothing enriched")
  }
  if(!("HALLMARK" %in% names(ORA_cluster_results_interaction))){
    ORA_cluster_results_interaction$HALLMARK <- c("Nothing enriched")
  }
  if(!("GO" %in% names(ORA_cluster_results_interaction))){
    ORA_cluster_results_interaction$GO <- c("Nothing enriched")
  }
  
  openxlsx::write.xlsx(ORA_cluster_results_interaction,file = paste0("../results/mothers_",i,".xlsx"))
  
}

# Progeny ----
## Preparation data ----
dds_kallisto_progeny <- readRDS("../data/DDS_kallisto_all.rds")

se_fly_all <- dds_kallisto
write.table(file ="../data/countMatrix_progeny.csv", assay(dds_kallisto),sep = ",")

se_fly_all$Gender <- as.factor(se_fly_all$Gender)
se_fly_all$Condition <- as.factor(se_fly_all$Condition)
se_fly_all$Day <- factor(se_fly_all$Day,levels=c("0","10","50"))
se_fly_all$Diet <- as.factor(se_fly_all$Diet)

# remove Outliers identifed via PCA
# If you want to see the before PCA plot uncomment the following line
se_fly_all <- se_fly_all[,!colnames(se_fly_all)%in%c("M8","F12")]

# create sex separated summarised Experiment objects
se_fly_male <- se_fly_all[,se_fly_all$Gender=="Male"]
se_fly_female <- se_fly_all[,se_fly_all$Gender=="Female"]

# Do inital filtering 
ddsSE <- DESeqDataSet(se_fly_all, design = ~ Diet + Day + Diet:Day)
ddsSE$Day <- relevel(ddsSE$Day, ref = "0")
ddsSE$Diet <- relevel(ddsSE$Diet, ref = "CD")

design(ddsSE)
ddsSE_preFilter <- ddsSE
ddsSE <- preprocessing(ddsSE_preFilter,10,
                        protCodingOnly=T,
                        removeConstRows=T,
                        filterPerSample=T)

## DESeq2 & PCA both ----

de_seq_all <- DESeq(ddsSE) 
resultsNames(de_seq_all)
summary(de_seq_all)
res <- results(de_seq_all)

# PCA all 
de_seq_all_vst <- vst(de_seq_all,blind=T)

de_seq_all_vst$Gender = factor(de_seq_all_vst$Gender ,
                               levels = c("Female","Male"),
                               ordered = T)
#reorder columns for plotting purposes when outlier present
idx_2 <- which(grepl("M8|F12",colnames(de_seq_all_vst)))
idx_1 <- which(!grepl("M8|F12",colnames(de_seq_all_vst)))
de_seq_all_vst <- de_seq_all_vst[,c(idx_1,idx_2)]

doPCA(
  de_seq_all_vst,
  colorTheme = colorTheme,
  shapeVar = "Diet", # one of colnames in colData(dds)
  colorVar = "Gender" # one of colnames in colData(dds)
)

## PCA male only ----
ddsSE_male <- DESeqDataSet(se_fly_male, design = ~ Diet + Day + Diet:Day)
ddsSE_male$Day <- relevel(ddsSE_male$Day, ref = "0")
ddsSE_male$Diet <- relevel(ddsSE_male$Diet, ref = "CD")

design(ddsSE_male)
ddsSE_preFilter_male <- ddsSE_male
ddsSE_male <- preprocessing(ddsSE_preFilter_male,10,
                       protCodingOnly=T,
                       removeConstRows=T,
                       filterPerSample=T)

de_seq_all_vst_male <- vst(ddsSE_male,blind=T)
PCA_Data <- doPCA(
  de_seq_all_vst_male,
  #xPC="PC2",
  #yPC="PC3",
  colorTheme = colorTheme,
  shapeVar = "Diet", # one of colnames in colData(dds)
  colorVar = "Day" # one of colnames in colData(dds)
)

#remove age signature for plotting within PCA only
removedBatch_Male <- limma::removeBatchEffect(assay(de_seq_all_vst_male), 
                                              covariates = model.matrix(~de_seq_all_vst_male$Day + de_seq_all_vst_male$Diet:de_seq_all_vst_male$Day),
                                              design = model.matrix(~de_seq_all_vst_male$Diet))


removedBatch_Male_SE <- de_seq_all_vst_male
assay(removedBatch_Male_SE) <- removedBatch_Male

PCA_Data<-doPCA(
  removedBatch_Male_SE,
  colorTheme = colorTheme,
  shapeVar = "Day", # one of colnames in colData(dds)
  colorVar = "Diet" # one of colnames in colData(dds)
)

## PCA female only ----
ddsSE_female <- DESeqDataSet(se_fly_female, design = ~ Diet + Day + Diet:Day)
ddsSE_female$Day <- relevel(ddsSE_female$Day, ref = "0")
ddsSE_female$Diet <- relevel(ddsSE_female$Diet, ref = "CD")

design(ddsSE_female)
ddsSE_preFilter_male <- ddsSE_female
ddsSE_female <- preprocessing(ddsSE_preFilter_male,10,
                            protCodingOnly=T,
                            removeConstRows=T,
                            filterPerSample=T)


de_seq_all_vst_female <- vst(ddsSE_female,blind=T)
PCA_Data<-doPCA(
  de_seq_all_vst_female,
  xPC="PC1",
  yPC="PC2",
  colorTheme = colorTheme,
  shapeVar = "Diet",
  colorVar = "Day" 
)

#remove age signature for plotting within PCA only
removedBatch_Female <- limma::removeBatchEffect(assay(de_seq_all_vst_female), 
                                              covariates = model.matrix(~de_seq_all_vst_female$Day + de_seq_all_vst_female$Diet:de_seq_all_vst_female$Day),
                                              design = model.matrix(~de_seq_all_vst_female$Diet))


removedBatch_Female_SE <- de_seq_all_vst_female
assay(removedBatch_Female_SE) <- removedBatch_Female

PCA_Data<-doPCA(
  removedBatch_Female_SE,
  colorTheme = colorTheme,
  shapeVar = "Day", 
  colorVar = "Diet"
)

# DE-Overlap Bar graphs----
colors2use <- annoCol
names(colors2use$Day) <- paste0("Day",names(colors2use$Day))
DE_SE_female <- DESeq(ddsSE_female)
DE_SE_male <- DESeq(ddsSE_male)

## female ----
# Day
res_list <- list(main_Day10_female = results(DE_SE_female, contrast=list(c("Day_10_vs_0"))),
                 main_Day50_female = results(DE_SE_female, contrast=list(c("Day_50_vs_0"))))

mainEffectVis(res_list,color=colors2use)

# Diet
res_list <- list(main_HFD_female = results(DE_SE_female, contrast=list(c("Diet_HFD_vs_CD"))),
                 main_HSD_female = results(DE_SE_female, contrast=list(c("Diet_HSD_vs_CD"))))

mainEffectVis(res_list,color=colors2use)

# total
res_list <- list(totelEffect_Day10_HFD = results(DE_SE_female, contrast=list(c("Diet_HFD_vs_CD","DietHFD.Day10"))),
                 totelEffect_Day50_HFD = results(DE_SE_female, contrast=list(c("Diet_HFD_vs_CD","DietHFD.Day50"))))

mainEffectVis(res_list,color=colors2use)

res_list <- list(totelEffect_Day10_HSD = results(DE_SE_female, contrast=list(c("Diet_HSD_vs_CD","DietHSD.Day10"))),
                 totelEffect_Day50_HSD = results(DE_SE_female, contrast=list(c("Diet_HSD_vs_CD","DietHSD.Day50"))))

mainEffectVis(res_list,color=colors2use)

## male ----
# Day
res_list <- list(main_Day10_male = results(DE_SE_male, contrast=list(c("Day_10_vs_0"))),
                 main_Day50_male = results(DE_SE_male, contrast=list(c("Day_50_vs_0"))))

mainEffectVis(res_list,color=colors2use)
# Diet
res_list <- list(main_HFD_male = results(DE_SE_male, contrast=list(c("Diet_HFD_vs_CD"))),
                 main_HSD_male = results(DE_SE_male, contrast=list(c("Diet_HSD_vs_CD"))))

mainEffectVis(res_list,color=colors2use)

# total
res_list <- list(totelEffect_Day10_HFD = results(DE_SE_male, contrast=list(c("Diet_HFD_vs_CD","DietHFD.Day10"))),
                 totelEffect_Day50_HFD = results(DE_SE_male, contrast=list(c("Diet_HFD_vs_CD","DietHFD.Day50"))))

mainEffectVis(res_list,color=colors2use)

res_list <- list(totelEffect_Day10_HSD = results(DE_SE_male, contrast=list(c("Diet_HSD_vs_CD","DietHSD.Day10"))),
                 totelEffect_Day50_HSD = results(DE_SE_male, contrast=list(c("Diet_HSD_vs_CD","DietHSD.Day50"))))

mainEffectVis(res_list,color = colors2use)


# Set comparison ----
resultsListAll_withinGender <- list()
namesOnly_withinGender <- list()

counter = 0
for(i in c("HSD","HFD")){
  for(j in c("Female","Male")){
    if(j=="Female"){
      de_seq_tmp <- DE_SE_female
    }else{
      de_seq_tmp <- DE_SE_male
    }
    
    counter = counter + 1
    
    dif0 = paste0("pAdjOnly_Diet_",i,"_CD_Day0_",j)
    dif10 = paste0("pAdjOnly_Diet_",i,"_CD_Day10_",j)
    dif50 = paste0("pAdjOnly_Diet_",i,"_CD_Day50_",j)
    
    comp0 = paste0("Diet_",i,"_vs_CD")
    comp10 = paste0("Diet",i,".Day10")
    comp50 = paste0("Diet",i,".Day50")
    
    # Diet
    #Day0
    res_tmp <- results(de_seq_tmp, contrast = list(c(comp0)))
    resultsListAll_withinGender[[dif0]] <- res_tmp
    namesOnly_withinGender[[paste0(dif0,"_UP")]] <- rownames(subset(
      results(de_seq_tmp,contrast = list(c(comp0))),
      padj < 0.05 & log2FoldChange > 0)) 
    namesOnly_withinGender[[paste0(dif0,"_DOWN")]] <- rownames(subset(
      results(de_seq_tmp,contrast = list(c(comp0))),
      padj < 0.05 & log2FoldChange < 0))
    
    #Day10
    res_tmp <- results(de_seq_tmp, contrast = list(comp0,comp10))
    resultsListAll_withinGender[[dif10]] <- res_tmp
    namesOnly_withinGender[[paste0(dif10,"_UP")]] <- rownames(subset(
      results(de_seq_tmp,contrast = list(c(comp0,comp10))),
      padj < 0.05 & log2FoldChange > 0)) 
    namesOnly_withinGender[[paste0(dif10,"_DOWN")]] <- rownames(subset(
      results(de_seq_tmp,contrast = list(c(comp0,comp10))),
      padj < 0.05 & log2FoldChange < 0))
    
    #Day50
    res_tmp <- results(de_seq_tmp, contrast = list(c(comp0,comp50)))
    resultsListAll_withinGender[[dif50]] <- res_tmp
    namesOnly_withinGender[[paste0(dif50,"_UP")]] <- rownames(subset(
      results(de_seq_tmp,contrast = list(c(comp0,comp50))),
      padj < 0.05 & log2FoldChange > 0)) 
    namesOnly_withinGender[[paste0(dif50,"_DOWN")]] <- rownames(subset(
      results(de_seq_tmp,contrast = list(c(comp0,comp50))),
      padj < 0.05 & log2FoldChange < 0))
    
    print(paste0("Done: ",counter,"/9"))
  }
}


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


library(msigdbr)
library(org.Dm.eg.db)
# take union of males & female identified genes
universe_entrez_all <- clusterProfiler::bitr(rownames(ddsSE_split),
                                            fromType = "FLYBASE",
                                            toType = "ENTREZID",
                                            OrgDb = org.Dm.eg.db)$ENTREZID

for(diet in c("HFD","HSD")){
  for(sex in c("Male","Female")){
    data_day0 <- openxlsx::read.xlsx(paste0("../results/DEG_lists_withinDay_DietDiff_",sex,".xlsx"),sheet = paste0(diet,"_CD_Day_0_",sex,"_DEG"))
    data_day10 <- openxlsx::read.xlsx(paste0("../results/DEG_lists_withinDay_DietDiff_",sex,".xlsx"),sheet = paste0(diet,"_CD_Day_10_",sex,"_DEG"))
    data_day50 <- openxlsx::read.xlsx(paste0("../results/DEG_lists_withinDay_DietDiff_",sex,".xlsx"),sheet = paste0(diet,"_CD_Day_50_",sex,"_DEG"))
    
    gg_day_venn_UP <- ggVennDiagram(list(data_day0_UP = subset(data_day0,padj < 0.1 & log2FoldChange > 0)$GeneName,
                                      data_day10_UP = subset(data_day10,padj < 0.1 & log2FoldChange > 0)$GeneName,
                                      data_day50_UP = subset(data_day50,padj < 0.1 & log2FoldChange > 0)$GeneName
    ),
    label = "count")+ 
      scale_x_continuous(expand = expansion(mult = .2))+ 
      scale_fill_distiller(palette = "Reds", direction = 1)
    
    gg_day_venn_DOWN <- ggVennDiagram(list(data_day0_DOWN =subset(data_day0,padj<0.1 & log2FoldChange <0)$GeneName,
                                         data_day10_DOWN =subset(data_day10,padj<0.1 & log2FoldChange <0)$GeneName,
                                         data_day50_DOWN = subset(data_day50,padj<0.1 & log2FoldChange <0)$GeneName
    ),
    label = "count")+ 
      scale_x_continuous(expand = expansion(mult = .2))+ 
      scale_fill_distiller(palette = "Reds", direction = 1)
    
    ggsave(paste0("../results/ggVenn_designGroup_",sex,"_",diet,"_DOWN.svg"),plot=gg_day_venn_DOWN)
    ggsave(paste0("../results/ggVenn_designGroup_",sex,"_",diet,"_UP.svg"),plot=gg_day_venn_UP)

    # do enrichment for all Diet (not the interactions)
    unique_UP_diet <- subset(data_day50,padj<0.1 & log2FoldChange >0)$GeneID
    unique_DOWN_diet <- subset(data_day50,padj<0.1 & log2FoldChange <0)$GeneID
    
    sets <- list(unique_UP = unique_UP_diet,
                 unique_DOWN_HSD = unique_DOWN_diet
                 )
    
    for(i in names(sets)){
      geneSetChoice_tranlsated <- clusterProfiler::bitr(sets[[i]],
                                                        fromType="FLYBASE",
                                                        toType="ENTREZID",
                                                        OrgDb=org.Dm.eg.db)$ENTREZID
      
      filename <- paste0("../results/",sex,"_",diet,"_",i)
      ORA_cluster_results_interaction <- doOra(
        geneSetChoice_tranlsated, # ENSEBML
        type=c("GO","KEGG","HALLMARK"),
        levelGOTerms=8,
        universe_entrez_all,
        GoFilter=T,
        filename = filename # will be ORA_[filename][type].png
      )
      
      openxlsx::write.xlsx(ORA_cluster_results_interaction,file = paste0(filename,".xlsx"))
      
    }
  }
}



# Enrichment ----
library(msigdbr)
library(org.Dm.eg.db)
# take union of males & female identified genes
universe_entrez_all <- clusterProfiler::bitr(rownames(ddsSE_split),
                                             fromType = "FLYBASE",
                                             toType = "ENTREZID",
                                             OrgDb = org.Dm.eg.db)$ENTREZID

for(diet in c("HFD","HSD")){
  for(sex in c("Male","Female")){
    data_day0 <- openxlsx::read.xlsx(paste0("../results/DEG_lists_withinDay_DietDiff_",sex,".xlsx"),sheet = paste0(diet,"_CD_Day_0_",sex,"_DEG"))
    data_day10 <- openxlsx::read.xlsx(paste0("../results/DEG_lists_withinDay_DietDiff_",sex,".xlsx"),sheet = paste0(diet,"_CD_Day_10_",sex,"_DEG"))
    data_day50 <- openxlsx::read.xlsx(paste0("../results/DEG_lists_withinDay_DietDiff_",sex,".xlsx"),sheet = paste0(diet,"_CD_Day_50_",sex,"_DEG"))
    
    gg_day_venn_UP <- ggVennDiagram(list(data_day0_UP = subset(data_day0,padj < 0.1 & log2FoldChange > 0)$GeneName,
                                         data_day10_UP = subset(data_day10,padj < 0.1 & log2FoldChange > 0)$GeneName,
                                         data_day50_UP = subset(data_day50,padj < 0.1 & log2FoldChange > 0)$GeneName
    ),
    label = "count")+ 
      scale_x_continuous(expand = expansion(mult = .2))+ 
      scale_fill_distiller(palette = "Reds", direction = 1)
    
    gg_day_venn_DOWN <- ggVennDiagram(list(data_day0_DOWN =subset(data_day0,padj<0.1 & log2FoldChange <0)$GeneName,
                                           data_day10_DOWN =subset(data_day10,padj<0.1 & log2FoldChange <0)$GeneName,
                                           data_day50_DOWN = subset(data_day50,padj<0.1 & log2FoldChange <0)$GeneName
    ),
    label = "count")+ 
      scale_x_continuous(expand = expansion(mult = .2))+ 
      scale_fill_distiller(palette = "Reds", direction = 1)
    
    ggsave(paste0("../results/ggVenn_designGroup_",sex,"_",diet,"_DOWN.svg"),plot=gg_day_venn_DOWN)
    ggsave(paste0("../results/ggVenn_designGroup_",sex,"_",diet,"_UP.svg"),plot=gg_day_venn_UP)
    
    # do enrichment for all Diet (not the interactions)
    unique_UP_diet <- subset(data_day50,padj<0.1 & log2FoldChange >0)$GeneID
    unique_DOWN_diet <- subset(data_day50,padj<0.1 & log2FoldChange <0)$GeneID
    
    sets <- list(unique_UP = unique_UP_diet,
                 unique_DOWN_HSD = unique_DOWN_diet
    )
    
    for(i in names(sets)){
      geneSetChoice_tranlsated <- clusterProfiler::bitr(sets[[i]],
                                                        fromType="FLYBASE",
                                                        toType="ENTREZID",
                                                        OrgDb=org.Dm.eg.db)$ENTREZID
      
      filename <- paste0("../results/",sex,"_",diet,"_",i)
      ORA_cluster_results_interaction <- doOra(
        geneSetChoice_tranlsated, # ENSEBML
        type=c("GO","KEGG","HALLMARK"),
        levelGOTerms=8,
        universe_entrez_all,
        GoFilter=T,
        filename = filename # will be ORA_[filename][type].png
      )
      
      openxlsx::write.xlsx(ORA_cluster_results_interaction,file = paste0(filename,".xlsx"))
      
    }
  }
}


## Note that from provided all enrichments, selected terms were specified and further shown in the following section
# was placed in results section
# vis selected terms 
sheetNames <- openxlsx::getSheetNames("../results/Selected_terms_males_EM.xlsx")
selectedTerms_males_list <- lapply(sheetNames, function(x) openxlsx::read.xlsx("../results/Selected_terms_males_EM.xlsx",sheet = x))
names(selectedTerms_males_list) <- sheetNames

sheetNames <- openxlsx::getSheetNames("../results/Selected_terms_females_EM.xlsx")
selectedTerms_female_list <- lapply(sheetNames, function(x) openxlsx::read.xlsx("../results/Selected_terms_females_EM.xlsx",sheet = x))
names(selectedTerms_female_list) <- sheetNames

for(i in names(selectedTerms_males_list)){
  print(nrow(selectedTerms_males_list[[i]]))
  selectedTerms_males_list[[i]]$category <- i
  selectedTerms_males_list[[i]]$sex <- "male"
}

for(i in names(selectedTerms_female_list)){
  print(nrow(selectedTerms_female_list[[i]]))
  selectedTerms_female_list[[i]]$category <- i
  selectedTerms_female_list[[i]]$sex <- "female"
}

#append both lists to big df
selectedTerms_m <- rlist::list.rbind(selectedTerms_males_list)
selectedTerms_f <- rlist::list.rbind(selectedTerms_female_list)

selectedTerms <- rbind(selectedTerms_m,selectedTerms_f)

# get gene ratio right
selectedTerms$GeneRatio_perc <- unlist(lapply(strsplit(selectedTerms$GeneRatio,"\\/"),function(x) as.numeric(x[1])/as.numeric(x[2])*100))


# plot the selected terms
# plot males:

selectedTerms$Day <- as.factor(selectedTerms$Day)
selectedTerms <- selectedTerms[order(selectedTerms$Day,decreasing = F),]

selectedTerms$Description <- gsub("KEGG","",selectedTerms$Description)
selectedTerms$Description <- gsub("_"," ",selectedTerms$Description)
selectedTerms$Description <- tolower(selectedTerms$Description)
#selectedTerms$Count <- cut(selectedTerms$Count, breaks = c(0,10,20,max(selectedTerms$Count)), labels = c("<10","10-20",">20"))

selectedTerms$origin <- "GO"
selectedTerms[grepl("KEGG",selectedTerms$ID),"origin"]<-"KEGG"

# Define colors for each category
category_colors <- c("KEGG" = "blue", "GO" = "green")
selectedTerms_m <- subset(selectedTerms, sex == "male")
selectedTerms_m$Description <- factor(selectedTerms_m$Description, levels = unique(selectedTerms_m$Description))

svglite::svglite("../results/selectedTerms_male.svg",width = 11,height = 10)
ggplot(selectedTerms_m, aes(x = Day, y = Description, color=p.adjust, size=Count)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  #scale_size_binned_area(breaks = c(0,20),max_size = 7) +
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  facet_wrap(~category, scales = "free_y") +
  theme_bw(base_size = 14) +
  labs(title = "Male - Selected terms", y = "Term")

dev.off()


selectedTerms_f <- subset(selectedTerms, sex=="female")
# need to sort describtion within each category
selectedTerms_f <- selectedTerms_f[order(selectedTerms_f$category,selectedTerms_f$Day,decreasing = F),]
selectedTerms_f$fixOrder <- 1:nrow(selectedTerms_f)
selectedTerms_f$Description_new <- paste0(selectedTerms_f$Description,"__",selectedTerms_f$fixOrder)
selectedTerms_f$Description_new <- factor(selectedTerms_f$Description_new, 
                                          levels = selectedTerms_f$Description_new)


library(stringr)
svglite::svglite("../results/selectedTerms_female.svg",width = 11,height = 10)
ggplot(selectedTerms_f, aes(x=Day, y=Description_new, color=p.adjust, size=Count)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  #scale_size_binned_area(breaks = c(0,20),max_size = 7) +
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  facet_wrap(~category, scales = "free_y") +
  scale_y_discrete(labels = function(y) str_replace_all(y, "__.*$", ""))+
  theme_bw(base_size = 14) +
  labs(title = "Female - Selected terms", y = "Term")
dev.off()



# Visualise expression from genes of interest ----
path <- "../data/geneLists_Seyhmus/Selected_GENES_for_heatmap_drosophila.xlsx"
Selected_list_provided <- list()
for(i in readxl::excel_sheets(path)){
  Selected_list_provided[[  gsub(" ","_",i)]] <- 
    read_excel(path,sheet=i)
}

# plot the given genes in a pheatmap to check clustering
dfToPlot <- as.matrix(assay(de_seq_all_vst)) 
rownames(dfToPlot) <- rowData(de_seq_all_vst)$SYMBOL
translatedSubset <- unlist(Selected_list_provided$Is_r_genes,use.names = F)

#check which are still present
keep <- translatedSubset[translatedSubset %in% rownames(dfToPlot)]
colData(de_seq_all_vst)$Group <- as.factor(colData(de_seq_all_vst)$Group)

# separate by gender the heatmaps
samples2keep <- colData(de_seq_all_vst)$Gender %in% "Male" 
svglite::svglite("Is_r_genes_Male.svg",width=10,height=10)
pheatmap::pheatmap(dfToPlot[keep[!apply(dfToPlot[keep,samples2keep],1,sd)==0],samples2keep],
                   color = myColor,
                   scale = "row",
                   fontsize_number = 15,
                   cluster_rows = T,
                   cluster_cols = T,
                   show_rownames = T,
                   show_colnames = F,
                   cutree_rows = 3,
                   cutree_cols = 4,
                   annotation_colors = annoCol_tmp,
                   annotation_col  = col_annp_df)
dev.off()

samples2keep <- colData(de_seq_all_vst)$Gender %in% "Female" 
svglite::svglite("Is_r_genes_Female.svg",width=10,height=10)
pheatmap::pheatmap(dfToPlot[keep[!apply(dfToPlot[keep,samples2keep],1,sd)==0],samples2keep],
                   color = myColor,
                   scale = "row",
                   fontsize_number = 15,
                   cluster_rows = T,
                   cluster_cols = T,
                   show_rownames = T,
                   show_colnames = F,
                   cutree_rows = 3,
                   cutree_cols = 4,
                   annotation_colors = annoCol_tmp,
                   annotation_col  = col_annp_df)
dev.off()

