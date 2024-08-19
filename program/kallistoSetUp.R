# Setting Up with kallisto input

## Offspring
allSampleData <- read.csv("../data/SampleData_prepped.csv",row.names = 1)
allSampleData$ID <- allSampleData$SampleID
allSampleData$Merged <- as.factor(allSampleData$Condition)

annot=rtracklayer::import("/Volumes/My_Book/Seep_Lea/current/fly_sugar/data/index/annotations.gtf")
annot=as.data.frame(annot)
c("GENEID", "SYMBOL", "GENETYPE")
annot <- annot[,c("gene_id","transcript_id","gene_name","transcript_biotype")]
tx_anno_df <- annot
colnames(tx_anno_df) <- c("GENEID", "TXNAME", "SYMBOL", "GENETYPE")

#remove all NAs
tx_anno_df<-tx_anno_df[!is.na(tx_anno_df$TXNAME),]

dds_kallisto <- doKallisto2dds(kallisto_path = "/Volumes/My_Book/Seep_Lea/current/fly_sugar/data/output/",
               sampleTable = allSampleData, # one column needs to be called ID!
               tx_anno_path = NULL,
               tx_anno_df = tx_anno_df)

rowData(dds_kallisto)$GENETYPE <- "none"

for(i in rownames(rowData(dds_kallisto))){
  tmp <- unique(subset(tx_anno_df,GENEID==i)$GENETYPE)
  rowData(dds_kallisto)[i,"GENETYPE"] <- tmp
}

saveRDS(dds_kallisto,"../data/DDS_kallisto_all.rds")


## Mothers

motherSampleData <- read.csv("../data/SampleTableMothers.csv",row.names = 1)
motherSampleData$ID <- motherSampleData$SampleID
motherSampleData$Condition <- motherSampleData$Diet
motherSampleData$Merged <- as.factor(motherSampleData$Condition)


dds_kallisto_mothers <- doKallisto2dds(kallisto_path = "/Volumes/My_Book/Seep_Lea/current/fly_sugar/data/output/",
                               sampleTable = motherSampleData, # one column needs to be called ID!
                               tx_anno_path = NULL,
                               tx_anno_df = tx_anno_df)

rowData(dds_kallisto_mothers)$GENETYPE <- "none"

for(i in rownames(rowData(dds_kallisto_mothers))){
  tmp <- unique(subset(tx_anno_df,GENEID==i)$GENETYPE)
  rowData(dds_kallisto_mothers)[i,"GENETYPE"] <- tmp
}

saveRDS(dds_kallisto_mothers,"../data/DDS_kallisto_mothers.rds")
