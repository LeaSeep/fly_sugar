# Does the Over-Representation Analysis for the KEGG, GO and Hallmark sets 
# (murine sets)
# Input:
#  - gene set to test
#  - universe to test against
# Output:
#  - Enrichment results as a list


doOra <- function(
    geneSetChoice_tranlsated,
    type=NULL,
    levelGOTerms=6,
    universe_entrez=NULL,
    filename,
    qvalue_cutoff=0.05,
    GoFilter=F){
  
  result = list()


  if("GO" %in% type){
    EnrichmentRes_GO <- clusterProfiler::enrichGO(gene = geneSetChoice_tranlsated,
                                                  OrgDb="org.Dm.eg.db",
                                                  pvalueCutoff = 0.1,
                                                  keyType = 'ENTREZID',
                                                  ont = "BP",
                                                  universe = universe_entrez,
                                                  qvalueCutoff = qvalue_cutoff,
                                                  readable=T)

    tryCatch(
      {
        # The following could be used to combine at a certain level
        EnrichmentRes_GO_filter <- EnrichmentRes_GO
        if(GoFilter){
          EnrichmentRes_GO_filter = clusterProfiler::gofilter(EnrichmentRes_GO, level = levelGOTerms)
        }
        EnrichmentRes_GO_filter_simple <- clusterProfiler::simplify(EnrichmentRes_GO_filter,
                                                                    cutoff = 0.8,
                                                                    by = "p.adjust",
                                                                    select_fun=min)
        },
      error=function(e){
        print("Simplyfing Go w.r.t to levels did not work, try to lower the number; the original object is returned and further processed")
      }
    )

    if(!exists("EnrichmentRes_GO_filter_simple")){
      EnrichmentRes_GO_filter_simple <- EnrichmentRes_GO
    }

    if(all(EnrichmentRes_GO_filter_simple@result$p.adjust>=qvalue_cutoff)){
      print("GO nothing enriched - check result object")
    }else{
      png(paste0(filename,"_ORA","_GO.png"))
      print(clusterProfiler::dotplot(EnrichmentRes_GO_filter_simple))
      dev.off()
      svg(paste0(filename,"_ORA","_GO.svg"))
      print(clusterProfiler::dotplot(EnrichmentRes_GO_filter_simple))
      dev.off()
    }
    
    # add a column to calc cluster gene ration
    EnrichmentRes_GO_filter_simple@result$GeneRatio_perc <- 
      unlist(lapply(strsplit(EnrichmentRes_GO_filter_simple@result$GeneRatio,"/"),function(x){
        tmp <- as.numeric(x)
        new <- round(tmp[1]/tmp[2],4)*100
        paste0(new,"%")
    }))
    
    # add a column to get genes in term , done in GO
    EnrichmentRes_GO_filter_simple@result$geneNames <- EnrichmentRes_GO_filter_simple@result$geneID

    result[["GO"]] <- EnrichmentRes_GO_filter_simple@result

  }

  if("KEGG" %in% type){
    KEGGset <- msigdbr(
      species = "Drosophila melanogaster",
      category = "C2",
      subcategory = "KEGG"
    ) %>% dplyr::select(gs_name, entrez_gene)
    EnrichmentRes_Kegg <- clusterProfiler::enricher(
      gene = geneSetChoice_tranlsated,
      pvalueCutoff = qvalue_cutoff, # otherwise taken this for adj value to determine if enriched
      qvalueCutoff = qvalue_cutoff,
      pAdjustMethod = "fdr",
      #universe = universe_entrez,
      minGSSize = 10,
      TERM2GENE = KEGGset
    )
    
    if(!is.null(EnrichmentRes_Kegg)){
      if(all(EnrichmentRes_Kegg@result$p.adjust >= qvalue_cutoff)){
        print("Kegg nothing enriched - check result object")
      }else{
        png(paste0(filename,"_ORA","_KEGG.png"))
        print(clusterProfiler::dotplot(EnrichmentRes_Kegg))
        dev.off()
        svg(paste0(filename,"_ORA","_KEGG.svg"))
        print(clusterProfiler::dotplot(EnrichmentRes_Kegg))
        dev.off()
      }
      
      # add a column to calc cluster gene ration
      EnrichmentRes_GO_filter_simple@result$GeneRatio_perc <- 
        unlist(lapply(strsplit(EnrichmentRes_GO_filter_simple@result$GeneRatio,"/"),function(x){
          tmp <- as.numeric(x)
          new <- round(tmp[1]/tmp[2],4)*100
          paste0(new,"%")
        }))
      
      
      # add a column to get genes in term
      EnrichmentRes_Kegg@result$geneNames <- 
       unlist(lapply(strsplit(EnrichmentRes_Kegg@result$geneID,"/"),function(x){
          tmp <- clusterProfiler::bitr(x,
                                       fromType="ENTREZID",
                                       toType="SYMBOL",
                                       OrgDb="org.Dm.eg.db")$SYMBOL
          tmp <- paste0(as.character(tmp),collapse = "/")
        }))
      


      
      result[["KEGG"]] <- EnrichmentRes_Kegg@result
    }else{
      print("consider changing 'minGSSize' within function")
    }

  }
  
  if("HALLMARK" %in% type){
    Hallmarkset <- msigdbr(
      species ="Drosophila melanogaster",
      category = "H",
    ) %>% dplyr::select(gs_name, entrez_gene)
    EnrichmentRes_Hallmarks <- clusterProfiler::enricher(
      gene = geneSetChoice_tranlsated,
      pvalueCutoff = qvalue_cutoff, # otherwise taken this for adj value to determine if enriched
      qvalueCutoff = qvalue_cutoff,
      minGSSize = 10,
      pAdjustMethod = "fdr",
      universe = universe_entrez,
      TERM2GENE = Hallmarkset
    )

    if( !is.null(EnrichmentRes_Hallmarks) ){    
      if( all(EnrichmentRes_Hallmarks@result$p.adjust >= qvalue_cutoff) ){
        print("Hallmark nothing enriched - check result object")
      }else{
        png(paste0(filename,"_ORA","_HALLMARK.png"))
        print(clusterProfiler::dotplot(EnrichmentRes_Hallmarks))
        dev.off()
        svg(paste0(filename,"_ORA","_HALLMARK.svg"))
        print(clusterProfiler::dotplot(EnrichmentRes_Hallmarks))
        dev.off()
      }
      
      # add a column to calc cluster gene ration
      EnrichmentRes_Hallmarks@result$GeneRatio_perc <- 
        unlist(lapply(strsplit(EnrichmentRes_Hallmarks@result$GeneRatio,"/"),function(x){
          tmp <- as.numeric(x)
          new <- round(tmp[1]/tmp[2],4)*100
          paste0(new,"%")
        }))
      
      
      # add a column to get genes in term
      EnrichmentRes_Hallmarks@result$geneNames <- 
        unlist(lapply(strsplit(EnrichmentRes_Hallmarks@result$geneID,"/"),function(x){
          tmp <- clusterProfiler::bitr(x,
                                       fromType="ENTREZID",
                                       toType="SYMBOL",
                                       OrgDb="org.Dm.eg.db")$SYMBOL
          tmp <- paste0(as.character(tmp),collapse = "/")
        }))
    
      result[["HALLMARK"]] <- EnrichmentRes_Hallmarks@result
    }else{
      print("consider changing 'minGSSize' within function")
    }
  }
  
  return(result)
}
