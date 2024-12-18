mainEffectVis <- function(res_list,color=annoCol){
  # show bar plot fo each object of the list split in up and down
  plot_df <- data.frame(catgeory = c(),entities=c(),direction = c())
  for(i in 1:length(res_list)){
   tmp <- subset(res_list[[i]],res_list[[i]]$padj < 0.1 & res_list[[i]]$log2FoldChange > 0)
   if(nrow(tmp) > 0){
     plot_df <- rbind(plot_df,data.frame(category = names(res_list)[i],entities=rownames(tmp),direction = "up"))
   }      
   tmp <- subset(res_list[[i]],res_list[[i]]$padj < 0.1 & res_list[[i]]$log2FoldChange < 0)
   if(nrow(tmp) > 0){
     plot_df <- rbind(plot_df,data.frame(category = names(res_list)[i],entities=rownames(tmp),direction = "down"))
   }   
  }
  # add number of count on tob of the bars
  # Calculate counts per direction per category
  library(dplyr)
  count_df <- plot_df %>%
    group_by(direction, category) %>%
    summarise(count = n()) %>%
    ungroup()
  plot_df <- left_join(plot_df, count_df, by = c("direction", "category"))
  plot_df$category <- factor(plot_df$category, levels = unique(plot_df$category))
  plot_df$direction <- factor(plot_df$direction, levels = c("up", "down"))
  
  # Function to get the color
  get_color <- function(category, color) {
    colors <- unlist(color)
    names(colors) <- gsub(".*\\.","",names(colors))
    matched_colors <- sapply(names(colors), function(name) {
      if (grepl(name, category)) {
        return(colors[name])
      }
      return(NA)
    })
    matched_colors <- na.omit(matched_colors)
    if (length(matched_colors) > 0) {
      return(matched_colors[1])
    } else {
      return("#000000") # Default color if no match found
    }
  }
  
  plot_df <- plot_df %>%
    mutate(color = sapply(category, get_color, color))
  
  plot <- ggplot(plot_df,aes(x=direction,fill=color)) + 
    geom_bar()+
    scale_fill_identity()+
    geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.5) + 
    facet_wrap(~category)+
    theme_bw()
  
  print(plot)
  }

