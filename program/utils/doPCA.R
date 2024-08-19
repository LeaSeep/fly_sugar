# Do PCA
# Input:
#   - DESeq2 Object
#   - name of factor for shape
#   - name of factor for color

# Output:
#   - PCA plot
doPCA <- function(dds,
                  colorTheme=NULL,
                  shapeVar,
                  colorVar,
                  xPC = "PC1",
                  yPC = "PC2"
                  ){
  data <- as.data.frame(assay(dds))
  pca <- prcomp(
    as.data.frame(t(data)),
    center = T,
    scale. = T
  )
  
  # how much variance is explained by each PC
  explVar <- pca$sdev^2/sum(pca$sdev^2)
  names(explVar) <- colnames(pca$x)
  # transform variance to percent
  percentVar <- round(100 * explVar, digits = 1)
  
  # Define data for plotting
  pcaData <- data.frame(pca$x,colData(dds))
  if(is.numeric(pcaData[,colorVar])){
    if(is.factor(pcaData[,colorVar])){
      #assume correct ordering of levels
    }else{
      pcaData[,colorVar] = factor(pcaData[,colorVar],
                                  levels = as.character(unique(pcaData[,colorVar])),
                                  ordered = T)
    }
  }else{
    
  }

if(is.null(colorTheme) & !is.null(shapeVar)){
  pca_plot <- ggplot(pcaData, aes(x = pcaData[,xPC],
                                  y = pcaData[,yPC])) +
    geom_point(size =3,aes(shape=pcaData[,shapeVar],fill=pcaData[,colorVar]))+
    scale_shape_manual(values=c(21,24,23,22))+
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    xlab(paste0(names(percentVar[xPC]),": ",percentVar[xPC], "% variance")) +
    ylab(paste0(names(percentVar[yPC]),": ", percentVar[yPC], "% variance")) +
    coord_fixed()+
    theme_classic()+
    theme(text=element_text(size = 21),aspect.ratio = 1)
}else if(!is.null(colorTheme) & is.null(shapeVar)){
  pca_plot <- ggplot(pcaData, aes(x = pcaData[,xPC],
                                  y = pcaData[,yPC])) +
    geom_point(size =3,aes(color = pcaData[,colorVar]))+
    scale_color_manual(values = colorTheme,
                       name = colorVar) +
    xlab(paste0(names(percentVar[xPC]),": ",percentVar[xPC], "% variance")) +
    ylab(paste0(names(percentVar[yPC]),": ", percentVar[yPC], "% variance")) +
    coord_fixed()+
    theme_classic()+
    theme(text=element_text(size = 21),aspect.ratio = 1)
}else{
  pca_plot <- ggplot(pcaData, aes(x = pcaData[,xPC],
                                  y = pcaData[,yPC])) +
    geom_point(size =3,aes(fill = pcaData[,colorVar],shape = pcaData[,shapeVar]))+
    scale_shape_manual(values = c(21,24,23,22)) +
    scale_fill_manual(values = colorTheme,
                       name = colorVar) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    xlab(paste0(names(percentVar[xPC]),": ",percentVar[xPC], "% variance")) +
    ylab(paste0(names(percentVar[yPC]),": ", percentVar[yPC], "% variance")) +
    coord_fixed()+
    theme_classic()+
    theme(text=element_text(size = 21),aspect.ratio = 1)
}

  return(pca_plot)
}
