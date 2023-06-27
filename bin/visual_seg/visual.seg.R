suppressWarnings(suppressMessages(library(tibble)))
suppressWarnings(suppressMessages(library(plyr)))
plot_segment_label <- function(filepath=NULL, filename=NULL, data,assembly='hg38',no_label=F,ylim=NULL,...){
  colnames(data)[1] <- 'sampleID'
  data <- data %>% add_column(arm=rep('.',dim(data)[1]),.after = 2)
  if (no_label){
    col = '#d69f7e'
  } else{
    col <- c('red','#f8961e','#90be6d','#8ecae6','#014f86')
    col_idx <- revalue(as.character(data$label),replace = c('0'=3,'+1'=2,'+2'=1,'-1'=4,'-2'=5),warn_missing=F)
    col <- col[as.numeric(col_idx)]
  }

  
  if (!is.null(filepath) & !is.null(filename)){
    pdf(paste0(filepath,'/',filename), width = 10, height = 5)
  } else if (!is.null(filepath) & is.null(filename)){
    pdf(paste0(filepath,data[,1][1],".pdf"), width = 10, height = 5)
  }
  plotGenome(segments = data,connect=FALSE,seg.col=col,assembly = assembly,ylim=ylim,...)
  if (!is.null(filepath)){
    dev.off()
  }
}
