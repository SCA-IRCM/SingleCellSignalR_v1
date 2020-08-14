#' @title data_prepare_pdx
#' @description Prepares the data for further analysis
#'
#'
#' @details `file_human` and `file_mouse` is the path to the file containing the read or UMI count
#' matrix the user wants to analyze.
#' @details
#' `lower` and `upper` are used to remove the genes whose average counts are outliers.
#' The values of these arguments are fractions of the total number of genes and
#' hence must be between 0 and 1. Namely, if `lower = 0.05`, then the function
#' removes the 5% less expressed genes and if `upper = 0.05`, then the function
#' removes the 5% most expressed genes.
#' @details
#' If `normalize` is FALSE, then the function skips the 99th percentile
#' normalization and the log transformation.
#' @details
#' If `write` is TRUE, then the function writes two text files. One for the
#' normalized and gene thresholded read counts table and another one for the
#' genes that passed the lower and upper threshold. Note that the length of the
#' genes vector written in the *genes.txt* file is equal to the number of rows
#' of the table of read counts written in the *data.txt* file.
#'
#' @param lower a number in [0,1], low quantile threshold
#' @param upper a number in [0,1], high quantile threshold
#' @param normalize a logical, if TRUE, then computes 99th percentile normalization
#' @param write a logical
#' @param verbose a logical
#' @param file.human a string for the scRNAseq data file
#' @param file.mouse a string for the scRNAseq data file
#'
#' @return The function returns a data frame of filtered and/or normalized data
#' with genes as row names.
#'
#' @export
#'
#' @import data.table
#' @importFrom stats quantile
#' @import graphics
#' @importFrom stats sd quantile
#'
#'
#' @examples
#' file.human <- system.file("scRNAseq_dataset.txt",package = "SingleCellSignalR")
#' file.mouse <- system.file("mouse_scRNAseq_dataset.txt",package = "SingleCellSignalR")
#' data <- data_prepare_pdx(file.human = file.human, file.mouse = file.mouse)

data_prepare_pdx = function(file.human, file.mouse, lower=0, upper=0, normalize=TRUE, write=FALSE, verbose=TRUE){
  if (dir.exists("data")==FALSE & write==TRUE){
    dir.create("data")
  }

  data.human = fread(file.human,data.table = FALSE)
  data.mouse = fread(file.mouse,data.table = FALSE)

  genes.human = as.character(data.human[,1])
  genes.mouse = as.character(data.mouse[,1])

  data.human = data.human[,-1]
  data.mouse = data.mouse[,-1]

  o.human = order(rowSums(data.human))
  o.mouse = order(rowSums(data.mouse))

  data.human = data.human[o.human,]
  genes.human = genes.human[o.human]
  data.mouse = data.mouse[o.mouse,]
  genes.mouse = genes.mouse[o.mouse]

  data.human = data.human[!duplicated(genes.human),]
  genes.human = genes.human[!duplicated(genes.human)]
  data.mouse = data.mouse[!duplicated(genes.mouse),]
  genes.mouse = genes.mouse[!duplicated(genes.mouse)]

  rownames(data.human) = genes.human
  rownames(data.mouse) = genes.mouse

  if (sum(genes.human %in% c(mm2Hs$`Gene name`))==0){
    stop("Please convert human dataset gene ID's (Ensembl or NCBI ID) to official HUGO
        gene symbols")
  }
  if (sum(genes.mouse %in% c(mm2Hs$`Mouse gene name`))==0){
    stop("Please convert mouse dataset gene ID's (Ensembl or NCBI ID) to official HUGO
        gene symbols")
  }

  data.human = data.human[rowSums(data.human)>0,]
  data.human = data.frame(data.human[,apply(data.human,2,function(x) quantile(x,0.99))>0])
  if (normalize==TRUE){
    cat("log-Normalization data human",fill=TRUE)
    q = apply(data.human,2,quantile,0.99)
    data.human = log(1+sweep(data.human,2,q/median(q),"/"))
  }
  data.human = data.human[rowSums(data.human)>0,]
  data.human = data.human[rowSums(data.human)<quantile(rowSums(data.human),1-upper) & rowSums(data.human)>quantile(rowSums(data.human),lower),]

  data.mouse = data.mouse[rowSums(data.mouse)>0,]
  data.mouse = data.frame(data.mouse[,apply(data.mouse,2,function(x) quantile(x,0.99))>0])
  if (normalize==TRUE){
    cat("log-Normalization data mouse",fill=TRUE)
    q = apply(data.mouse,2,quantile,0.99)
    data.mouse = log(1+sweep(data.mouse,2,q/median(q),"/"))
  }
  data.mouse = data.mouse[rowSums(data.mouse)>0,]
  data.mouse = data.mouse[rowSums(data.mouse)<quantile(rowSums(data.mouse),1-upper) & rowSums(data.mouse)>quantile(rowSums(data.mouse),lower),]

  app.hum = matrix(0,ncol=ncol(data.human),nrow = nrow(data.mouse))
  colnames(app.hum) = colnames(data.human)
  app.mus = matrix(0,ncol=ncol(data.mouse),nrow = nrow(data.human))
  colnames(app.mus) = colnames(data.mouse)

  mix.data = cbind(rbind(data.human,app.hum),rbind(app.mus,data.mouse))

  colnames(mix.data) = c(paste0(colnames(data.human),"_human"),paste0(colnames(data.mouse),"_mouse"))
  rownames(mix.data) = c(paste0(rownames(data.human),"_human"),paste0(rownames(data.mouse),"_mouse"))


  if (write==TRUE){
    fwrite(data.frame(mix.data),"./data/mixed_data.txt",sep="\t")
    fwrite(data.frame(rownames(mix.data)),"./data/mixed_genes.txt",sep="\t")
  }

  res = mix.data
  return(res)
}
