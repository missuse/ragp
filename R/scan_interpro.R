#' Find the distribution of inter-proline distances within each sequence
#'
#' 
#'@param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#'@param truncate An integer defining the inter-proline distance above which all distances are combined into a single bin (default = 20)
#'
#'@return matrix of inter-proline distances for each sequence (with overflow bin for long inter-proline distances)
#' 
#'@details
#'@seealso \code{\link[ragp]{scan_ag} \link[ragp]{scan_nglc} \link[ragp]{plot_interpro_heatmap}}
#'
#'@examples
#' interP  <- scan_interpro(sequences)
#' heatmap <- plot_interpro_heatmap(interP)
#' heatmap$newick
#' 
#'@export

scan_interpro <- function(sequences, truncate = 20){
  
  interP.dist <- lapply(stringi::stri_split(sequences,regex = "P"),
                        stringi::stri_length)
  
  interP.tidy <- lapply (interP.dist,
                         hist,
                         breaks = -1:max(unlist(interP.dist)),plot=FALSE)
  interP.tidy <- lapply(interP.tidy, `[[`, 2)
  interP.tidy <- t(array(unlist(interP.tidy),
                         dim = c(max(unlist(interP.dist))+1,
                                 length(sequences))))
  
  colnames(interP.tidy)<- 0:max(unlist(interP.dist))
  rownames(interP.tidy)<- names(sequences)
  
  if(is.numeric(truncate)){
    interP <- cbind(interP.tidy[,1:(truncate+1)],
                    rowSums(interP.tidy[,(truncate+2):max(unlist(interP.dist))]))
    colnames(interP)[truncate+2]<-paste0(">",truncate)
  }else{
    interP <- interP.tidy
  }
  
  interP
}

#' Visualise the distribution of inter-proline distances within each sequence
#'
#' 
#'@param interP A matrix of inter-proline distances for each sequence
#'@param normalise whether to normalise interproline distances for the heatmap
#'
#'@return heatmap of inter-proline distances for each sequences
#' 
#'@details
#'@seealso \code{\link[ragp]{scan_interpro} \link[ragp]{plot_domains}}
#'
#'@examples
#' interP  <- scan_interpro(sequences)
#' heatmap <- plot_interpro_heatmap(interP)
#' heatmap$newick
#' 
#'@export

plot_interpro_heatmap <- function(interP, normalise = TRUE){
  if(normalise){
    interP <- interP/rowSums(interP)
  }
  heatmap <- heatmap(interP,
                     Colv = NA,
                     keep.dendro = TRUE)
  newick <- ape::write.tree(
    ape::as.phylo(
      as.hclust(
        heatmap$Rowv
      )))
  list(heatmap = heatmap,
       newick = newick)
}
