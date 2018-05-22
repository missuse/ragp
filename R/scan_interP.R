scan_intropro <- function(sequences){
  
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
  
  interP <- cbind(interP.tidy[,1:21],
                        rowSums(interP.tidy[,22:max(unlist(interP.dist))]))
  colnames(interP)[22]<-">20"
  
  interP
}

plot_interP_heatmap <- function(interP, normalise = TRUE){
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

interP  <- scan_intropro(sequences)
heatmap <- plot_interP_heatmap(interP)
heatmap$newick
