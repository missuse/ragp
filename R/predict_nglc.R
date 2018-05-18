predict_nglc <- function(sequences){
  
  nglc <- gregexpr('N[^P][ST]', sequences, ignore.case = FALSE)
  countnglc <- unlist(lapply(nglc,length))
  
  nglc2 <- data.frame(id          = rep(names(sequences),countnglc),
                      align_start = unlist(nglc))
  
  nglc2
}
