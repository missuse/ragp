#' Generate domain diagram of domains and AG regions
#'
#' 
#'@param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#'@return a ggplot2 visualisation of domain structure
#' 
#'@details 
#'
#'@references 
#'
#'@seealso \code{\link[ragp]{maab} \link[ragp]{scan_ag}
#'
#'@examples
#'
#'@export
#'
plot_domains <- function(sequences,
                         annotations){
  
  domcount <- NULL
  for(i in names(sequences)){
    domcount[i] <- sum(annotations[,1]==i)
  }
  annotations <- data.frame(y = rep(1:(length(sequences)),domcount),
                            annotations)
  
  agregions  <- scan_ag(sequence = sequences,
                        id = names(sequences),
                        dim = 3,
                        div = 6,
                        type = "extended", simplify = FALSE)$locations
  agregions2 <- matrix(unlist(lapply(agregions,t)), ncol = 2, byrow = TRUE)
  agregcount <- unlist(lapply(agregions, nrow))
  
  agregions3 <- data.frame(y           = rep(1:(length(sequences)),agregcount),
                           id          = rep(names(sequences),agregcount),
                           align_start = agregions2[,1],
                           align_end   = agregions2[,2])
  
  # backbones
  y <- rep(1:(length(sequences)),2)
  x <- c(rep(0,length(sequences)),
         nchar(sequences))
  
  d.backbone=data.frame(x,y)
  
  # AG regions
  y1 <- agregions3$y - 0.2
  y2 <- agregions3$y + 0.2
  x1 <- agregions3$align_start
  x2 <- agregions3$align_end
  domain  <- rep("AGP",nrow(agregions3))
  
  d.agregions=data.frame(x1,x2,y1,y2,domain)
  
  # domains
  y1 <- annotations$y - 0.2
  y2 <- annotations$y + 0.2
  x1 <- annotations$align_start
  x2 <- annotations$align_end
  domain  <- annotations$name
  
  d.domains=data.frame(x1,x2,y1,y2,domain)
  d.domains <- d.domains[!is.na(d.domains[,1]),]
  
  
  # plot
  d <- rbind(d.backbone,
             d.agregions,
             d.domains)
  p <- ggplot() + 
    scale_x_continuous(name="Length") + 
    scale_y_continuous(name="Protein", breaks = 1:length(sequences), labels=names(sequences)) +
    geom_line(data=d.backbone,  mapping=aes(x=x, y=y, group=y)) +
    geom_rect(data=d.domains, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=domain), color="black", alpha=0.7) +
    geom_rect(data=d.agregions, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="White", color="black") +
    geom_text(data=d.domains, aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=domain), size=3) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank())
  
  p
}   