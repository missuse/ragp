#' Generate domain diagram of domains and AG regions
#'
#' 
#'@param sequences   A list with elements of class "SeqFastaAA"
#'@param annotations A data frame of annotations produced by \code{\link[ragp]{get_hmm}}
#'@return a ggplot2 visualisation of domain structure
#' 
#'@details 
#'
#'@references 
#'
#'@seealso \code{\link[ragp]{maab} \link[ragp]{scan_ag}}
#'
#'@examples
#'annotations <- get_hmm(sequence = sequences,
#'                       id = names(sequences),
#'                       verbose = FALSE,
#'                       sleep = 0.1)
#'gpis <- ragp::get_big_pi(sequence = sequences,
#'                         id = names(sequences),
#'                         verbose = FALSE,
#'                         sleep = 0.1)   
#'plot_domains(sequences   = sequences,
#'             annotations = annotations,
#'             gpis        = gpis) 
#'      
#'@export
#'
plot_domains <- function(sequences,
                         annotations,
                         gpis       = NULL,
                         signal     = NULL,
                         nglc       = TRUE,
                         labels     = FALSE,
                         dom_limits = TRUE,
                         width      = 0.8){
  
  # Domains ---------
  domcount <- NULL
  for(i in names(sequences)){
    domcount[i] <- sum(annotations[,1]==i, na.rm = 1)
  }
  annotations <- data.frame(y = rep(1:(length(sequences)),domcount),
                            annotations)
  annotations <- annotations[!is.na(annotations$name),]
  
  # Domains data to plot
  y1 <- annotations$y - width/2
  y2 <- annotations$y + width/2
  x1 <- annotations$align_start
  x2 <- annotations$align_end
  domain  <- annotations$name
  
  d.domains <- data.frame(x1,x2,y1,y2,domain)
  
  # Domains (full-length) -------
  dommax <- NULL
  for(i in unique(annotations$name)){
    dommax[i] <- max(annotations$model_end[annotations$name==i],na.rm = 1)
  }
  annotations$full_seq_length <- nchar(sequences[annotations$id])
  annotations$max_start <- annotations$align_start - annotations$model_start
  annotations$max_end   <- annotations$align_end   + (dommax[annotations$name] - annotations$model_end)
  annotations$max_start[annotations$max_start <=0] <- 0
  annotations$max_end[annotations$max_end >= annotations$full_seq_length] <-
             annotations$full_seq_length[annotations$max_end >= annotations$full_seq_length]

  # Full-length domains data to plot
  y1 <- annotations$y - width/2
  y2 <- annotations$y + width/2
  x1 <- annotations$max_start
  x2 <- annotations$max_end
  domain  <- annotations$name
  
  d.domains.max <- data.frame(x1,x2,y1,y2,domain)

  # AGP regions  -----
  agregions  <- ragp::scan_ag(sequence = sequences,
                              id       = names(sequences),
                              dim      = 3,
                              div      = 6,
                              type     = "extended",
                              simplify = FALSE) $locations
  
  agregions2 <- matrix(unlist(lapply(agregions,t)), ncol = 2, byrow = TRUE)
  agregcount <- unlist(lapply(agregions, nrow))
  agregions3 <- data.frame(y           = rep(1:(length(sequences)),agregcount),
                           id          = rep(names(sequences),agregcount),
                           align_start = agregions2[,1],
                           align_end   = agregions2[,2])

  # AGP regions data to plot
  y1 <- agregions3$y - width/2
  y2 <- agregions3$y + width/2
  x1 <- agregions3$align_start
  x2 <- agregions3$align_end
  domain <- rep("AGP",nrow(agregions3))
  
  d.agregions <- data.frame(x1,x2,y1,y2,domain)
  
  ## Hyp predition ---------------------
  hyp <- ragp::predict_hyp(sequence  = sequences,
                           id        = names(sequences),
                           tprob     = 0.3) $prediction
  
  hyp <- hyp[hyp$HYP=="Yes",]
  hyp <- hyp[!is.na(hyp$HYP),]
  hypcount <- NULL
  for(i in names(sequences)){
    hypcount[i] <- sum(hyp$id==i, na.rm = 1)
  }
  hyp <- data.frame(y = rep(1:(length(sequences)),hypcount),
                    hyp)

  # Hyp data to plot
  y1 <- hyp$y - width/2
  y2 <- hyp$y + width/2
  x1 <- hyp$P_pos
  x2 <- hyp$P_pos
  
  y <- c(y1,y2)
  x <- c(x1,x2)
  P <- c(1:length(x1),1:length(x1))
  
  d.hyp=data.frame(x,y,P)
  
  ## GPI anchors ------
  if(!is.null(gpis)){
    gpi <- data.frame(y=1:nrow(gpis), gpis)
    gpi <- gpi[gpi$is.bigpi,]
    
    # GPI data to plot
    y <- gpi$y
    x <- as.numeric(gpi$omega_site)
    d.gpi=data.frame(x,y)
  }
  # N-Glycosylation ---------
  if (!is.null(nglc)){
    nglc <- pred_nglc (sequences)
    nglccount <- NULL
    for(i in names(sequences)){
      nglccount[i] <- sum(nglc$id==i, na.rm = 1)
    }
    nglc <- data.frame(y = rep(1:(length(sequences)),nglccount),
                      nglc)
  
    # N-glc data to plot
    y <- nglc$y
    x <- as.numeric(nglc$align_start)
    d.nglc <- data.frame(x,y)
  }
  
  # Signal peptide -----------------
  
  if(!is.null(signal)){
    signal <- data.frame(y = 1:nrow(signal),
                         signal)
    signal <- signal[!is.na(signal$cut_site),]
    
    # Signal data to plot
    y1 <- signal$y - width/2
    y2 <- signal$y + width/2
    x1 <- rep(0,nrow(signal))
    x2 <- as.numeric(signal$cut_site)
    domain  <- rep("signal",nrow(signal))
    
    d.signal <- data.frame(x1,x2,y1,y2,domain)
  }
  
  # Plot data ---------------------------
  # Backbones data to plot
  y1 <- 1:(length(sequences))
  y2 <- 1:(length(sequences))
  x1 <- rep(0,length(sequences))
  x2 <- nchar(sequences)

  y <- c(y1,y2)
  x <- c(x1,x2)

  d.backbone <- data.frame(x,y)
 
  # Generate plot with ggplot2
  p <- ggplot2::ggplot() 
  p <- p + ggplot2::scale_x_continuous(name="Length", breaks = seq(0,max(nchar(sequences)), by = 100))
  p <- p + ggplot2::scale_y_continuous(name="Protein", breaks = 1:length(sequences), labels=names(sequences))                 
  p <- p + ggplot2::geom_line (data=d.backbone,  mapping=ggplot2::aes(x=x, y=y, group=y))                                               
 
  if(dom_limits){
    p <- p + ggplot2::geom_rect(data=d.domains.max, mapping=ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=domain),
                                color="black", alpha = 0.3)    
  }
  
  p <- p + ggplot2::geom_rect (data=d.agregions, mapping=ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),fill="White", color="black")
  p <- p + ggplot2::geom_rect (data=d.domains,   mapping=ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=domain), color="black")
  p <- p + ggplot2::geom_line (data=d.hyp,       mapping=ggplot2::aes(x=x, y=y, group=P), color="#333333")
  p <- p + ggplot2::geom_rect (data=d.agregions, mapping=ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="black", alpha=0)
  p <- p + ggplot2::geom_point(d.gpi,            mapping=ggplot2::aes(x, y),shape = 21, colour = "black", fill = "yellow", size = width*2.5)
  p <- p + ggplot2::geom_point(d.nglc,           mapping=ggplot2::aes(x, y),shape = 21, colour = "black", fill = "black", size = width*1.5)

  if(!is.null(signal)){
    p <- p + ggplot2::geom_rect (data=d.signal,   mapping=ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="darkblue", color="black")
  }
  
  p <- p + ggplot2::theme_bw()                                                                                                
  p <- p + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                          panel.grid.minor.y = ggplot2::element_blank(),
                          panel.grid.minor.x = ggplot2::element_blank())
  
  if(labels){
    p <- p + ggplot2::geom_text(data=d.domains, ggplot2::aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=domain), size=3)
  }
  
  p
}   
