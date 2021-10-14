#' Protein structure diagram.
#' 
#' Plots a diagram of protein structure based on hmmscan domain annotation and several types of predictions.
#' 
#' @param sequence String representing a protein amino acid sequence.
#' @param id String representing a protein identifier. Will be converted using \code{\link[base]{make.names}}.
#' @param hyp_col Plotting color of predicted hydroxyproline positions. At default set to: '#868686FF'.
#' @param gpi_col Plotting color of the predicted omega site (glycosylphosphatidylinositol attachment). At default set to: '#0073C2FF'.
#' @param nsp_col Plotting color of the N-terminal signal peptide. At default set to: '#CD534CFF'.
#' @param ag_col Plotting color of the AG glycomodul spans. At default set to: '#E5E5E5FF'.
#' @param tm_col Plotting color of the transmembrane regions. At default set to: '#EFC000FF'.
#' @param hyp Boolean, should hydroxyprolines be plotted.
#' @param gpi A string indicating if \code{\link[ragp]{get_big_pi}} (gpi = "bigpi"), \code{\link[ragp]{get_pred_gpi}} (gpi = "predgpi") or \code{\link[ragp]{get_netGPI}} (gpi = "netgpi") should be called when predicting omega sites. To turn off omega site prediction use gpi = "none". At default set to "netgpi". Alternatively the output data frame of the mentioned functions (called with simplify = TRUE) can be supplied.
#' @param nsp A string indicating if \code{\link[ragp]{get_signalp5}} (nsp = "signalp5") or \code{\link[ragp]{get_signalp}} (nsp = "signalp") should be used to obtain N-sp predictions. Alternatively a data frame containing three columns: a character column "id" indicating the protein id as from input, a logical column "is.signalp" and an integer column "sp.length". See \code{\link[ragp]{get_signalp5}} or \code{\link[ragp]{get_signalp}} for details.
#' @param ag Boolean, should the AG glycomodul spans be plotted.
#' @param tm A string indicating if \code{\link[ragp]{get_phobius}} (tm = "phobius") or \code{\link[ragp]{get_tmhmm}} (tm = "tmhmm") should be used to obtain transmembrane region predictions. Alternatively a data frame with two columns:  a character column "id" indicating the protein id as from input and a "prediction" column containing the topology of the transmembrane regions (example "42o81-101i108-126o"). To turn off tm prediction use tm = "none".
#' @param domain Boolean, should the domain predictions obtained using \code{\link[ragp]{get_hmm}}  be plotted. Alternatively the output data frame from \code{\link[ragp]{get_hmm}} can be supplied.
#' @param disorder Boolean, should disordered region predictions obtained using \code{\link[ragp]{get_espritz}} be plotted. Alternatively the output data frame from \code{\link[ragp]{get_espritz}} (called with simplify = TRUE) can be supplied.
#' @param dom_sort One of c("ievalue", "abc", "cba"), defaults to "abc". Domain plotting order. If 'ievalue' domains with the lowest ievalue as determined by hmmscan will be plotted above. If 'abc' or 'cba' the order is determined by domain Names.
#' @param progress Boolean, whether to show the progress bar, at default set to FALSE.
#' @param gpi_size Integer, the size of the gpi symbol. Appropriate values are 1 - 10.
#' @param gpi_shape Integer, the shape of the gpi symbol. Appropriate values are 0 - 25
#' @param hyp_scan Boolean, if ag = TRUE, should \code{\link[ragp]{scan_ag}} be performed on \code{\link[ragp]{predict_hyp}} output thus scanning only arabinogalactan motifs which contain predicted hydroxyprolines.
#' @param ... Appropriate arguments passed to \code{\link[ragp]{get_signalp5}}, \code{\link[ragp]{get_signalp}}, \code{\link[ragp]{get_espritz}}, \code{\link[ragp]{predict_hyp}}, \code{\link[ragp]{get_hmm}} and \code{\link[ragp]{scan_ag}}.
#'
#'
#' @return A ggplot2 plot object
#' 
#' @seealso \code{\link[ragp]{get_signalp5}} \link[ragp]{get_signalp}} \code{\link[ragp]{get_phobius}} \link[ragp]{get_tmhmm} \code{\link[ragp]{get_hmm}} \code{\link[ragp]{get_espritz}} \code{\link[ragp]{predict_hyp}} \code{\link[ragp]{scan_ag}}
#'
#' @examples
#' library(ragp)
#' library(ggplot2)
#' ind <- c(23, 5, 80, 81, 345)
#' pred <- plot_prot(sequence = at_nsp$sequence[ind],
#'                   id = at_nsp$Transcript.id[ind],
#'                   bitscore = 30) #passed to get_hmm
#' pred +
#'   theme(legend.position = "bottom",
#'         legend.direction = "vertical")
#'         
#' #alternatively:      
#' nsp <- get_signalp(data = at_nsp[ind,],
#'                    id = Transcript.id,
#'                    sequence = sequence)
#'                     
#' hmm <- get_hmm(data = at_nsp[ind,],
#'                id = Transcript.id,
#'                sequence = sequence)
#'                
#' gpi <- get_netGPI(data = at_nsp[ind,],
#'                  id = Transcript.id,
#'                  sequence = sequence)                                    
#'
#' tm <- get_phobius(data = at_nsp[ind,],
#'                   id = Transcript.id,
#'                   sequence = sequence)                                                        
#'  
#' disorder <- get_espritz(data = at_nsp[ind,],	
#'                         id = Transcript.id,	
#'                         sequence = sequence)
#'         
#' pred2 <- plot_prot(sequence = at_nsp$sequence[ind],
#'                    id = at_nsp$Transcript.id[ind],
#'                    tm = tm,
#'                    nsp = nsp,
#'                    gpi = gpi,
#'                    domain = hmm,
#'                    disorder = disorder,
#'                    bitscore = 30)
#'                    
#'                    
#' pred2 +
#'   theme(legend.position = "bottom",
#'         legend.direction = "vertical")
#'         
#' #mixing both methods is also a possibility
#'         
#' @export


plot_prot <- function(sequence,
                      id,
                      hyp_col = "#868686FF",
                      gpi_col = "#0073C2FF",
                      nsp_col = "#CD534CFF",
                      ag_col = "#E5E5E5FF",
                      tm_col = "#EFC000FF",
                      hyp = TRUE,
                      gpi = c("bigpi", "predgpi", "netgpi", "none"),
                      nsp = c("signalp", "signalp5", "none"),
                      ag = TRUE,
                      tm = c('phobius', "tmhmm", "none"),
                      domain = TRUE,
                      disorder = FALSE,
                      hyp_scan = if(ag == TRUE && hyp == TRUE) TRUE else FALSE,
                      dom_sort = c("ievalue", "abc", "cba"),
                      progress = FALSE,
                      gpi_size = 4,
                      gpi_shape = 18,
                      ...) {
  
  if (missing(sequence)){
    stop("protein sequence must be provided to obtain predictions",
         call. = FALSE)
  }
  if (missing(id)){
    stop("protein id must be provided to obtain predictions",
         call. = FALSE)
  }
  if (missing(progress)) {
    progress <- FALSE
  }
  if (length(progress) > 1){
    progress <- FALSE
    warning("progress should be of length 1, setting to default: progress = FALSE",
            call. = FALSE)
  }
  if (!is.logical(progress)){
    progress <- as.logical(progress)
    warning("progress is not logical, converting using 'as.logical'",
            call. = FALSE)
  }
  if (is.na(progress)){
    progress <- FALSE
    warning("progress was set to NA, setting to default: progress = FALSE",
            call. = FALSE)
  }
  id <- as.character(id)
  id_label <- id
  id <- make.names(id)
  sequence <- toupper(as.character(sequence))
  if (length(sequence) != length(id)){
    stop("id and sequence vectors are not of same length",
         call. = FALSE)
  }
  sequence <- sub("\\*$", "", sequence)
  aa_regex <- "[^ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv]"
  
  if (any(grepl(aa_regex, sequence))){
    warning(paste("sequences: ",
                  paste(id[grepl(aa_regex, sequence)],
                        collapse = ", "), 
                  " contain symbols not corresponding to amino acids",
                  sep = ""),
            call. = FALSE)
  }
  
  if (missing(dom_sort)){
    dom_sort <- "abc"
  }
  if (!dom_sort %in% c("ievalue", "abc", "cba")) {
    stop("dom_sort should be one of: 'ievalue', 'abc', 'cba",
         call. = FALSE)
  }
  
  areColors <- function(x) {
    sapply(x, function(X) {
      tryCatch(is.matrix(grDevices::col2rgb(X)), 
               error = function(e) FALSE)
    })
  }
  
  if (length(hyp_col) > 1){
    hyp_col <- "#868686FF"
    warning("One color should be provided for hyp_col. Using default: '#868686FF'",
            call. = FALSE)
  }
  
  if (!areColors(hyp_col)){
    hyp_col <- "#868686FF"
    warning("hyp_col provided is not a valid color, default will be used: '#868686FF'",
            call. = FALSE)
  }
  
  if (length(gpi_col) > 1){
    gpi_col <- "#0073C2FF"
    warning("One color should be provided for gpi_col. Using default: '#0073C2FF'",
            call. = FALSE)
  }
  
  if (!areColors(gpi_col)){
    gpi_col <- "#0073C2FF"
    warning("gpi_col provided is not a valid color, default will be used: '#0073C2FF'",
            call. = FALSE)
  }
  
  if (length(nsp_col) > 1){
    nsp_col <- "#CD534CFF"
    warning("One color should be provided for nsp_col. Using default: '#CD534CFF'",
            call. = FALSE)
  }
  
  if (!areColors(nsp_col)){
    nsp_col <- "#CD534CFF"
    warning("nsp_col provided is not a valid color, default will be used: '#CD534CFF'",
            call. = FALSE)
  }
  
  if (length(ag_col) > 1){
    ag_col <- "#E5E5E5FF"
    warning("One color should be provided for ag_col. Using default: '#E5E5E5FF'",
            call. = FALSE)
  }
  
  if (!areColors(ag_col)){
    ag_col <- "#E5E5E5FF"
    warning("ag_col provided is not a valid color, default will be used: '#E5E5E5FF'",
            call. = FALSE)
  }
  
  if (length(tm_col) > 1){
    tm_col <- "#EFC000FF"
    warning("One color should be provided for tm_col. Using default: '#EFC000FF'",
            call. = FALSE)
  }
  
  if (!areColors(tm_col)){
    tm_col <- "#EFC000FF"
    warning("tm_col provided is not a valid color, default will be used: '#EFC000FF'",
            call. = FALSE)
  }
  
  if(missing(tm)){
    tm <-  "phobius"
  }
  if(is.character(tm)){
    if(!tm %in% c("phobius", "tmhmm", "none")){
      tm <-  "phobius"
      warning(paste("tm should be one of",
                    "'phobius'", "'tmhmm'", "'none';",
                    "setting to default: tm = 'phobius'"),
              call. = FALSE)
    }
    if (length(tm) > 1){
      tm <-  "phobius"
      warning("tm should be of length 1, setting to default: tm = 'phobius'",
              call. = FALSE)
    }
  }
  
  if (missing(disorder)){	
    disorder <- FALSE	
  }	
  
  if (is.logical(disorder)){	
    if (length(disorder) > 1){	
      disorder <- FALSE	
      warning("disorder should be of length 1, setting to default: disorder = FALSE",	
              call. = FALSE)	
    }	
  }
  
  if(missing(nsp)){
    nsp <-  "signalp5"
  }
  if(is.character(nsp)){
    if(!nsp %in% c("signalp5", "signalp", "none")){
      nsp <-  "signalp5"
      warning(paste("nsp should be one of",
                    "'signalp5'", "'signalp'", "'none';",
                    "setting to default: nsp = 'signalp5'"),
              call. = FALSE)
    }
    if (length(nsp) > 1){
      nsp <-  "signalp5"
      warning("nsp should be of length 1, setting to default: nsp = 'signalp5'",
              call. = FALSE)
    }
  }
  
  if (missing(ag)){
    ag <- TRUE
  }
  if (length(ag) > 1){
    ag <- TRUE
    warning("ag should be of length 1, setting to default: ag = TRUE",
            call. = FALSE)
  }
  if (!is.logical(ag)){
    ag <- as.logical(ag)
    warning("ag is not logical, converting using 'as.logical'",
            call. = FALSE)
  }
  if (is.na(ag)){
    ag <- TRUE
    warning("ag was set to NA, setting to default: ag = TRUE",
            call. = FALSE)
  }
  
  if (missing(hyp)){
    hyp <- TRUE
  }
  if (length(hyp) > 1){
    hyp <- TRUE
    warning("hyp should be of length 1, setting to default: hyp = TRUE",
            call. = FALSE)
  }
  if (!is.logical(hyp)){
    hyp <- as.logical(hyp)
    warning("hyp is not logical, converting using 'as.logical'",
            call. = FALSE)
  }
  if (is.na(hyp)){
    hyp <- TRUE
    warning("hyp was set to NA, setting to default: hyp = TRUE",
            call. = FALSE)
  }
  
  if (length(hyp_scan) > 1){
    hyp_scan <- if(ag == TRUE && hyp == TRUE) TRUE else FALSE
    warning("hyp_scan should be of length 1, setting to default: hyp_scan = if(ag == TRUE && hyp == TRUE) TRUE else FALSE",
            call. = FALSE)
  }
  if (!is.logical(hyp_scan)){
    hyp_scan <- as.logical(hyp_scan)
    warning("hyp_scan is not logical, converting using 'as.logical'",
            call. = FALSE)
  }
  if (is.na(hyp_scan)){
    hyp_scan <- if(ag == TRUE && hyp == TRUE) TRUE else FALSE
    warning("hyp_scan was set to NA, setting to default: hyp_scan = if(ag == TRUE && hyp == TRUE) TRUE else FALSE",
            call. = FALSE)
  }
  
  if (!hyp && isTRUE(hyp_scan)) {
    stop("hyp must be TRUE if scan_ag is to use hydroxyproline predictions for motif scan")
  }
  
  if (!ag && isTRUE(hyp_scan)) {
    warning("hyp_scan = TRUE and ag = FALSE, arabinogalactan motif scan will not be performed")
  }
  
  if (missing(domain)){
    domain <- TRUE
  }
  
  if (is.logical(domain)){
    if (length(domain) > 1){
      domain <- TRUE
      warning("domain should be of length 1, setting to default: domain = TRUE",
              call. = FALSE)
    }
  }
  
  if(missing(gpi)){
    gpi <- 'netgpi'
  }
  if(is.character(gpi)){
    if(!gpi %in% c("bigpi", "predgpi", "netgpi", "none")){
      gpi <- 'netgpi'
      warning(paste("gpi should be one of",
                    "'bigpi', 'predgpi', 'netgpi', 'none',",
                    "setting to default: gpi = 'netgpi'"),
              call. = FALSE)
    }
    if (length(gpi) > 1){
      gpi <- 'netgpi'
      warning("gpi should be of length 1, setting to default: gpi = 'netgpi'",
              call. = FALSE)
    }
  }
  
  if(length(gpi_size) != 1){
    gpi_size <- 4
    warning("gpi_size should be of length 1, setting to default: gpi_size = 4",
            call. = FALSE)
  }
  if(!gpi_size %in% 1:10) {
    gpi_size <- 4
    warning("gpi_size can have values from 1 - 10, setting to default: gpi_size = 4",
            call. = FALSE)
  }
  
  if(length(gpi_shape) != 1){
    gpi_shape <- 18
    warning("gpi_shape should be of length 1, setting to default: gpi_shape = 18",
            call. = FALSE)
  }
  if(!gpi_shape %in% 0:25) {
    gpi_shape <- 18
    warning("gpi_shape can have values from 0 - 25, setting to default: gpi_shape = 18",
            call. = FALSE)
  }
  
  args_signalp <- names(formals(get_signalp.character))[2:8]
  args_signalp5 <- names(formals(get_signalp5.character))[2]
  args_scanag <- names(formals(scan_ag.default))[4:7]
  args_espritz <- names(formals(get_espritz.default))[4:5]
  args_hmm <- names(formals(get_hmm.default))[9:10]
  dots <- list(...)
  
  dat <- data.frame(sequence = sequence,
                    id = id,
                    nchar = nchar(sequence),
                    stringsAsFactors = FALSE)
  
  dat$id <- factor(dat$id,
                   levels = unique(dat$id))
  
  dat$id_num <- as.numeric(dat$id)
  
  seq_hmm <- NULL
  if (isTRUE(domain)) {
    if(progress){
      message("querying hmmscan")
    }
    seq_hmm <- do.call(ragp::get_hmm,
                       c(list(data = dat,
                              sequence = "sequence",
                              id = "id",
                              progress = progress),
                         dots[names(dots) %in% args_hmm]))
    
    seq_hmm <- seq_hmm[seq_hmm$reported,]
    
    if(nrow(seq_hmm) != 0){
      seq_hmm <- seq_hmm[!is.na(seq_hmm$align_start),]
      
      seq_hmm$id <- factor(seq_hmm$id,
                           levels = unique(dat$id))
      
      seq_hmm$id_num <- as.numeric(seq_hmm$id)
      
      seq_hmm$domain <- with(seq_hmm,
                             paste(desc,
                                   " (",
                                   acc,
                                   ") ",
                                   sep = ""))
      if (dom_sort == "ievalue"){
        seq_hmm <- seq_hmm[with(seq_hmm, order(id_num,
                                               as.numeric(ievalue),
                                               decreasing = c(FALSE, TRUE),
                                               method = "radix")),]
      }
      
      if (dom_sort == "abc"){
        seq_hmm <- seq_hmm[with(seq_hmm, order(id_num, name)),]
      }
      if (dom_sort == "cba"){
        seq_hmm <- seq_hmm[with(seq_hmm, order(id_num,
                                               name,
                                               decreasing = c(FALSE, TRUE),
                                               method = "radix")),]
      }
    } else {
      seq_hmm <- NULL
    }
  } else {
    seq_hmm <- NULL
  }
  
  
  if(is.data.frame(domain)){
    seq_hmm <- domain
    if(any(!c("id",
              "name",
              "acc",
              "desc",
              "align_start",
              "align_end",
              "ievalue",
              "bitscore") %in% colnames(seq_hmm))){
      stop("domain is not the output from get_hmm function")
    }
    seq_hmm$id <- make.names(seq_hmm$id)
    if(!all(seq_hmm$id %in% id)){
      stop("protein ids from domain do not match with id argument")
    }
    seq_hmm <- seq_hmm[seq_hmm$reported,]
    if(any(names(dots) == "ievalue")){
      seq_hmm <- seq_hmm[seq_hmm$ievalue <= dots[names(dots) == "ievalue"],]
    }
    if(any(names(dots) == "bitscore")){
      seq_hmm <- seq_hmm[seq_hmm$bitscore >= dots[names(dots) == "bitscore"],]
    }
    
    if(nrow(seq_hmm) != 0){
      seq_hmm <- seq_hmm[!is.na(seq_hmm$align_start),]
      
      seq_hmm$id <- factor(seq_hmm$id,
                           levels = unique(dat$id))
      
      seq_hmm$id_num <- as.numeric(seq_hmm$id)
      
      seq_hmm$domain <- with(seq_hmm,
                             paste(desc,
                                   " (",
                                   acc,
                                   ") ",
                                   sep = ""))
      if (dom_sort == "ievalue"){
        seq_hmm <- seq_hmm[with(seq_hmm, order(id_num,
                                               as.numeric(ievalue),
                                               decreasing = c(FALSE, TRUE),
                                               method = "radix")),]
      }
      
      if (dom_sort == "abc"){
        seq_hmm <- seq_hmm[with(seq_hmm, order(id_num, name)),]
      }
      if (dom_sort == "cba"){
        seq_hmm <- seq_hmm[with(seq_hmm, order(id_num,
                                               name,
                                               decreasing = c(FALSE, TRUE),
                                               method = "radix")),]
      }
    } else {
      seq_hmm <- NULL
    }
  }
  
  seq_signalp <- NULL
  if (is.character(nsp)){
    if (nsp == "signalp5") {
      if(progress){
        message("querying signalp5")
      }
      seq_signalp <- do.call(ragp::get_signalp5,
                             c(list(data = dat,
                                    sequence = "sequence",
                                    id = "id",
                                    progress = progress),
                               dots[names(dots) %in% args_signalp5]))
      
      seq_signalp <- seq_signalp[seq_signalp$is.signalp,]
      
      if (nrow(seq_signalp) != 0){
        seq_signalp$sp.length <- as.numeric(
          as.character(
            seq_signalp$sp.length
          )
        )
      } else {
        seq_signalp <- NULL
      }
    }
    
    if (nsp == "signalp") {
      if(progress){
        message("querying signalp")
      }
      seq_signalp <- do.call(ragp::get_signalp,
                             c(list(data = dat,
                                    sequence = "sequence",
                                    id = "id",
                                    progress = progress),
                               dots[names(dots) %in% args_signalp]))
      
      seq_signalp <- seq_signalp[seq_signalp$is.signalp,]
      
      if (nrow(seq_signalp) != 0){
        seq_signalp$sp.length <- as.numeric(
          as.character(
            seq_signalp$sp.length
          )
        )
      } else {
        seq_signalp <- NULL
      }
    }
  }
  
  

  if(is.data.frame(nsp)){
    seq_signalp <- nsp
    if(any(!c("id", "is.signalp", "sp.length") %in% colnames(seq_signalp))){
      stop(paste("nsp does not contain columns", "'id'", "'is.signalp'",  "and", "'sp.length'"))
    }
    seq_signalp$id <- make.names(seq_signalp$id)
    if(!all(seq_signalp$id %in% id)){
      stop("protein ids from nsp do not match with id argument")
    }
    if(!is.logical(seq_signalp$is.signalp)){
      stop("is.signalp column is not logical")
    }
    if(!is.numeric(seq_signalp$sp.length)){
      stop("sp.length column is not numeric")
    }
    seq_signalp <- seq_signalp[seq_signalp$is.signalp,]
    if(nrow(seq_signalp) != 0){
      seq_signalp$sp.length <- as.numeric(
        as.character(
          seq_signalp$sp.length
        )
      )
    } else {
      seq_signalp <- NULL
    }
  }
  
  phobius_seq <- NULL
  if (is.character(tm)){
    if (tm == "phobius") {
      if(progress){
        message("querying phobius")
      }
      phobius_seq <- ragp::get_phobius(data = dat,
                                       sequence = sequence,
                                       id = id,
                                       progress = progress)
    }
    
    if (tm == "tmhmm") {
      if(progress){
        message("querying tmhmm")
      }
      phobius_seq <- ragp::get_tmhmm(data = dat,
                                     sequence = sequence,
                                     id = id,
                                     progress = progress)
    }
  }
  
  
  
  if(is.data.frame(tm)){
    phobius_seq <- tm
    if(any(!c("id",
              "prediction") %in% colnames(phobius_seq))){
      stop("tm does not contain the columns 'id' and 'prediction'")
    }
    phobius_seq$id <- make.names(phobius_seq$id)
    if(!all(phobius_seq$id %in% id)){
      stop("protein ids from tm do not match with id argument")
    }
  }
  
  if (!is.null(phobius_seq)){
    phobius_seq <- merge(dat,
                         phobius_seq,
                         by.x = "id",
                         by.y = "id",
                         sort = FALSE)
    
    phobius_seq_pred <- unlist(
      lapply(
        strsplit(
          phobius_seq$prediction,
          "/"),
        function(x){
          if (length(x) == 2) {
            x[2]
          } else {
            x[1]
          }
        }
      )
    )
    phobius_seq <- data.frame(id = phobius_seq$id, 
                              id_num = phobius_seq$id_num,
                              pred = phobius_seq_pred,
                              is.tm = grepl("\\d+-\\d+",
                                            phobius_seq_pred),
                              seq = phobius_seq$sequence,
                              stringsAsFactors = FALSE)
    
    phobius_seq <- phobius_seq[phobius_seq$is.tm,]
    
    if (nrow(phobius_seq) == 0){
      phobius_seq <- NULL
    } else {
      tm <- stringr::str_extract_all(phobius_seq$pred, 
                                     "\\d+-\\d+")
      names(tm) <- phobius_seq$id
      tm <- lapply(tm, cbind)
      id_tm <- rep(phobius_seq$id_num, sapply(tm, nrow))
      tm <- do.call(rbind, tm)
      tm <- data.frame(tm, id_tm)
      tm$tm_start <- as.numeric(stringr::str_extract(tm$tm, 
                                                     "\\d+(?=-)"))
      tm$tm_end <- as.numeric(stringr::str_extract(tm$tm, 
                                                   "(?<=-)\\d+"))
      out <- stringr::str_extract_all(phobius_seq$pred, 
                                      "\\d+o\\d+|\\d+o|o\\d+")
      out <- lapply(out, cbind)
      id_out <- rep(phobius_seq$id_num, sapply(out, nrow))
      seq_out <- nchar(rep(phobius_seq$seq, sapply(out, 
                                                   nrow)))
      out <- do.call(rbind, out)
      out <- data.frame(out, id_out)
      out$out_start <- as.numeric(stringr::str_extract(out$out, 
                                                       "\\d+(?=o)"))
      out$out_start <- ifelse(is.na(out$out_start), 1, 
                              out$out_start)
      out$out_end <- as.numeric(stringr::str_extract(out$out, 
                                                     "(?<=o)\\d+"))
      out$out_end <- ifelse(is.na(out$out_end), seq_out, 
                            out$out_end)
      inside <- stringr::str_extract_all(phobius_seq$pred, 
                                         "\\d+i\\d+|\\d+i|i\\d+")
      inside <- lapply(inside, cbind)
      id_inside <- rep(phobius_seq$id_num, sapply(inside, 
                                                  nrow))
      seq_inside <- nchar(rep(phobius_seq$seq, sapply(inside, 
                                                      nrow)))
      inside <- do.call(rbind, inside)
      inside <- data.frame(inside, id_inside)
      inside$inside_start <- as.numeric(stringr::str_extract(inside$inside, 
                                                             "\\d+(?=i)"))
      inside$inside_start <- ifelse(is.na(inside$inside_start), 
                                    1, inside$inside_start)
      
      inside$inside_end <- as.numeric(stringr::str_extract(inside$inside, 
                                                           "(?<=i)\\d+"))
      inside$inside_end <- ifelse(is.na(inside$inside_end), 
                                  seq_inside,
                                  inside$inside_end)
    }
  }
  
  seq_gpi <- NULL
  if(is.character(gpi)){
    if (gpi == 'bigpi') {
      if(progress){
        message("querying big pi")
      }
      seq_gpi <- ragp::get_big_pi(data = dat,
                                  sequence = "sequence",
                                  id = "id",
                                  progress = progress)
      
      seq_gpi$omega_site <- as.numeric(seq_gpi$omega_site)
      seq_gpi <- seq_gpi[seq_gpi$is.gpi,]
      if (nrow(seq_gpi) == 0) {
        seq_gpi <- NULL
      }
    }
    
    if (gpi == 'predgpi') {
      if(progress){
        message("querying predGPI")
      }
      seq_gpi <- ragp::get_pred_gpi(dat,
                                    sequence = "sequence",
                                    id = "id",
                                    progress = progress)
      seq_gpi <- seq_gpi[seq_gpi$is.gpi,]
      
      if (nrow(seq_gpi) == 0) {
        seq_gpi <- NULL
      }
    }
    
    if (gpi == 'netgpi') {
      if(progress){
        message("querying NetGPI")
      }
      seq_gpi <- ragp::get_netGPI(dat,
                                  sequence = "sequence",
                                  id = "id",
                                  progress = progress)
      seq_gpi <- seq_gpi[seq_gpi$is.gpi,]
      
      if (nrow(seq_gpi) == 0) {
        seq_gpi <- NULL
      }
    }
  }
  
  if(is.data.frame(gpi)){
    seq_gpi <- gpi
    if(any(!c("id",
              "omega_site",
              "is.gpi") %in% colnames(seq_gpi))){
      stop("gpi is not the output from get_big_pi, get_pred_gpi or get_netGPI functions")
    }
    seq_gpi$id <- make.names(seq_gpi$id)
    if(!all(seq_gpi$id %in% id)){
      stop("protein ids from gpi do not match with id argument")
    }
    if(!is.logical(seq_gpi$is.gpi)){
      stop("is.gpi column is not logical")
    }
    seq_gpi <- seq_gpi[seq_gpi$is.gpi,]
    
    if (nrow(seq_gpi) == 0) {
      seq_gpi <- NULL
    }
  } 
  
  
  if (hyp) {
    seq_hyp <- do.call(ragp::predict_hyp,
                       c(list(data = dat,
                              sequence = "sequence",
                              id = "id"),
                         dots[names(dots) == "tprob"],
                         dots[names(dots) == "version"])
    )
    
    if (hyp_scan) hyp_ag_scan <- seq_hyp$sequence
    
    seq_hyp <- seq_hyp$prediction
    
    if (!is.null(seq_hyp)) {
      seq_hyp <- seq_hyp[seq_hyp$HYP == "Yes",]
    }
    
    if (!is.null(seq_hyp) && nrow(seq_hyp) == 0) {
      seq_hyp <- NULL
    }
  } else {
    seq_hyp <- NULL
  }
  
  if (ag && hyp_scan) {
    seq_scan <- do.call(ragp::scan_ag,
                        c(list(data = hyp_ag_scan,
                               sequence = "sequence",
                               id = "id",
                               simplify = FALSE,
                               tidy = TRUE),
                          dots[names(dots) %in% args_scanag]))
    
    seq_scan <- seq_scan[, -5]
    seq_scan <- seq_scan[!is.na(seq_scan$location.start),]
    
    if (nrow(seq_scan) == 0) {
      seq_scan <- NULL
    }
  } else if (ag && !hyp_scan) {
    seq_scan <- do.call(ragp::scan_ag,
                        c(list(data = dat,
                               sequence = "sequence",
                               id = "id",
                               simplify = FALSE,
                               tidy = TRUE),
                          dots[names(dots) %in% args_scanag]))
    
    seq_scan <- seq_scan[, -5]
    seq_scan <- seq_scan[!is.na(seq_scan$location.start),]
    
    if (nrow(seq_scan) == 0) {
      seq_scan <- NULL
    }
  } else {
    seq_scan <- NULL
  }
  
  seq_espritz <- NULL	
  if (isTRUE(disorder)) {	
    if(progress){	
      message("querying espritz")	
    }	
    seq_espritz <- do.call(get_espritz,	
                           c(list(data = dat,	
                                  sequence = "sequence",	
                                  id = "id",	
                                  progress = progress),	
                             dots[names(dots) %in% args_espritz]))	
    
    seq_espritz$id <- factor(seq_espritz$id,	
                             levels = unique(dat$id))	
    
    seq_espritz$id_num <- as.numeric(seq_espritz$id)	
    seq_espritz$start[is.na(seq_espritz$start)] <- 0	
    seq_espritz$end[is.na(seq_espritz$end)] <- 0	
    
    seq_espritz <- merge(seq_espritz,	
                         dat[,c(2,3)],	
                         by.x = "id",	
                         by.y = "id",	
                         all.y = TRUE,	
                         all.x = TRUE,	
                         sort = FALSE) 	
    
    seq_espritz <- split(seq_espritz ,	
                         seq_espritz$id)	
    
    seq_espritz <- lapply(seq_espritz, function(x){	
      start <-  x$start	
      end <- x$end	
      id <- unique(x$id)	
      id_num <- unique(x$id_num)	
      resD <- lapply(seq_along(start), function(i){	
        start[i]:end[i]	
      })	
      groups <- rep(seq_along(resD), times = lengths(resD))	
      data.frame(residue = unlist(resD),	
                 groups = groups,	
                 id = id,	
                 id_num = id_num,	
                 nchar = unique(x$nchar))	
    })	
    
    rle_predO <- lapply(seq_espritz, function(x){	
      n <- unique(x$nchar)	
      n <- 1:n	
      n <- setdiff(n, x$residue)	
      if(length(n) != 0){	
        s <- split(n, cumsum(c(0, diff(n) != 1)))	
        groups <- rep(seq_along(s), times = lengths(s))	
      } else {	
        n <- NA	
        groups <- 1	
      }	
      rle_dat <- data.frame(residue = n,	
                            id = unique(x$id),	
                            id_num = unique(x$id_num),	
                            group = groups)	
      
      rle_dat <- split(rle_dat, rle_dat$group)	
      rle_dat <- lapply(rle_dat, function(x){	
        data.frame(start = min(x$residue),	
                   end = max(x$residue),	
                   id = unique(x$id),	
                   id_num = unique(x$id_num))	
      })	
      rle_dat <- do.call(rbind, rle_dat)	
      rle_dat	
    })	
    
    rle_predO <- do.call(rbind, rle_predO)	
    rownames(rle_predO) <- NULL	
    
    plot_segments <- rle_predO[!is.na(rle_predO$start),]	
    
    seq_espritz <- do.call(rbind, seq_espritz) 	
    
    seq_espritz <- seq_espritz[seq_espritz$residue != 0,]	
    if (nrow(seq_espritz) != 0){	
      rownames(seq_espritz) <- NULL	
      seq_espritz$y <- rep(c(0, -0.1, 0, 0.1),	
                           length.out = nrow(seq_espritz))	
      seq_espritz$y_plot <- seq_espritz$y + seq_espritz$id_num	
      seq_espritz$inter <- interaction(seq_espritz$group,	
                                       seq_espritz$id_num)
    } else {	
      seq_espritz <- NULL	
    }	
  }	
  
  if(is.data.frame(disorder)){	
    seq_espritz <- disorder	
    if(any(!c("start",	
              "end",	
              "id") %in% colnames(seq_espritz))){	
      stop("disorder is not the output from get_espritz function")	
    }	
    seq_espritz$id <- make.names(seq_espritz$id)	
    if(!all(seq_espritz$id %in% id)){	
      stop("protein ids from disorder do not match with id argument")	
    }	
    seq_espritz$id <- factor(seq_espritz$id,	
                             levels = unique(dat$id))	
    seq_espritz$id_num <- as.numeric(seq_espritz$id)	
    seq_espritz$start[is.na(seq_espritz$start)] <- 0	
    seq_espritz$end[is.na(seq_espritz$end)] <- 0	
    
    seq_espritz <- merge(seq_espritz,	
                         dat[,c(2,3)],	
                         by.x = "id",	
                         by.y = "id",	
                         all.y = TRUE,	
                         all.x = TRUE,	
                         sort = FALSE) 	
    
    seq_espritz <- split(seq_espritz ,	
                         seq_espritz$id)	
    
    seq_espritz <- lapply(seq_espritz, function(x){	
      start <-  x$start	
      end <- x$end	
      id <- unique(x$id)	
      id_num <- unique(x$id_num)	
      resD <- lapply(seq_along(start), function(i){	
        start[i]:end[i]	
      })	
      groups <- rep(seq_along(resD), times = lengths(resD))	
      data.frame(residue = unlist(resD),	
                 groups = groups,	
                 id = id,	
                 id_num = id_num,	
                 nchar = unique(x$nchar))	
    })	
    
    rle_predO <- lapply(seq_espritz, function(x){	
      n <- unique(x$nchar)	
      n <- 1:n	
      n <- setdiff(n, x$residue)	
      if(length(n) != 0){	
        s <- split(n, cumsum(c(0, diff(n) != 1)))	
        groups <- rep(seq_along(s), times = lengths(s))	
      } else {	
        n <- NA	
        groups <- 1	
      }	
      rle_dat <- data.frame(residue = n,	
                            id = unique(x$id),	
                            id_num = unique(x$id_num),	
                            group = groups)	
      
      rle_dat <- split(rle_dat, rle_dat$group)	
      rle_dat <- lapply(rle_dat, function(x){	
        data.frame(start = min(x$residue),	
                   end = max(x$residue),	
                   id = unique(x$id),	
                   id_num = unique(x$id_num))	
      })	
      rle_dat <- do.call(rbind, rle_dat)	
      rle_dat	
    })	
    
    rle_predO <- do.call(rbind, rle_predO)	
    rownames(rle_predO) <- NULL	
    
    plot_segments <- rle_predO[!is.na(rle_predO$start),]	
    
    seq_espritz <- do.call(rbind, seq_espritz) 	
    
    seq_espritz <- seq_espritz[seq_espritz$residue != 0,]	
    if (nrow(seq_espritz) != 0){	
      rownames(seq_espritz) <- NULL	
      seq_espritz$y <- rep(c(0, -0.1, 0, 0.1),	
                           length.out = nrow(seq_espritz))	
      seq_espritz$y_plot <- seq_espritz$y + seq_espritz$id_num	
      seq_espritz$inter <- interaction(seq_espritz$group,	
                                       seq_espritz$id_num)
    } else {	
      seq_espritz <- NULL	
    }	
  }
  
  if (!is.null(seq_gpi)){
    for_plot <- merge(x = dat,
                      y = seq_gpi,
                      by = "id",
                      all.x = TRUE,
                      all.y = FALSE,
                      sort = FALSE)
  } else {
    for_plot <- dat
  }
  
  if (!is.null(seq_signalp)){
    for_plot <- merge(for_plot,
                      seq_signalp,
                      by = "id",
                      all.x = TRUE,
                      all.y = FALSE,
                      sort = FALSE)
  } else {
    for_plot <- for_plot
  }
  
  if (!is.null(dim(seq_hyp))){
    for_plot <- merge(for_plot,
                      seq_hyp,
                      by = "id",
                      all.x = TRUE,
                      all.y = FALSE,
                      sort = FALSE)
  } else {
    for_plot <- for_plot 
  }
  
  if (!is.null(dim(seq_scan))){
    for_plot <- merge(for_plot,
                      seq_scan,
                      by = "id",
                      all.x = TRUE,
                      all.y = FALSE,
                      sort = FALSE)
  } else {
    for_plot <- for_plot
  }
  
  lims <- c("N-sp", "omega site (gpi)", "hyp", "AG span", "TM")
  vals <- c(nsp_col, gpi_col, hyp_col, ag_col, tm_col)
  labs <- c("N-sp", "omega site (gpi)", "hyp", "AG span", "TM")
  
  subs <- c(!is.null(seq_signalp),
            !is.null(seq_gpi),
            !is.null(seq_hyp),
            !is.null(seq_scan),
            !is.null(phobius_seq))
  
  if (any(subs)){
    lims <- lims[subs]
    vals <- vals[subs]
    labs <- labs[subs]
  } else {
    lims <- NA
    vals <- NA
    labs <- NULL
  }
  rat <- max(nchar(sequence))
  rat <- round(rat/100, 0)*100
  rat <- rat/20
  
  if(!is.null(seq_espritz)){
    p <- ggplot2::ggplot(plot_segments)+
      ggplot2::geom_segment(ggplot2::aes_(y = ~id_num,
                                          yend = ~id_num,
                                          x = ~start,
                                          xend = ~end))
  } else {
    p <- ggplot2::ggplot(for_plot)+
      ggplot2::geom_segment(ggplot2::aes_(y = ~id_num,
                                          yend = ~id_num,
                                          x = 1,
                                          xend = ~nchar))
  }

  p <- p  +
    ggplot2::xlab("residue") +
    ggplot2::ylab("id") +
    ggplot2::scale_y_continuous(breaks = seq_along(id),
                                labels = id_label,
                                limits = c(0.5, length(id) + 0.5)) +
    ggplot2::theme_bw()+
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_line()) +
    ggplot2::coord_equal(ratio = rat,
                         expand = FALSE)
  
  if(!is.null(phobius_seq)){
    p <- p +
      ggplot2::geom_rect(data = tm,
                         ggplot2::aes_(ymin = ~id_tm - 0.35,
                                       ymax = ~id_tm + 0.35,
                                       xmin = ~tm_start,
                                       xmax = ~tm_end, color = "TM"),
                         fill = tm_col,
                         na.rm = TRUE) +
      ggplot2::geom_segment(data = inside,
                            ggplot2::aes_(y = ~id_inside - 0.35,
                                          yend = ~id_inside - 0.35,
                                          x = ~inside_start,
                                          xend = ~inside_end),
                            lty = 2,
                            na.rm = TRUE)+
      ggplot2::geom_segment(data = out,
                            ggplot2::aes_(y = ~id_out + 0.35,
                                          yend = ~id_out + 0.35,
                                          x = ~out_start,
                                          xend = ~out_end),
                            lty = 2,
                            na.rm = TRUE)
  }
  
  if(!is.null(seq_hmm)) {
    p <- p + 
      ggplot2::geom_rect(data = seq_hmm,
                         ggplot2::aes_(ymin = ~id_num - 0.25,
                                       ymax = ~id_num + 0.25,
                                       xmin = ~align_start,
                                       xmax = ~align_end,
                                       fill = ~domain),
                         colour = "grey20",
                         na.rm = TRUE)
  }
  
  
  if(!is.null(seq_signalp)) {
    p <- p +
      ggplot2::geom_segment(data = for_plot,
                            ggplot2::aes_(y = ~id_num,
                                          yend = ~id_num,
                                          x = 1,
                                          xend = ~sp.length,
                                          color = "N-sp"),
                            size = 2,
                            na.rm = TRUE)
  } 
  
  if(!is.null(seq_scan)) {
    p <- p +
      ggplot2::geom_rect(data = for_plot,
                         ggplot2::aes_(ymin = ~id_num - 0.18,
                                       ymax = ~id_num + 0.18,
                                       xmin = ~location.start,
                                       xmax = ~location.end,
                                       color = "AG span"),
                         fill = ag_col,
                         na.rm = TRUE) 
  } 
  
  if(!is.null(seq_hyp)) {
    p <- p +
      ggplot2::geom_errorbar(data = for_plot,
                             ggplot2::aes_(ymin = ~id_num - 0.18,
                                           ymax = ~id_num + 0.18,
                                           x  = ~P_pos,
                                           color = "hyp"),
                             size = 0.2,
                             width = 0,
                             na.rm = TRUE)
  }
  
  if(!is.null(seq_espritz)){	
    p <- p + ggplot2::geom_path(data = seq_espritz,	
                                ggplot2::aes_(x = ~residue,	
                                              y = ~y_plot,	
                                              group = ~inter),	
                                size = 0.7,	
                                color = "grey25") 	
  }

  if(!is.null(seq_gpi)) {
    p <- p +
      ggplot2::geom_point(data = for_plot,
                          ggplot2::aes_(y = ~id_num,
                                        x = ~omega_site,
                                        color = "omega site (gpi)"),
                          shape = gpi_shape,
                          size = gpi_size,
                          na.rm = TRUE) 
    
  } 
  
  if (any(subs)){
    p <- p +
      ggplot2::scale_color_manual("feature",
                                  limits = lims,
                                  values = vals,
                                  labels = labs)
    p <- p + 
      ggplot2::guides(color = ggplot2::guide_legend(keywidth = 1,
                                                    keyheight = 1,
                                                    override.aes = list(size = 5,
                                                                        shape = 15,
                                                                        linetype = 1)))
  }
  return(p)
}
