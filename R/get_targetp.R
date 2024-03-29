#' Query TargetP web server.
#'
#' TargetP 1.1 predicts the subcellular location of eukaryotic proteins. The location assignment is based on the predicted presence of any of the N-terminal presequences: chloroplast transit peptide (cTP), mitochondrial targeting peptide (mTP) or secretory pathway signal peptide (SP). TargetP uses ChloroP and SignalP to predict cleavage sites for cTP and SP, respectively. For the sequences predicted to contain an N-terminal presequence a potential cleavage site is also predicted.
#'
#' @aliases get_targetp get_targetp.default get_targetp.character get_targetp.data.frame get_targetp.list get_targetp.AAStringSet
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class \code{\link[seqinr]{SeqFastaAA}} resulting from \code{\link[seqinr]{read.fasta}} call. Alternatively an \code{\link[Biostrings]{AAStringSet}} object. Should be left blank if vectors are provided to sequence and id arguments.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param org_type One of c("non_plant", "plant"), defaults to "plant". Which models should be used for prediction.
#' @param cutoffs One of c("winner_takes_all", "spec95", "spec90", "custom"), defaults to "winner_takes_all". If "winner_takes_all" no cutoffs are specified, if "spec95" or "spec90" are selected, predefined set of cutoffs that yielded >0.95 or >0.9 specificity on the TargetP test sets. If "custom" specified user defined cutoffs should be specified in "tcut", "pcut", "scut", ocut".
#' @param tcut A numeric value, with range 0 - 1, defaults to 0 (cutoff = "winner_takes_all"). mTP user specified cutoff.
#' @param pcut A numeric value, with range 0 - 1, defaults to 0 (cutoff = "winner_takes_all"). cTP user specified cutoff.
#' @param scut A numeric value, with range 0 - 1, defaults to 0 (cutoff = "winner_takes_all"). SP user specified cutoff.
#' @param ocut A numeric value, with range 0 - 1, defaults to 0 (cutoff = "winner_takes_all"). User specified cutoff for "other" (not with mTP, cTP, SP).
#' @param splitter An integer indicating the number of sequences to be in each .fasta file that is to be sent to the server. Defaults to 1000. Change only in case of a server side error. Accepted values are in range of 1 to 2000.
#' @param attempts Integer, number of attempts if server unresponsive, at default set to 2.
#' @param progress Boolean, whether to show the progress bar, at default set to FALSE.
#' @param ... currently no additional arguments are accepted apart the ones documented bellow.
#'
#' @return  A data frame with columns:
#' \describe{
#'   \item{Name}{Character, sequence name truncated to 20 characters}
#'   \item{Len}{Integer, length of analyzed sequence}
#'   \item{cTP}{Numeric, final NN sequence score to contain a chloroplast transit peptide (cTP). ChloroP is used to predict cleavage sites for cTP}
#'   \item{mTP}{Numeric, final NN sequence score to contain a mitochondrial targeting peptide (mTP).}
#'   \item{SP}{Numeric, final NN sequence score to contain a secretory pathway signal peptide (SP). SignalP is used to predict cleavage sites for SP}
#'   \item{other}{Numeric, final NN sequence score of a sequence not to contain mTP, SP or cTP}
#'   \item{Loc}{Character, one of C (chloroplast), M (mitochondrion), S (secretory pathway), - (any other location) or * ("don't know"; indicates that cutoff restrictions were set and the winning network output score was below the requested cutoff for that category.) }
#'   \item{TPlen}{Integer, predicted presequence length}
#'   \item{RC}{Integer, reliability class, from 1 to 5, where 1 indicates the strongest prediction. RC is a measure of the size of the difference ('diff') between the highest (winning) and the second highest output scores. There are 5 reliability classes, defined as follows: 1 : diff > 0.800, 2 : 0.800 > diff > 0.600, 3 : 0.600 > diff > 0.400, 4 : 0.400 > diff > 0.200, 5 : 0.200 > diff. Thus, the lower the value of RC the safer the prediction.}
#'   \item{is.targetp}{Logical, did TargetP predict the presence of a signal peptide}
#'   }
#'
#' @note This function creates temporary files in the working directory. Protein ids should be shorter then 20 characters due to server side truncation.
#'
#' @source \url{https://services.healthtech.dtu.dk/service.php?TargetP-1.1}
#' @references Emanuelsson O, Nielsen H, Brunak S,von Heijne G. (2000) Predicting subcellular localization of proteins based on their N-terminal amino acid sequence. J. Mol. Biol.300: 1005-1016
#'
#' @seealso \code{\link[ragp]{get_signalp}} \code{\link[ragp]{get_signalp5}}
#'
#' @examples
#' library(ragp)
#' data(at_nsp)
#' targetp_pred <- get_targetp(at_nsp[1:20,],
#'                             sequence,
#'                             Transcript.id)
#' targetp_pred     
#'                       
#' @import seqinr
#' @import httr
#' @import xml2
#' @export

get_targetp <- function (data, ...){
  if (missing(data) || is.null(data)) get_targetp.default(...)
  else UseMethod("get_targetp")
}

#' @rdname get_targetp
#' @method get_targetp character
#' @export

get_targetp.character <-   function(data,
                                    org_type = c("non_plant", "plant"),
                                    cutoffs = c("winner_takes_all", "spec95", "spec90", "custom"),
                                    tcut = NULL,
                                    pcut = NULL,
                                    scut = NULL,
                                    ocut = NULL,
                                    splitter = 1000,
                                    attempts = 2,
                                    progress = FALSE,
                                    ...){
  if (missing(org_type)){
    org_type <- "plant"
  }
  if (!org_type %in% c("non_plant", "plant")) {
    stop("org_type should be one of: 'non_plant', 'plant'",
         call. = FALSE)
  }
  if (length(org_type) > 1){
    stop("org_type should be one of: 'non_plant', 'plant'",
         call. = FALSE)
  }
  if (missing(tcut)){
    tcut <- 0
  }
  if (missing(pcut)){
    pcut <- 0
  }
  if (missing(scut)){
    scut <- 0
  }
  if (missing(ocut)){
    ocut <- 0
  }
  if (missing(splitter)){
    splitter <- 1000
  }
  if (length(splitter) > 1){
    splitter <- 1000
    warning("splitter should be of length 1, setting to default: splitter = 1000",
            call. = FALSE)
  }
  if (!is.numeric(splitter)){
    splitter <- as.numeric(splitter)
    warning("splitter is not numeric, converting using 'as.numeric'",
            call. = FALSE)
  }
  if (is.na(splitter)){
    splitter <- 1000
    warning("splitter was set to NA, setting to default: splitter = 1000",
            call. = FALSE)
  }
  if (is.numeric(splitter)) {
    splitter <- floor(splitter)
  }
  if (!(splitter %in% 1:2000)){
    splitter <- 1000
    warning("Illegal splitter input, splitter will be set to 1000",
            call. = FALSE)
  }
  if (length(attempts) > 1){
    attempts <- 2L
    warning("attempts should be of length 1, setting to default: attempts = 2",
            call. = FALSE)
  }
  if (!is.numeric(attempts)){
    attempts<- as.numeric(attempts)
    warning("attempts is not numeric, converting using 'as.numeric'",
            call. = FALSE)
  }
  if (is.na(attempts)){
    attempts <- 2L
    warning("attempts was set to NA, setting to default: attempts = 2",
            call. = FALSE)
  }
  if (is.numeric(attempts)) {
    attempts <- floor(attempts)
  }
  if (attempts < 1) {
    attempts <- 2L
    warning("attempts was set to less then 1, setting to default: attempts = 2",
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
  if (missing(cutoffs)){
    cutoffs <- "winner_takes_all"
  }
  if (!cutoffs %in% c("winner_takes_all", "spec95", "spec90", "custom")){
    stop ("cutoffs can be one of: 'winner_takes_all', 'spec95', 'spec90', 'custom'",
          call. = FALSE)
  }
  if (length(cutoffs) > 1){
    stop("cutoffs can be one of: 'winner_takes_all', 'spec95', 'spec90', 'custom'",
         call. = FALSE)
  }
  if (cutoffs != "custom"){
    if (!is.null(tcut)){
      if (tcut != 0) {
        warning ("cutoffs were not set to custom, pcut will be held at 0",
                 call. = FALSE)
      }
    }
    if (!is.null(pcut)){
      if (pcut != 0) {
        warning ("cutoffs were not set to custom, pcut will be held at 0",
                 call. = FALSE)
      }
    }
    if (!is.null(scut)){
      if (scut != 0) {
        warning ("cutoffs were not set to custom, scut will be held at 0",
                 call. = FALSE)
      }
    }
    if (!is.null(ocut)){
      if (ocut != 0) {
        warning ("cutoffs were not set to custom, ocut will be held at 0",
                 call. = FALSE)
      }
    }
  }
  if (cutoffs != "custom"){
    tcut <- 0
    pcut <- 0
    scut <- 0
    ocut <- 0
  }
  if (cutoffs == "custom"){
    if (missing(tcut)){
      tcut <- 0
    }
    if (length(tcut) > 1){
      tcut <- 0
      warning("tcut should be of length 1, setting to default: tcut = 0",
              call. = FALSE)
    }
    if (!is.numeric(tcut)){
      tcut <- as.numeric(tcut)
      warning("tcut is not numeric, converting using 'as.numeric'",
              call. = FALSE)
    }
    if (is.na(tcut)){
      tcut <- 0
      warning("tcut was set to NA, setting to default: tcut = 0",
              call. = FALSE)
    }
    if (tcut < 0){
      tcut <- 0
      warning(paste("tcut should be in 0 - 1 range,",
                    "tcut was set to the default 0"),
              call. = FALSE)
    }
    if (tcut > 1){
      tcut <- 1
      warning(paste("tcut should be in 0 - 1 range,",
                    "tcut was set to 1"),
              call. = FALSE)
    }
    if (missing(pcut)){
      pcut <- 0
    }
    if (length(pcut) > 1){
      pcut <- 0
      warning("pcut should be of length 1, setting to default: pcut = 0",
              call. = FALSE)
    }
    if (!is.numeric(pcut)){
      pcut <- as.numeric(pcut)
      warning("pcut is not numeric, converting using 'as.numeric'",
              call. = FALSE)
    }
    if (is.na(pcut)){
      pcut <- 0
      warning("pcut was set to NA, setting to default: pcut = 0",
              call. = FALSE)
    }
    if (pcut < 0){
      pcut <- 0
      warning(paste("pcut should be in 0 - 1 range,",
                    "pcut was set to the default 0"),
              call. = FALSE)
    }
    if (pcut > 1){
      pcut <- 1
      warning(paste("pcut should be in 0 - 1 range,",
                    "pcut was set to 1"),
              call. = FALSE)
    }
    if (missing(scut)){
      scut <- 0
    }
    if (length(scut) > 1){
      scut <- 0
      warning("scut should be of length 1, setting to default: scut = 0",
              call. = FALSE)
    }
    if (!is.numeric(scut)){
      scut <- as.numeric(scut)
      warning("scut is not numeric, converting using 'as.numeric'",
              call. = FALSE)
    }
    if (is.na(scut)){
      scut <- 0
      warning("scut was set to NA, setting to default: scut = 0",
              call. = FALSE)
    }
    if (scut < 0){
      scut <- 0
      warning(paste("scut should be in 0 - 1 range,",
                    "scut was set to the default 0"),
              call. = FALSE)
    }
    if (scut > 1){
      scut <- 1
      warning(paste("scut should be in 0 - 1 range,",
                    "scut was set to 1"),
              call. = FALSE)
    }
    if (missing(ocut)){
      ocut <- 0
    }
    if (length(ocut) > 1){
      ocut <- 0
      warning("ocut should be of length 1, setting to default: ocut = 0",
              call. = FALSE)
    }
    if (!is.numeric(ocut)){
      ocut <- as.numeric(ocut)
      warning("ocut is not numeric, converting using 'as.numeric'",
              call. = FALSE)
    }
    if (is.na(ocut)){
      ocut <- 0
      warning("ocut was set to NA, setting to default: ocut = 0",
              call. = FALSE)
    }
    if (ocut < 0){
      ocut <- 0
      warning(paste("ocut should be in 0 - 1 range,",
                    "ocut was set to the default 0"),
              call. = FALSE)
    }
    if (ocut > 1){
      ocut <- 1
      warning(paste("ocut should be in 0 - 1 range,",
                    "ocut was set to 1"),
              call. = FALSE)
    }
  }
  if (cutoffs == "winner_takes_all"){
    spec <- 0
  }
  if (cutoffs == "spec95"){
    spec <- 1
  }
  if (cutoffs == "spec90"){
    spec <- 2
  }
  if (cutoffs == "custom"){
    spec <- 3
  }
  
  if (file.exists(data)){
    file_name <- data
  } else {
    stop("cannot find file in the specified path",
         call. = FALSE)
  }
  
  file_list <- ragp::split_fasta(path_in = file_name,
                                 path_out = "tmp_targetp_",
                                 num_seq = splitter,
                                 trim = TRUE)
  if(grepl("temp_", file_name)){
    unlink(file_name)
  }
  
  if(progress){
    for_pb <- length(file_list)
    pb <- utils::txtProgressBar(min = 0,
                                max = for_pb,
                                style = 3)
  }

  output <- vector("list", length(file_list))
  url <- "https://services.healthtech.dtu.dk/cgi-bin/webface2.fcgi"
  cfg <- "/var/www/html/services/TargetP-1.1/webface.cf"
  for(k in seq_along(file_list)){
    x <- file_list[[k]]
    file_up <-  httr::upload_file(x)
    
    res <- httr::POST(
      url = url,
      encode = "multipart",
      body = list(
        `configfile` = cfg,
        `SEQSUB` =  file_up,
        `orgtype` = org_type,
        `cleavsite` = "on",
        `spec` = spec,
        `tcut` = tcut,
        `pcut` = pcut,
        `scut` = scut,
        `ocut` = ocut
      ))
    
    if(!grepl("jobid=", res$url)){
      stop("something went wrong on server side")
    }
    
    res <- sub("https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?jobid=",
               "",
               res$url,
               fixed = TRUE)
    
    res <- sub("&wait=20",
               "",
               res,
               fixed = TRUE)
    
    jobid <- res

    if(progress){
      utils::setTxtProgressBar(pb,
                               k)
    }
      
    Sys.sleep(2)

    time1 <- Sys.time()
    
    repeat {
      res2 <- httr::GET(
        url = url,
        query = list(
          jobid = jobid,
          wait = "20"
        ))
      
      code <- res2$status_code
      
      if(code != 200){
        res2_split <- NULL
        warning(paste0( ". Problem in file: ",
                        x))
      } else {
        res2 <- as.character(
          xml2::xml_find_all(
            httr::content(res2,
                          as = "parsed"),
            ".//pre")
        )
        res2_split <- unlist(
          strsplit(res2,
                   "\n")
        )
      }
    
      Sys.sleep(1)
      
      if (any(grepl("mTP", res2_split))){
        
        break
      }

      time2 <- Sys.time()
      
      max.time <- as.difftime(pmax(50, splitter),
                              units = "secs")
      
      
      if ((time2 - time1) > max.time) {
        res2_split <- NULL
        if(progress) message(
          "file",
          x,
          "took longer then expected")
        break
      }
    }
    if (is.null(res2_split)) {
      tms <- 0
      while(tms < attempts && is.null(res2_split)){
        if(progress) message(
          "reattempting file",
          x)
        file_up <-  httr::upload_file(x)
        res <- httr::POST(
          url = url,
          encode = "multipart",
          body = list(
            `configfile` = cfg,
            `SEQSUB` =  file_up,
            `orgtype` = org_type,
            `cleavsite` = "on",
            `spec` = spec,        
            `tcut` = tcut,
            `pcut` = pcut,
            `scut` = scut,
            `ocut` = ocut
          ))
        
        if(!grepl("jobid=", res$url)){
          stop("something went wrong on server side")
        }
        
        res <- sub("https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?jobid=",
                   "",
                   res$url,
                   fixed = TRUE)
        
        res <- sub("&wait=20",
                   "",
                   res,
                   fixed = TRUE)
        
        jobid <- res
        
        time1 <- Sys.time()
        
        repeat {
          res2 <- httr::GET(
            url = url,
            query = list(
              jobid = jobid,
              wait = "20"
            ))
          
          code <- res2$status_code
          
          if(code != 200){
            res2_split <- NULL
            warning(paste0( ". Problem in file: ",
                            x))
          } else {
            res2 <- as.character(
              xml2::xml_find_all(
                httr::content(res2,
                              as = "parsed"),
                ".//pre")
            )
            res2_split <- unlist(
              strsplit(res2,
                       "\n")
            )
          }
          Sys.sleep(1)
          if (any(grepl("mTP", res2_split))){
            break
          }
          
          time2 <- Sys.time()
          
          max.time <- as.difftime(pmax(100, splitter * 1.5),
                                  units = "secs")
          
          if ((time2 - time1) > max.time) {
            res2_split <- NULL
            break
          }
        }
        
        tms <- tms + 1
      }
    }
    
    if (is.null(res2_split)){
      output <- do.call(rbind,
                        output)
      output$is.targetp <- output$Loc == "S"
      if(progress){
        utils::setTxtProgressBar(pb,
                                 for_pb)
        close(pb)
      }
      warning(
        "maximum attempts reached at",
        x,
        "returning finished queries",
        call. = FALSE)
      return(output)
    }
    unlink(x)
    res2_split <- res2_split[(which(grepl("mTP",
                                          res2_split))[1]+2):(which(grepl("cutoff",
                                                                          res2_split))[1] - 2)]
    res2_split <- strsplit(res2_split,
                           " +")
    
    res2_split <- do.call(rbind,
                          res2_split)
    
    res2_split <- as.data.frame(res2_split,
                                stringsAsFactors = FALSE)
    
    if(org_type == "plant"){
      colnames(res2_split) <- c("Name",
                                "Len",
                                "cTP",
                                "mTP",
                                "SP",
                                "other",
                                "Loc",
                                "RC",
                                "TPlen")
    }
    
    if(org_type == "non_plant"){
      colnames(res2_split) <- c("Name",
                                "Len",
                                "mTP",
                                "SP",
                                "other",
                                "Loc",
                                "RC",
                                "TPlen")
    }
    if(progress){
      utils::setTxtProgressBar(pb,
                               k)
    }
    output[[k]] <- res2_split
  }
  
  if(progress){
    utils::setTxtProgressBar(pb,
                             for_pb)
    close(pb)
  }
  output <- do.call(rbind,
                    output)
  output$is.targetp <- output$Loc == "S"
  return(output)
}

#' @rdname get_targetp
#' @method get_targetp data.frame
#' @export

get_targetp.data.frame <- function(data,
                                   sequence,
                                   id,
                                   ...){
  if(missing(sequence)){
    stop("the column name with the sequences must be specified",
         call. = FALSE)
  }
  if(missing(id)){
    stop("the column name with the sequence id's must be specified",
         call. = FALSE)
  }
  id <- as.character(substitute(id))
  sequence <- as.character(substitute(sequence))
  if (length(id) != 1L){
    stop("only one column name for 'id' must be specifed",
         call. = FALSE)
  }
  if (length(sequence) != 1L){
    stop("only one column name for 'sequence' must be specifed",
         call. = FALSE)
  }
  id <- if(id %in% colnames(data)){
    data[[id]]
  } else {
    stop("specified 'id' not found in data",
         call. = FALSE)
  }
  id <- as.character(id)
  sequence  <- if(sequence %in% colnames(data)){
    data[[sequence]]
  } else {
    stop("specified 'sequence' not found in data",
         call. = FALSE)
  }
  sequence <- toupper(as.character(sequence))
  sequence <- sub("\\*$",
                  "",
                  sequence)
  aa_regex <- "[^ARNDCQEGHILKMFPSTWYV]"
  if (any(grepl(aa_regex, sequence))){
    warning(paste("sequences: ",
                  paste(id[grepl(aa_regex,
                                 sequence)],
                        collapse = ", "),
                  " contain symbols not corresponding to amino acids",
                  sep = ""),
            call. = FALSE)
  }
  file_name <- paste("temp_",
                     gsub("^X",
                          "",
                          make.names(Sys.time())),
                     ".fasta",
                     sep = "")
  seqinr::write.fasta(sequence = strsplit(sequence, ""),
                      name = id,
                      file = file_name)
  res <- get_targetp.character(data = file_name, ...)
  return(res)
}

#' @rdname get_targetp
#' @method get_targetp list
#' @export


get_targetp.list <- function(data,
                             ...){
  if(class(data[[1]]) ==  "SeqFastaAA"){
    dat <- lapply(data,
                  paste0,
                  collapse ="")
    id <- names(dat)
    sequence <- toupper(as.character(unlist(dat)))
    sequence <- sub("\\*$",
                    "",
                    sequence)
    aa_regex <- "[^ARNDCQEGHILKMFPSTWYV]"
    if (any(grepl(aa_regex, sequence))){
      warning(paste("sequences: ",
                    paste(id[grepl(aa_regex,
                                   sequence)],
                          collapse = ", "),
                    " contain symbols not corresponding to amino acids",
                    sep = ""),
              call. = FALSE)
    }
    file_name <- paste("temp_",
                       gsub("^X",
                            "",
                            make.names(Sys.time())),
                       ".fasta",
                       sep = "")
    seqinr::write.fasta(sequence = strsplit(sequence, ""),
                        name = id,
                        file = file_name)
  } else {
    stop("only lists containing objects of class SeqFastaAA are supported")
  }
  res <- get_targetp.character(data = file_name, ...)
  return(res)
}

#' @rdname get_targetp
#' @method get_targetp default
#' @export

get_targetp.default <- function(data = NULL,
                                sequence,
                                id,
                                ...){
  if (missing(sequence)){
    stop("protein sequence must be provided to obtain predictions",
         call. = FALSE)
  }
  if (missing(id)){
    stop("protein id must be provided to obtain predictions",
         call. = FALSE)
  }
  id <- as.character(id)
  sequence <- toupper(as.character(sequence))
  if (length(sequence) != length(id)){
    stop("id and sequence vectors are not of same length",
         call. = FALSE)
  }
  sequence <- sub("\\*$",
                  "",
                  sequence)
  aa_regex <- "[^ARNDCQEGHILKMFPSTWYV]"
  if (any(grepl(aa_regex, sequence))){
    warning(paste("sequences: ",
                  paste(id[grepl(aa_regex,
                                 sequence)],
                        collapse = ", "),
                  " contain symbols not corresponding to amino acids",
                  sep = ""),
            call. = FALSE)
  }
  file_name <- paste("temp_",
                     gsub("^X",
                          "",
                          make.names(Sys.time())),
                     ".fasta",
                     sep = "")
  seqinr::write.fasta(sequence = strsplit(sequence, ""),
                      name = id,
                      file = file_name)
  res <- get_targetp.character(data = file_name, ...)
  return(res)
}

#' @rdname get_targetp
#' @method get_targetp AAStringSet
#' @export

get_targetp.AAStringSet <-  function(data,
                                     ...){
  sequence <- as.character(data)
  id <- names(sequence)
  sequence <- unname(sequence)
  sequence <- toupper(sequence)
  sequence <- sub("\\*$",
                  "",
                  sequence)
  
  res <- get_targetp.default(sequence = sequence,
                             id = id,
                             ...)
  return(res)
}

