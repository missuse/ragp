#' Query SignalP web server.
#'
#' SignalP 4.1 server predicts the presence and location of signal peptide cleavage sites in amino acid sequences from different organisms: Gram-positive prokaryotes, Gram-negative prokaryotes, and eukaryotes. The method incorporates a prediction of cleavage sites and a signal peptide/non-signal peptide prediction based on a combination of several artificial neural networks.
#'
#' @aliases get_signalp get_signalp.default get_signalp.character get_signalp.data.frame get_signalp.list
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class "SeqFastaAA" resulting from \code{\link[seqinr]{read.fasta}} call. Should be left blank if vectors are provided to sequence and id arguments.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param org_type One of c("euk", "gram-", "gram+"), defaults to "euk". Which model should be used for prediction.
#' @param Dcut_type One of c("default", "sensitive", "user"), defaults to "default". The default cutoff values for SignalP 4 are chosen to optimize the performance measured as Matthews Correlation Coefficient (MCC). This results in a lower sensitivity (true positive rate) than SignalP 3.0 had. Setting this argument to "sensitive" will yield the same sensitivity as SignalP 3.0. This will make the false positive rate slightly higher, but still better than that of SignalP 3.0.
#' @param Dcut_noTM A numeric value, with range 0 - 1, defaults to 0.45. For experimenting with cutoff values.
#' @param Dcut_TM A numeric value, with range 0 - 1, defaults to 0.5. For experimenting with cutoff values.
#' @param method One of c("best", "notm"), defaults to "best". Signalp 4.1 contains two types of neural networks. SignalP-TM has been trained with sequences containing transmembrane segments in the data set, while SignalP-noTM has been trained without those sequences. Per default, SignalP 4.1 uses SignalP-TM as a preprocessor to determine whether to use SignalP-TM or SignalP-noTM in the final prediction (if 4 or more positions are predicted to be in a transmembrane state, SignalP-TM is used, otherwise SignalP-noTM). An exception is Gram-positive bacteria, where SignalP-TM is used always. If you are confident that there are no transmembrane segments in your data, you can get a slightly better performance by choosing "Input sequences do not include TM regions", which will tell SignalP 4.1 to use SignalP-noTM always.
#' @param minlen An integer value corresponding to the minimal predicted signal peptide length, at default set to 10. SignalP 4.0 could, in rare cases, erroneously predict signal peptides shorter than 10 residues. These errors have in SignalP 4.1 been eliminated by imposing a lower limit on the cleavage site position (signal peptide length). The minimum length is by default 10, but you can adjust it. Signal peptides shorter than 15 residues are very rare. If you want to disable this length restriction completely, enter 0 (zero).
#' @param trunc An integer value corresponding to the N-terminal truncation of input sequence, at default set to 70. By default, the predictor truncates each sequence to max. 70 residues before submitting it to the neural networks. If you want to predict extremely long signal peptides, you can try a higher value, or disable truncation completely by entering 0 (zero).
#' @param splitter An integer indicating the number of sequences to be in each .fasta file that is to be sent to the server. Defaults to 500. Change only in case of a server side error. Accepted values are in range of 1 to 2000.
#' @param sleep A numeric indicating the pause in seconds between POST and GET server calls, at default set to 1s. Decreasing is not recommended.
#' @param attempts Integer, number of attempts if server unresponsive, at default set to 2.
#' @param progress Boolean, whether to show the progress bar, at default set to FALSE.
#' @param ... currently no additional arguments are accepted apart the ones documented bellow.
#'
#' @return  A data frame with columns:
#' \describe{
#'   \item{id}{Character, as from input}
#'   \item{Cmax}{Numeric, C-score (raw cleavage site score). The output from the CS networks, which are trained to distinguish signal peptide cleavage sites from everything else. Note the position numbering of the cleavage site: the C-score is trained to be high at the position immediately after the cleavage site (the first residue in the mature protein).}
#'   \item{Cmax.pos}{Integer, position of Cmax. position immediately after the cleavage site (the first residue in the mature protein).}
#'   \item{Ymax}{Numeric, Y-score (combined cleavage site score), A combination (geometric average) of the C-score and the slope of the S-score, resulting in a better cleavage site prediction than the raw C-score alone. This is due to the fact that multiple high-peaking C-scores can be found in one sequence, where only one is the true cleavage site. The Y-score distinguishes between C-score peaks by choosing the one where the slope of the S-score is steep.}
#'   \item{Ymax.pos}{Integer, position of Ymax}
#'   \item{Smax}{Numeric, S-score (signal peptide score). The output from the SP networks, which are trained to distinguish positions within signal peptides from positions in the mature part of the proteins and from proteins without signal peptides.}
#'   \item{Smax.pos}{Integer, position of Smax}
#'   \item{Smean}{Numeric, The average S-score of the possible signal peptide (from position 1 to the position immediately before the maximal Y-score)}
#'   \item{Dmean}{Numeric, D-score (discrimination score). A weighted average of the mean S and the max. Y scores. This is the score that is used to discriminate signal peptides from non-signal peptides.}
#'   \item{is.sp}{Character, does the sequence contain a N-sp}
#'   \item{Dmaxcut}{Numeric, as from input, Dcut_noTM if SignalP-noTM network used and Dcut_TM if SignalP-TM network used}
#'   \item{Networks.used}{Character, which network was used for the prediction: SignalP-noTM or SignalP-TM}
#'   \item{is.signalp}{Logical, did SignalP predict the presence of a signal peptide}
#'   }
#'
#' @note This function creates temporary files in the working directory.
#'
#' @source \url{http://www.cbs.dtu.dk/services/SignalP-4.1/}
#' @references Petersen TN. Brunak S. Heijne G. Nielsen H. (2011) SignalP 4.0: discriminating signal peptides from transmembrane regions. Nature Methods 8: 785-786
#'
#' @seealso \code{\link[ragp]{get_phobius}} \code{\link[ragp]{get_targetp}}
#'
#' @examples
#' library(ragp)
#' signalp_pred <- get_signalp(data = at_nsp[1:10,],
#'                             sequence,
#'                             Transcript.id)
#' signalp_pred
#'
#' @import seqinr
#' @import httr
#' @import xml2
#' @export 

get_signalp <- function (data, ...){
  if (missing(data) || is.null(data)) get_signalp.default(...)
  else UseMethod("get_signalp")
}

#' @rdname get_signalp
#' @method get_signalp character
#' @export

get_signalp.character <- function(data,
                                  org_type = c("euk", "gram-", "gram+"),
                                  Dcut_type = c("default", "sensitive", "user"),
                                  Dcut_noTM = 0.45,
                                  Dcut_TM = 0.5,
                                  method = c("best", "notm"),
                                  minlen = NULL,
                                  trunc = 70L,
                                  splitter = 500L,
                                  sleep = 3,
                                  attempts = 2,
                                  progress = FALSE,
                                  ...){
  if (missing(splitter)) {
    splitter <- 500L
  }
  if (length(splitter) > 1){
    splitter <- 500L
    warning("splitter should be of length 1, setting to default: splitter = 500",
            call. = FALSE)
  }
  if (!is.numeric(splitter)){
    splitter <- as.numeric(splitter)
    warning("splitter is not numeric, converting using 'as.numeric'",
            call. = FALSE)
  }
  if (is.na(splitter)){
    splitter <- 500L
    warning("splitter was set to NA, setting to default: splitter = 500",
            call. = FALSE)
  }
  if (is.numeric(splitter)) {
    splitter <- floor(splitter)
  }
  if (!(splitter %in% 1:2000)) {
    splitter <- 500L
    warning("Illegal splitter input, splitter will be set to 500",
            call. = FALSE)
  }
  if (!missing(trunc)){
    if (length(trunc) > 1){
      stop("trunc should be of length 1.",
           call. = FALSE)
    }
    if (!is.numeric(trunc)){
      stop("trunc is not numeric.",
           call. = FALSE)
    }
    if (is.na(trunc)){
      stop("trunc was set to NA.",
           call. = FALSE)
    }
    if (is.numeric(trunc)){
      trunc <- floor(trunc)
    }
    if (trunc < 0){
      stop("trunc was set to a negative number.",
           call. = FALSE)
    }
    if (trunc == 0){
      trunc <- 1000000L
    }
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
  if (missing(sleep)) {
    sleep <- 3
  }
  if (length(sleep) > 1){
    sleep <- 3
    warning("sleep should be of length 1, setting to default: sleep = 3",
            call. = FALSE)
  }
  if (!is.numeric(sleep)){
    sleep <- as.numeric(sleep)
    warning("sleep is not numeric, converting using 'as.numeric'",
            call. = FALSE)
  }
  if (is.na(sleep)){
    sleep <- 3
    warning("sleep was set to NA, setting to default: sleep = 3",
            call. = FALSE)
  }
  if (sleep < 2){
    warning("setting sleep to less than 2s can cause problems when fetching results from the server",
            call. = FALSE)
  }
  if (missing(org_type)) {
    org_type <- "euk"
  }
  if (!org_type %in% c("euk", "gram-", "gram+")) {
    stop("org_type should be one of: 'euk', 'gram-', 'gram+'",
         call. = FALSE)
  }
  if (length(org_type) > 1){
    stop("org_type should be one of: 'euk', 'gram-', 'gram+'",
         call. = FALSE)
  }
  if (missing(Dcut_type)) {
    Dcut_type <- "default"
  }
  if (!Dcut_type %in% c("default", "sensitive", "user")) {
    stop("Dcut_type should be one of: 'default', 'sensitive', 'user'",
         call. = FALSE)
  }
  if (length(Dcut_type) > 1){
    stop("Dcut_type should be one of: 'default', 'sensitive', 'user'",
         call. = FALSE)
  }
  if (missing(Dcut_noTM)) {
    Dcut_noTM <- "0.45"
  }  else {
    Dcut_noTM <- as.character(Dcut_noTM)[1]
  }
  if (!is.numeric(as.numeric(Dcut_noTM))){
    Dcut_noTM <- "0.45"
    warning("Dcut_noTM could not be converted to numeric, setting to default: Dcut_noTM = '0.45'",
            call. = FALSE)
  }
  if (is.na(Dcut_noTM)) {
    Dcut_noTM <- "0.45"
    warning("Dcut_noTM was set to NA, setting to default: Dcut_noTM = '0.45'",
            call. = FALSE)
  }
  if (as.numeric(Dcut_noTM[1]) > 1) {
    Dcut_noTM <- "0.45"
    warning("Dcut_noTM must take values in the range 0 - 1,
            it was set to the default: Dcut_noTM = '0.45'",
            call. = FALSE)
  }
  if (as.numeric(Dcut_noTM[1]) < 0) {
    Dcut_noTM <- "0.45"
    warning("Dcut_noTM must take values in the range 0 - 1,
            it was set to the default: Dcut_noTM = '0.45'",
            call. = FALSE)
  }
  if (missing(Dcut_TM)) {
    Dcut_TM <- "0.5"
  } else {
    Dcut_TM <- as.character(Dcut_TM)[1]
  }
  if (!is.numeric(as.numeric(Dcut_TM))){
    Dcut_TM <- "0.5"
    warning("Dcut_TM could not be converted to numeric, setting to default: Dcut_TM = '0.5'",
            call. = FALSE)
  }
  if (is.na(Dcut_TM)) {
    Dcut_TM <- "0.5"
    warning("Dcut_noTM was set to NA, setting to default: Dcut_TM = '0.5'",
            call. = FALSE)
  }
  if (as.numeric(Dcut_TM[1]) > 1) {
    Dcut_TM <- "0.5"
    warning("Dcut_TM must take values in the range 0 - 1,
            it was set to the default: Dcut_TM = '0.5'",
            call. = FALSE)
  }
  if (as.numeric(Dcut_TM[1]) < 0) {
    Dcut_TM <- "0.5"
    warning("Dcut_TM must take values in the range 0 - 1,
            it was set to the default: Dcut_TM = '0.5'",
            call. = FALSE)
  }
  if (missing(method)) {
    method <- "best"
  }
  if (!method %in% c("best", "notm")){
    stop("method should be one of: 'best', 'notm'",
         call. = FALSE)
  }
  if (length(method) > 1){
    stop("method should be one of: 'best', 'notm'",
         call. = FALSE)
  }
  if (missing(minlen)) {
    minlen <- ""
  }  else {
    minlen <- as.character(minlen)[1]
  }
  if(length(data) > 1){
    stop("one fasta file per function call can be supplied",
         call. = FALSE)
  }
  if (file.exists(data)){
    file_name <- data
  } else {
    stop("cannot find file in the specified path",
         call. = FALSE)
  }
  url <- "http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi"
  cfg_file <- "/usr/opt/www/pub/CBS/services/SignalP-4.1/SignalP.cf"
  file_list <- ragp::split_fasta(path_in = file_name,
                                 path_out = "tmp_signalp_",
                                 num_seq = splitter,
                                 trunc = trunc)
  if(grepl("temp_", file_name)){
    unlink(file_name)
  }
  for_pb <- length(file_list)
  if(progress){
    pb <- utils::txtProgressBar(min = 0,
                                max = for_pb,
                                style = 3)
  }
  splt <- (seq_along(file_list) - 1) %/% 10
  file_list <- split(file_list,
                     splt)
  output <- vector("list", length(file_list)*10)
  for(k in seq_along(file_list)){
    x <- file_list[[k]]
    jobid <- vector("character", 10)
    for (i in seq_along(x)) {
      file_up <- httr::upload_file(x[i])
      if (trunc == 1000000L){
        trunc <- ""
      }
      res <- httr::POST(url = url,
                        encode = "multipart",
                        body = list(configfile = cfg_file,
                                    SEQSUB = file_up,
                                    orgtype = org_type,
                                    `Dcut-type` = Dcut_type,
                                    `Dcut-noTM` = Dcut_noTM,
                                    `Dcut-TM` = Dcut_TM,
                                    graphmode = NULL,
                                    format = "short",
                                    minlen = minlen,
                                    method = method,
                                    trunc = as.character(trunc)))
      if(!grepl("jobid=", res$url)){
        stop("something went wrong on server side")
      }
      res <- sub("http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi?jobid=",
                 "",
                 res$url,
                 fixed = TRUE)
      
      res <- sub("&wait=20",
                 "",
                 res,
                 fixed = TRUE)
      jobid[i] <- res
      if(progress){
        utils::setTxtProgressBar(pb,
                                 floor(i/2) + (10 * (k - 1)))
      }
      #Sys.sleep(sleep)
    }
    
    collected_res <- vector("list", length(x))
    for (i in seq_along(x)) {
      time1 <- Sys.time()
      repeat {
        res2 <- httr::GET(url = url,
                          query = list(jobid = jobid[i],
                                       wait = "20"))
        bad <- xml2::xml_text(
          xml2::xml_find_all(
            httr::content(res2,
                          as = "parsed"),
            "//head")
        )
        if (grepl("Illegal", bad)) {
          prt <- xml2::xml_text(
            xml2::xml_find_all(
              httr::content(res2,
                            as = "parsed"),
              "//li")
          )
          stop(paste0(prt, ". Problem in file: ", "temp_",
                      i, ".fa"),
               call. = FALSE)
        }
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
        Sys.sleep(1)
        if (any(grepl("Cmax", res2_split))) {
          break
        }
        
        time2 <- Sys.time()
        
        max.time <- as.difftime(pmax(50, splitter),
                                units = "secs")
        
        if ((time2 - time1) > max.time) {
          res2_split <- NULL
          if(progress) message(
            "file",
            x[i],
            "took longer then expected")
          break
        }
      }
      if (is.null(res2_split)) {
        tms <- 0
        while(tms < attempts && is.null(res2_split)){
          if(progress) message(
            "reattempting file",
            x[i])
          file_up <-  httr::upload_file(x[i])
          res <- httr::POST(url = url,
                            encode = "multipart",
                            body = list(configfile = cfg_file,
                                        SEQSUB = file_up,
                                        orgtype = org_type,
                                        `Dcut-type` = Dcut_type,
                                        `Dcut-noTM` = Dcut_noTM,
                                        `Dcut-TM` = Dcut_TM,
                                        graphmode = NULL,
                                        format = "short",
                                        minlen = minlen,
                                        method = method,
                                        trunc = as.character(trunc)))
          if(!grepl("jobid=", res$url)){
            stop("something went wrong on server side")
          }
          res <- sub("http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi?jobid=",
                     "",
                     res$url,
                     fixed = TRUE)
          
          res <- sub("&wait=20",
                     "",
                     res,
                     fixed = TRUE)
          jobidi <- res
          
          time1 <- Sys.time()
          
          repeat {
            res2 <- httr::GET(url = url,
                              query = list(jobid = jobidi,
                                           wait = "20"))
            bad <- xml2::xml_text(
              xml2::xml_find_all(
                httr::content(res2,
                              as = "parsed"),
                "//head")
            )
            if (grepl("Illegal", bad)) {
              prt <- xml2::xml_text(
                xml2::xml_find_all(
                  httr::content(res2,
                                as = "parsed"),
                  "//li")
              )
              stop(paste0(prt,
                          ". Problem in file: ",
                          "temp_",
                          i,
                          ".fa"),
                   call. = FALSE)
            }
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
            Sys.sleep(1)
            if (any(grepl("Cmax", res2_split))) {
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
        output$is.signalp <- output$is.sp == "Y"
        if(progress){
          utils::setTxtProgressBar(pb,
                                   for_pb)
          close(pb)
        }
        warning(
          "maximum attempts reached at",
          x[i],
          "returning finished queries",
          call. = FALSE)
        return(output)
      }
      unlink(x[i])
      
      res2_split <- res2_split[(which(grepl("name",
                                            res2_split))[1] +
                                  1):(which(grepl("/pre",
                                                  res2_split ))[1] - 1)]
      
      if(any(grepl("hr", res2_split))){
        res2_split <- res2_split[1:(which(grepl("<hr>",
                                                res2_split))[1] - 1)]
      }
      res2_split <- strsplit(res2_split,
                             " +")
      res2_split <- do.call(rbind,
                            res2_split)
      res2_split <- as.data.frame(res2_split,
                                  stringsAsFactors = F)
      colnames(res2_split) <- c("id",
                                "Cmax",
                                "Cmax.pos",
                                "Ymax",
                                "Ymax.pos",
                                "Smax",
                                "Smax.pos",
                                "Smean",
                                "Dmean",
                                "is.sp",
                                "Dmaxcut",
                                "Networks.used")
      if(progress){
        utils::setTxtProgressBar(pb,
                                 floor(i/2) + 5 + (10 * (k - 1)))
      }
      output[[((k*10)-10)+i]] <- res2_split
    }
  }
  if(progress){
    utils::setTxtProgressBar(pb,
                             for_pb)
    close(pb)
  }
  
  output <- do.call(rbind,
                    output)
  
  output$is.signalp <- output$is.sp == "Y"
  return(output)
}

#' @rdname get_signalp
#' @method get_signalp data.frame
#' @export

get_signalp.data.frame <- function(data,
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
  res <- get_signalp.character(data = file_name, ...)
  return(res)
}

#' @rdname get_signalp
#' @method get_signalp list
#' @export


get_signalp.list <- function(data,
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
  res <- get_signalp.character(data = file_name, ...)
  return(res)
}

#' @rdname get_signalp
#' @method get_signalp default
#' @export

get_signalp.default <- function(data = NULL,
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
  res <- get_signalp.character(data = file_name, ...)
  return(res)
}


