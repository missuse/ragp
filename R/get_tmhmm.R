#' Query TMHMM-2.0 web server.
#'
#' TMHMM server offers prediction of transmembrane helices in proteins
#'
#' @aliases get_tmhmm get_tmhmm.default get_tmhmm.character get_tmhmm.data.frame get_tmhmm.list
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class \code{\link[seqinr]{SeqFastaAA}} resulting from \code{\link[seqinr]{read.fasta}} call. Alternatively an \code{\link[Biostrings]{AAStringSet}} object. Should be left blank if vectors are provided to sequence and id arguments.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param splitter An integer indicating the number of sequences to be in each .fasta file that is to be sent to the server. Default is 2500. Change only in case of a server side error. Accepted values are in range of 1 to 10000.
#' @param attempts Integer, number of attempts if server unresponsive, at default set to 2.
#' @param progress Boolean, whether to show messages of the job id for each batch. Default is FALSE
#' @param ... currently no additional arguments are accepted apart the ones documented bellow.
#' 
#' @return  A data frame with columns:
#' \describe{
#'   \item{id}{Character, name of the submitted sequence.}
#'   \item{length}{Integer, length of the protein sequence}
#'   \item{ExpAA}{Numeric, the expected number of amino acids in transmembrane helices.}
#'   \item{First60}{Numeric, the expected number of amino acids in transmembrane helices in the first 60 amino acids of the protein.}
#'   \item{tm}{Integer, the number of predicted transmembrane segments.}
#'   \item{prediction}{Character string, predicted topology of the protein.}}
#'
#' @note This function creates temporary files in the working directory. If something goes wrong during communication with the server and progress was set to TRUE, predictions can be obtained using the web address `paste("https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?jobid=", jobid, "&wait=20", sep = "")`.
#'
#' @source \url{https://services.healthtech.dtu.dk/service.php?TMHMM-2.0}
#' @references Krogh A, Larsson B, von Heijne G, Sonnhammer EL (2001) Predicting transmembrane protein topology with a hidden Markov model: application to complete genomes. J Mol Biol 305(3):567-80.
# 
#' @seealso \code{\link[ragp]{get_phobius}} 
#'
#' @examples
#' library(ragp)
#' tmhmm_pred <- get_tmhmm(data = at_nsp[1:10,],
#'                         sequence,
#'                         Transcript.id)
#' tmhmm_pred
#'
#' @import seqinr
#' @import httr
#' @import xml2
#' @export

get_tmhmm <- function(data, ...){
  if (missing(data) || is.null(data)) get_tmhmm.default(...)
  else UseMethod("get_tmhmm")
}

#' @rdname get_tmhmm
#' @method get_tmhmm character
#' @export

get_tmhmm.character <- function(data,
                                splitter = 2500L,
                                attempts = 2,
                                progress = FALSE,
                                ...){
  if (missing(splitter)) {
    splitter <- 2500L
  }
  if (length(splitter) > 1){
    splitter <- 2500L
    warning("splitter should be of length 1, setting to default: splitter = 2500",
            call. = FALSE)
  }
  if (!is.numeric(splitter)){
    splitter <- as.numeric(splitter)
    warning("splitter is not numeric, converting using 'as.numeric'",
            call. = FALSE)
  }
  if (is.na(splitter)){
    splitter <- 2500L
    warning("splitter was set to NA, setting to default: splitter = 2500",
            call. = FALSE)
  }
  if (is.numeric(splitter)) {
    splitter <- floor(splitter)
  }
  if (!(splitter %in% 1:5000)) {
    splitter <- 2500L
    warning("Illegal splitter input, splitter will be set to 2500",
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
  
  url <- "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi"
  cfg_file <- "/var/www/html/services/TMHMM-2.0/webface.cf"
  file_list <- ragp::split_fasta(path_in = file_name,
                                 path_out = "tmp_tmhmm_",
                                 num_seq = splitter)
  
  if(grepl("temp_", file_name)){
    unlink(file_name)
  }
  output <- vector("list", length(file_list))
  
  for(k in seq_along(file_list)){
    x <- file_list[k]
    file_up <- httr::upload_file(x)
    
    res <- httr::POST(configfile = cfg_file,
                      url = url,
                      encode = "multipart",
                      body = list(seqfile = file_up,
                                  configfile = cfg_file,
                                  outform = "-short",
                                  SEQ = NULL,
                                  version = NULL))
  
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
      message(paste("batch", k, "jobid is:", jobid))
    }
    
    time1 <- Sys.time()
    
    repeat {
      res2 <- httr::GET(url = url,
                        query = list(jobid = jobid,
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
                    x, ".fa"),
             call. = FALSE)
      }
      res2 <- as.character(
        xml2::xml_find_all(
          httr::content(res2,
                        as = "parsed"),
          ".//pre")
      )
      
      Sys.sleep(1)
      if (length(res2) > 0) {
        if (grepl("pre", res2)){
        break
        }
      }
      
      time2 <- Sys.time()
      
      max.time <- as.difftime(pmax(50, splitter),
                              units = "secs")
      
      if ((time2 - time1) > max.time) {
        res2 <- NULL
        if(progress) message(
          "file",
          x[i],
          "took longer then expected")
        break
      }
    }
    if (is.null(res2)) {
      tms <- 0
      while(tms < attempts && is.null(res2)){
        if(progress) message(
          "reattempting file",
          x)
        file_up <- httr::upload_file(x)
        
        res <- httr::POST(configfile = cfg_file,
                          url = url,
                          encode = "multipart",
                          body = list(seqfile = file_up,
                                      configfile = cfg_file,
                                      outform = "-short",
                                      SEQ = NULL,
                                      version = NULL))
        
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
          message(paste("batch", k, "jobid is:", jobid))
        }
        
        time1 <- Sys.time()
        
        repeat {
          res2 <- httr::GET(url = url,
                            query = list(jobid = jobid,
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
            stop(paste0(prt, ". Problem in file: ",
                        x),
                 call. = FALSE)
          }
          res2 <- as.character(
            xml2::xml_find_all(
              httr::content(res2,
                            as = "parsed"),
              ".//pre")
          )
          
          Sys.sleep(1)
          if (length(res2) > 0) {
            if (grepl("pre", res2)){
              break
            }
          }
          
          time2 <- Sys.time()
          
          max.time <- as.difftime(pmax(100, splitter * 1.5),
                                  units = "secs")
          
          if ((time2 - time1) > max.time) {
            res2 <- NULL
            break
          }
        }
        tms <- tms + 1
      }
    }
    if(is.null(res2)){
      output <- do.call(rbind,
                        output)
      warning(
        "maximum attempts reached at",
        x,
        "returning finished queries",
        call. = FALSE)
      return(output)
    }
    unlink(x)
    
    res2 <- sub("^<pre>", "", res2)
    res2 <- sub("</pre>.*$", "", res2)
    
    res2 <- read.delim(text = res2,
                       sep = "\t",
                       header = FALSE,
                       stringsAsFactors = FALSE)
    
    colnames(res2) <- c("id",
                        "len",
                        "ExpAA",
                        "First60",
                        "PredHel",
                        "Topology")
    
    for(i in 2:6){
      res2[,i] <- sub(
        paste(
          colnames(res2)[i], "[.=]", sep = ""), "", res2[,i])
    }
    
    colnames(res2)[c(2, 5, 6)] <- c("length",
                                    "tm",
                                    "prediction")
    
    for(i in c(2,5)){
      res2[,i] <- as.integer(as.character(res2[,i]))
    }
    
    for(i in 3:4){
      res2[,i] <- as.numeric(as.character(res2[,i]))
    }
    output[[k]] <- res2
  }
  output <- do.call(rbind,
                    output)
  return(output)
}

#' @rdname get_tmhmm
#' @method get_tmhmm data.frame
#' @export

get_tmhmm.data.frame <- function(data,
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
  res <- get_tmhmm.character(data = file_name, ...)
  return(res)
}

#' @rdname get_tmhmm
#' @method get_tmhmm list
#' @export


get_tmhmm.list <- function(data,
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
  res <- get_tmhmm.character(data = file_name, ...)
  return(res)
}

#' @rdname get_tmhmm
#' @method get_tmhmm default
#' @export

get_tmhmm.default <- function(data = NULL,
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
  res <- get_tmhmm.character(data = file_name, ...)
  return(res)
}

#' @rdname get_tmhmm
#' @method get_tmhmm AAStringSet
#' @export

get_tmhmm.AAStringSet <-  function(data,
                                   ...){
  sequence <- as.character(data)
  id <- names(sequence)
  sequence <- unname(sequence)
  sequence <- toupper(sequence)
  sequence <- sub("\\*$",
                  "",
                  sequence)
  
  res <- get_tmhmm.default(sequence = sequence,
                           id = id,
                           ...)
  return(res)
}
