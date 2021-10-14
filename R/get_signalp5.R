#' Query SignalP 5.0 web server.
#'
#' The SignalP 5.0 server predicts the presence of signal peptides and the location of their cleavage sites in proteins from Archaea, Gram-positive Bacteria, Gram-negative Bacteria and Eukarya.
#'
#' @aliases get_signalp5 get_signalp5.default get_signalp5.character get_signalp5.data.frame get_signalp5.list get_signalp5.AAStringSet
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class "SeqFastaAA" resulting from \code{\link[seqinr]{read.fasta}} call. Alternatively an `AAStringSet` object. Should be left blank if vectors are provided to sequence and id arguments.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param org_type One of "euk", "gram-", "gram+" or "archea". Default is "euk". Are the protein sequences from Eukarya, Gram-negative Bacteria, Gram-positive Bacteria or Archaea.
#' @param splitter An integer indicating the number of sequences to be in each .fasta file that is to be sent to the server. Default is 2500. Change only in case of a server side error. Accepted values are in range of 1 to 5000.
#' @param attempts Integer, number of attempts if server unresponsive, at default set to 2.
#' @param progress Boolean, whether to show messages of the job id for each batch. Default is FALSE.
#' @param ... currently no additional arguments are accepted apart the ones documented bellow.
#' 
#' @return  if org_type is "euk" a data frame with columns:
#' \describe{
#'   \item{id}{Character, as from input}
#'   \item{Prediction}{Integer, The type of signal peptide predicted: Sec/SPI, Tat/SPI, Sec/SPII or Other if no signal peptide predicted}
#'   \item{SP.Sec.SPI}{Numeric, marginal probability that the protein contains a Sec  N-terminal signal peptide (Sec/SPI).}
#'   \item{Other}{Numeric, the probability that the sequence does not have any kind of signal peptid.}
#'   \item{CS_pos}{Character, cleavage site position.}
#'   \item{Pr}{Numeric, probability of the predicted cleavage site position.}
#'   \item{cleave.site}{Numeric,llocal amino acid sequence arround the predicted cleavage site.}
#'   \item{is.signalp}{Logical, did SignalP5 predict the presence of a signal peptide.}
#'   \item{sp.length}{Integer, length of the predicted signal peptide.}
#'   }
#'   
#' if org_type is one of "gram-", "gram+" or "archea" the returned data frame will have two additional columns between `SP.Sec.SPI` and `Other`:
#'  
#' \describe{
#'   \item{TAT.Tat.SPI}{Numeric, marginal probability that the protein contains a Tat N-terminal signal peptide (Tat/SPI).}
#'   \item{LIPO.Sec.SPII}{Numeric, marginal probability that the protein contains a Lipoprotein N-terminal signal peptide (Sec/SPII).}
#'   }   
#'   
#'
#' @note This function creates temporary files in the working directory. If something goes wrong during communication with the server and progress was set to TRUE, predictions can be obtained using `file.path("http://www.cbs.dtu.dk/services/SignalP-5.0/tmp", jobid, "output_protein_type.txt")` eg `read.delim(file.path(...), header = TRUE, skip = 1)`.
#'
#'
#' @source \url{http://www.cbs.dtu.dk/services/SignalP/}
#' @references Almagro Armenteros JJ, Tsirigos KD, Sønderby CK, Petersen TN, Winther O, Brunak S, von Heijne G, Nielsen H. (2019) SignalP 5.0 improves signal peptide predictions using deep neural networks. Nature Biotechnology, 37:420-423, doi:10.1038/s41587-019-0036-z
# 
#' @seealso \code{\link[ragp]{get_signalp}} \code{\link[ragp]{get_targetp}}
#'
#' @examples
#' 
#' library(ragp)
#' sp5_pred <- get_signalp5(data = at_nsp[1:10,],
#'                          sequence,
#'                          Transcript.id)
#' sp5_pred
#' 
#' @import seqinr
#' @import httr
#' @import xml2
#' @export

get_signalp5 <- function(data, ...){
  if (missing(data) || is.null(data)) get_signalp5.default(...)
  else UseMethod("get_signalp5")
}

#' @rdname get_signalp5
#' @method get_signalp5 character
#' @export

get_signalp5.character <- function(data,
                                   org_type = c("euk", "gram-", "gram+", "archea"),
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
  
  if (missing(org_type)) {
    org_type <- "euk"
  }
  if (!org_type %in% c("euk", "gram-", "gram+", "archea")) {
    stop("org_type should be one of: 'euk', 'gram-', 'gram+', 'archea'",
         call. = FALSE)
  }
  if (length(org_type) > 1){
    stop("org_type should be one of: 'euk', 'gram-', 'gram+', 'archea'",
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
  
  organism <- c("Eukarya",
                "Gram-negative bacteria",
                "Gram-positive bacteria",
                "Archaea")
  names(organism) <- c("euk", "gram-", "gram+", "archea")
  
  organism <- organism[names(organism) == org_type]
  
  url <- "http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi"
  cfg_file <- "/usr/opt/www/pub/CBS/services/SignalP-5.0/signalp5.cf"
  file_list <- ragp::split_fasta(path_in = file_name,
                                 path_out = "tmp_signalp5_",
                                 num_seq = splitter)
  if(grepl("temp_", file_name)){
    unlink(file_name)
  }


  output <- vector("list", length(file_list))
  
  for(k in seq_along(file_list)){
    x <- file_list[k]
    file_up <- httr::upload_file(x)
    res <-  httr::POST(url = url,
                       encode = "multipart",
                       body = list(configfile = cfg_file,
                                   uploadfile = file_up,
                                   organism = organism,
                                   format = "short"))
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
    
    jobid <- res
    
    if(progress){
      message(paste("batch", k, "jobid is:", jobid))
    }

    time1 <- Sys.time()
    
    repeat {
      res2 <- httr::GET(url = file.path("http://www.cbs.dtu.dk/services/SignalP-5.0/tmp",
                                        jobid,
                                        "output_protein_type.txt"))
      code <- res2$status_code
      
      Sys.sleep(5)
      if (code == 200L) {
        res2 <- read.delim(
          file.path("http://www.cbs.dtu.dk/services/SignalP-5.0/tmp",
                    jobid,
                    "output_protein_type.txt"),
          header = TRUE,
          skip = 1)
        break
      }
    
      time2 <- Sys.time()
      
      max.time <- as.difftime(pmax(50, splitter),
                              units = "secs")
      
      if ((time2 - time1) > max.time) {
        code <- 404L
        if(progress) message("file",
                             x,
                             "took longer then expected")
        break
      }
    }
    if (code == 404L) {
      tms <- 0
      while(tms < attempts && code == 404L){
        if(progress) message(
          "reattempting file",
          x)
        file_up <-  httr::upload_file(x)
        res <-  httr::POST(url = url,
                           encode = "multipart",
                           body = list(configfile = cfg_file,
                                       uploadfile = file_up,
                                       organism = organism,
                                       format = "short"))
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
        jobid <- res
        
        if(progress){
          message(paste("batch", k, "reatempt jobid is:", jobid))
        }

        time1 <- Sys.time() 
        
        repeat {
          res2 <- httr::GET(url = file.path("http://www.cbs.dtu.dk/services/SignalP-5.0/tmp",
                                            jobid,
                                            "output_protein_type.txt"))
          code <- res2$status_code
          Sys.sleep(5)
          if (code == 200L) {
            res2 <- read.delim(
              file.path("http://www.cbs.dtu.dk/services/SignalP-5.0/tmp",
                        jobid,
                        "output_protein_type.txt"),
              header = TRUE,
              skip = 1,
              stringsAsFactors = FALSE)
            break
          }
          time2 <- Sys.time()
          
          max.time <- as.difftime(pmax(100, splitter*1.5),
                                  units = "secs")
          
          if ((time2 - time1) > max.time) {
            code <- 404L
            if(progress) message(
              "file",
              x,
              "took longer then expected")
            break
          }
        }
        tms <- tms + 1
      }
    }
    
    if (code == 404L){
      output <- do.call(rbind,
                        output)
      output$is.signalp <- output$is.sp == "Y"

      warning(
        "maximum attempts reached at",
        x,
          "returning finished queries",
        call. = FALSE)
      return(output)
    }
    unlink(x)
    
    res2 <- data.frame(res2[,colnames(res2) != "CS.Position"],
                       CS_pos = gsub("CS pos: (\\d+)-(\\d+)\\..*$", "\\1-\\2", res2[,"CS.Position"]),
                       Pr = gsub(".*Pr: ", "",  res2[,"CS.Position"]),
                       substr = gsub("^.*?([ARNDCEQGHILKMFPSTWYV]{3}-[ARNDCEQGHILKMFPSTWYV]{2}).*$" ,
                                     "\\1", res2[,"CS.Position"]),
                       stringsAsFactors = FALSE)
    if(org_type == "euk"){
      colnames(res2) <- c("id",
                          "Prediction",
                          "SP.Sec.SPI",
                          "Other",
                          "CS_pos",
                          "Pr",
                          "cleave.site")
      res2$SP.Sec.SPI <- as.numeric(res2$SP.Sec.SPI)
      
    }
    if(org_type != "euk"){
      colnames(res2) <- c("id",
                          "Prediction",
                          "SP.Sec.SPI",
                          "TAT.Tat.SPI",
                          "LIPO.Sec.SPII",
                          "Other",
                          "CS_pos",
                          "Pr",
                          "cleave.site")
      res2$SP.Sec.SPI <- as.numeric(res2$SP.Sec.SPI)
      res2$TAT.Tat.SPI <- as.numeric(res2$TAT.Tat.SPI)
      res2$LIPO.Sec.SPII <- as.numeric(res2$LIPO.Sec.SPII)
    }
    

    res2$Other <- as.numeric(res2$Other)
    res2$Pr <- as.numeric(res2$Pr)
    output[[k]] <- res2
  }
  
  output <- do.call(rbind,
                    output)
  output$is.signalp <- !is.na(output$Pr)
  output$sp.length <- as.integer(gsub("(\\d+)-\\d+", "\\1", output[,"CS_pos"]))
  return(output)
}

#' @rdname get_signalp5
#' @method get_signalp5 data.frame
#' @export
          
get_signalp5.data.frame <- function(data,
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
  res <- get_signalp5.character(data = file_name, ...)
  return(res)
}

#' @rdname get_signalp5
#' @method get_signalp5 list
#' @export


get_signalp5.list <- function(data,
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
  res <- get_signalp5.character(data = file_name, ...)
  return(res)
}

#' @rdname get_signalp5
#' @method get_signalp5 default
#' @export

get_signalp5.default <- function(data = NULL,
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
  res <- get_signalp5.character(data = file_name, ...)
  return(res)
}

#' @rdname get_signalp5
#' @method get_signalp5 AAStringSet
#' @export

get_signalp5.AAStringSet <-  function(data,
                                     ...){
  sequence <- as.character(data)
  id <- names(sequence)
  sequence <- unname(sequence)
  sequence <- toupper(sequence)
  sequence <- sub("\\*$",
                  "",
                  sequence)
  
  res <- get_signalp5.default(sequence = sequence,
                              id = id,
                              ...)
  return(res)
}