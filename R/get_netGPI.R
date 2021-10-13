#' Query NetGPI - 1.1 web server.
#'
#' NetGPI server offers GPI Anchor predictions
#'
#' @aliases get_netGPI get_netGPI.default get_netGPI.character get_netGPI.data.frame get_netGPI.list
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class "SeqFastaAA" resulting from \code{\link[seqinr]{read.fasta}} call. Should be left blank if vectors are provided to sequence and id arguments.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param splitter An integer indicating the number of sequences to be in each .fasta file that is to be sent to the server. Defaults to 2500. Change only in case of a server side error. Accepted values are in range of 1 to 5000.
#' @param attempts Integer, number of attempts if server unresponsive, at default set to 2.
#' @param progress Boolean, whether to show the progress bar, at default set to FALSE.
#' @param ... currently no additional arguments are accepted apart the ones documented bellow.
#' 
#' @return  A data frame with columns:
#' \describe{
#'   \item{id}{Character, as from input}
#'   \item{length}{Integer, length of the protein sequence}
#'   \item{is.gpi}{Logical, is the protein predicted to be GPI anchored.}
#'   \item{omega_site}{Integer, indicating the sequence position of the omega-site.}
#'   \item{likelihood}{Numeric, likelihood of the prediction.}}
#'
#' @note This function creates temporary files in the working directory.
#'
#' @source \url{https://services.healthtech.dtu.dk/service.php?NetGPI-1.1}
#' @references Gislason MH. Nielsen H. Armenteros JA. AR Johansen AR. (2019) Prediction of GPI-Anchored proteins with pointer neural networks. bioRxiv. doi: https://doi.org/10.1101/838680
# 
#' @seealso \code{\link[ragp]{get_big_pi}} \code{\link[ragp]{get_pred_gpi}}
#'
#' @examples
#' \dontrun{
#' library(ragp)
#' netGPI_pred <- get_netGPI(data = at_nsp[1:10,],
#'                           sequence,
#'                           Transcript.id)
#' netGPI_pred
#' }
#' @import seqinr
#' @import httr
#' @import xml2
#' @export

get_netGPI <- function(data, ...){
  if (missing(data) || is.null(data)) get_netGPI.default(...)
  else UseMethod("get_netGPI")
}

#' @rdname get_netGPI
#' @method get_netGPI character
#' @export

get_netGPI.character <- function(data,
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
  url <- "https://services.healthtech.dtu.dk/cgi-bin/webface2.fcgi"
  cfg_file <- "/var/www/html/services/NetGPI-1.1/webface.cf"
  file_list <- ragp::split_fasta(path_in = file_name,
                                 path_out = "tmp_netGPI_",
                                 num_seq = splitter,
                                 id = TRUE)
  if(grepl("temp_", file_name)){
    unlink(file_name)
  }

  fasta_ids <- file_list$id
  file_list <- file_list$file_list

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
      res <- httr::POST(configfile = cfg_file,
                        url = url,
                        encode = "multipart",
                        body = list(configfile = cfg_file,
                                    uploadfile = file_up,
                                    format = "short"))

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
      jobid[i] <- res
      if(progress){
        utils::setTxtProgressBar(pb,
                                 floor(i/2) + (10 * (k - 1)))
      }
    }
    
    collected_res <- vector("list", length(x))
    for (i in seq_along(x)) {
      time1 <- Sys.time()
      repeat {
        jobidi <- jobid[i]
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
          stop(paste0(prt, ". Problem in file: ", "temp_",
                      i, ".fa"),
               call. = FALSE)
        }
        res2 <- as.character(
          xml2::xml_find_all(
            httr::content(res2,
                          as = "parsed"),
            xpath = "//div[@ng-controller='ResultsCtrl']")
        )
        res2_split <- unlist(
          strsplit(res2,
                   "\n")
        )
        Sys.sleep(1)
        if (any(grepl("Prediction summary", res2_split))) {
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

          res <- httr::POST(configfile = cfg_file,
                            url = url,
                            encode = "multipart",
                            body = list(configfile = cfg_file,
                                        uploadfile = file_up,
                                        format = "short"))

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
                xpath = "//div[@ng-controller='ResultsCtrl']")
            )
            res2_split <- unlist(
              strsplit(res2,
                       "\n")
            )
            Sys.sleep(1)
            if (any(grepl("Prediction summary", res2_split))) {
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
      url_dll <- paste0("https://services.healthtech.dtu.dk/services/NetGPI-1.1/tmp/",
                        jobidi,
                        "/output_protein_type.txt")

      res2_split <- read.table(file = url_dll,
                               header = FALSE,
                               stringsAsFactors = FALSE,
                               sep = "\t")
      res2_split <- res2_split[,1:5]
      
      colnames(res2_split) <- c("id",
                                "length",
                                "is.gpi",
                                "omega_site",
                                "likelihood")
      
      res2_split[,3] <- !grepl("Not",
                               res2_split[,3],
                               fixed = TRUE)
      res2_split[,2] <- as.integer(res2_split[,2])
      
      res2_split[grepl("-",
                       res2_split[,4],
                       fixed = TRUE),4] <- "0"
      res2_split[,4] <- as.integer(res2_split[,4])
      res2_split[res2_split[,4] == 0,4] <- NA_integer_
      res2_split[,5] <- as.numeric(res2_split[,5])
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
  
  if(all(fasta_ids %in% output$id)){
    output <- merge(data.frame(id = fasta_ids,
                               stringsAsFactors = FALSE),
                    output,
                    all.x = TRUE,
                    all.y = TRUE,
                    by = "id",
                    sort = FALSE)
    return(output)
  } else {
    warning("Server changed sequence id's because they contained special characters, returning servers output")
    return(output)
  }
}

#' @rdname get_netGPI
#' @method get_netGPI data.frame
#' @export

get_netGPI.data.frame <- function(data,
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
  res <- get_netGPI.character(data = file_name, ...)
  return(res)
}

#' @rdname get_netGPI
#' @method get_netGPI list
#' @export


get_netGPI.list <- function(data,
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
  res <- get_netGPI.character(data = file_name, ...)
  return(res)
}

#' @rdname get_netGPI
#' @method get_netGPI default
#' @export

get_netGPI.default <- function(data = NULL,
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
  res <- get_netGPI.character(data = file_name, ...)
  return(res)
}
