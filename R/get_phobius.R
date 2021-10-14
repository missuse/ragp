#' Query Phobius web server.
#'
#' Phobius web server is a combined transmembrane topology and signal peptide (N-sp) predictor. Currently only "normal prediction" of signal peptides is supported by the function.
#'
#' @aliases get_phobius get_phobius.default get_phobius.character get_phobius.data.frame get_phobius.list get_phobius.AAStringSet
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class "SeqFastaAA" resulting from \code{\link[seqinr]{read.fasta}} call. Alternatively an `AAStringSet` object. Should be left blank if vectors are provided to sequence and id arguments.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param progress Boolean, whether to show the progress bar, at default set to FALSE.
#' @param ... currently no additional arguments are accepted apart the ones documented bellow.
#'
#' @return  A data frame with columns:
#' \describe{
#' \item{id}{Character, name of the submitted sequence.}
#' \item{tm}{Integer, the number of predicted transmembrane segments.}
#' \item{sp}{Character, Y/0 indicator if a signal peptide was predicted or not.}
#' \item{prediction}{Character string, predicted topology of the protein.}
#' \item{cut_site}{Integer, first amino acid after removal of the signal peptide}
#' \item{is.phobius}{Logical, did Phobius predict the presence of a signal peptide}
#' \item{sp.length}{Integer, length of the predicted signal peptide.}
#'}
#'
#' @details
#' The topology (prediction column of the result) is given as the position of the transmembrane helices separated by 'i' if the loop is on the cytoplasmic or 'o' if it is on the non-cytoplasmic side. A signal peptide is given by the position of its h-region separated by a n and a c, and the position of the last amino acid in the signal peptide and the first of the mature protein separated by a /.
#'
#' @note This function creates temporary files in the working directory.
#'
#' @source \url{http://phobius.sbc.su.se/}
#'
#' @references Kall O. Krogh A. Sonnhammer E. L. L. (2004) A Combined Transmembrane Topology and Signal Peptide Prediction Method. Journal of Molecular Biology 338(5): 1027-1036
#' @seealso \code{\link[ragp]{get_signalp}}
#'
#' @examples
#' library(ragp)
#' data(at_nsp)
#'
#' phobius_pred <- get_phobius(at_nsp[1:20,],
#'                             sequence,
#'                             Transcript.id)
#' phobius_pred
#' 
#' @import seqinr
#' @import httr
#' @import stringr
#' @import xml2
#' @export get_phobius

get_phobius <- function (data, ...){
  if (missing(data) || is.null(data)) get_phobius.default(...)
  else UseMethod("get_phobius")
}

#' @rdname get_phobius
#' @method get_phobius character
#' @export


get_phobius.character <- function(data,
                                  progress = FALSE,
                                  ...){
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
    }  else {
      stop("cannot find file in the specified path",
           call. = FALSE)
    }
  file_list <- ragp::split_fasta(path_in = file_name,
                                 path_out = "temp_phob_",
                                 num_seq = 500)
  len <- length(file_list)
  if(grepl("temp_", file_name)){
    unlink(file_name)
  }
  if(progress){
    pb <- utils::txtProgressBar(min = 0,
                                max = len,
                                style = 3)
  }
  url <- "https://phobius.sbc.su.se/cgi-bin/predict.pl"
  collected_res = vector("list", len)
  for (i in 1 : len){
    file_up <- httr::upload_file(file_list[i])
    res <- httr::POST(url = url,
                      encode = "multipart",
                      body = list(`protfile` = file_up ,
                                  `format` = "short"))
    res <- httr::content(res,
                         as = "parsed")
    res <- xml2::xml_text(res,
                          "//pre")
    res <- unlist(strsplit(as.character(res),
                           "\n"))
    res <- res[(grep("SEQENCE",
                     res)+1) : (grep("_uacct",
                                     res)-1)]
    res <- strsplit(res,
                    " +")
    res <- do.call(rbind,
                   res)
    res <- as.data.frame(res,
                         stringsAsFactors = FALSE)
    colnames(res) <- c("id",
                       "tm",
                       "sp",
                       "prediction")
    collected_res[[i]] <- res
    unlink(file_list[i])
    if(progress){
      utils::setTxtProgressBar(pb, i)
    }
  }
  if(progress){
    close(pb)
  }
  collected_res <- do.call(rbind,
                           collected_res)
  collected_res$cut_site <- as.integer(
    stringr::str_extract(
      collected_res$prediction,
      "(?<=/)\\d+")
  )
  collected_res$is.phobius <- collected_res$sp == "Y"
  collected_res$sp.length <- collected_res$cut_site - 1
  return(collected_res)
}


#' @rdname get_phobius
#' @method get_phobius data.frame
#' @export

get_phobius.data.frame <- function(data,
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
  res <- get_phobius.character(file_name,
                               ...)
  return(res)
}

#' @rdname get_phobius
#' @method get_phobius list
#' @export

get_phobius.list <- function(data,
                             ...){
  file_name <- paste("temp_",
                     gsub("^X",
                          "",
                          make.names(Sys.time())),
                     ".fasta",
                     sep = "")
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
    seqinr::write.fasta(sequence = strsplit(sequence, ""),
                        name = id,
                        file = file_name)
  } else {
    stop("only lists containing objects of class SeqFastaAA are supported")
  }
  res <- get_phobius.character(file_name,
                               ...)
  return(res)
}

#' @rdname get_phobius
#' @method get_phobius default
#' @export

get_phobius.default <- function(data = NULL,
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
  res <- get_phobius.character(file_name,
                               ...)
  return(res)
}

#' @rdname get_phobius
#' @method get_phobius AAStringSet
#' @export

get_phobius.AAStringSet <-  function(data,
                                     ...){
  sequence <- as.character(data)
  id <- names(sequence)
  sequence <- unname(sequence)
  sequence <- toupper(sequence)
  sequence <- sub("\\*$",
                  "",
                  sequence)
  
  res <- get_phobius.default(sequence = sequence,
                             id = id,
                             ...)
  return(res)
}
