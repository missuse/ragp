#' Query PredGPI web server.
#'
#' PredGPI web server is a predictor of GPI modification sites.
#'
#' @aliases get_pred_gpi get_pred_gpi.default get_pred_gpi.character get_pred_gpi.data.frame get_pred_gpi.list
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class "SeqFastaAA" resulting from seqinr::read.fasta call.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param spec Numeric in the 0-1 range, indicating the threshold specificity.
#'
#' @return  A data frame with columns:
#' \describe{
#' \item{id}{Character, name of the submitted sequence.}
#' \item{omega_site}{Integer, indicating the sequence position of the omega-site.}
#' \item{specificity}{Numeric, indicating the specificity of the prediction.}
#' \item{is.gpi}{Logical, is the specificity of the prediction higher than the set threshold specificity}
#'}
#'
#' @note This function creates temporary files in the working directory.
#'
#' @source \url{http://gpcr.biocomp.unibo.it/predgpi/pred.htm}
#'
#' @references Pierleoni A. Martelli P. L. Casadio R. (2008) PredGPI: a GPI-anchor predictor. BMC Bioinformatics 9: 392
#' @seealso \code{\link[ragp]{get_big_pi}}
#'
#' @examples
#' library(ragp)
#' data(at_nsp)
#'
#' gpi_pred <- get_pred_gpi(at_nsp[1:20,],
#'                          sequence,
#'                          Transcript.id)
#'
#' @import seqinr
#' @import httr
#' @import xml2
#' @export

get_pred_gpi <- function (data, ...){
  if (missing(data) || is.null(data)) get_pred_gpi.default(...)
  else UseMethod("get_pred_gpi")
}

#' @rdname get_pred_gpi
#' @method get_pred_gpi character
#' @export

get_pred_gpi.character <- function(data,
                                   spec = 0.99){
  tmr <- paste("temp_",
               gsub("^X",
                    "",
                    make.names(Sys.time())),
               ".fasta",
               sep = "")
  if (!is.numeric(as.numeric(spec))){
    spec <- 0.99
    warning("spec could not be converted to numeric, setting to default: spec = 0.99",
            call. = FALSE)
  }
  if (is.na(spec)){
    spec <- 0.99
    warning("spec was set to NA, setting to default: spec = 0.99",
            call. = FALSE)
  }
  if(length(spec) > 1){
    spec <- spec[1]
    warning("spec has more than one element, using spec[1]",
            call. = FALSE)
  }
  if (as.numeric(spec) > 1) {
    spec <- 0.99
    warning("spec must take values in the range 0 - 1,
            it was set to the default: spec = 0.99",
            call. = FALSE)
  }
  if (as.numeric(spec) < 0) {
    spec <- 0.99
    warning("spec must take values in the range 0 - 1,
            it was set to the default: spec = 0.99",
            call. = FALSE)
  }
  if(class(data) == "character"){
    if (file.exists(data)){
      file_name <- data
    } else {
      stop("cannot find file in the specified path",
           call. = FALSE)
    }
  }
  file_list <- ragp::split_fasta(path_in = file_name,
                                 path_out = "temp_predgpi_",
                                 num_seq = 500,
                                 id = TRUE)
  if(grepl("temp_", file_name)){
    unlink(file_name)
  }
  id <- file_list$id
  file_list <- file_list$file_list
  len <- length(file_list)
  pb <- utils::txtProgressBar(min = 0,
                              max = len,
                              style = 3)
  url <- "http://gpcr.biocomp.unibo.it/cgi-bin/predictors/gpi/gpipe_1.4.cgi"
  collected_res = vector("list", len)
  for (i in 1 : len){
    file_up <- httr::upload_file(file_list[i])
    res <- httr::POST(url = url,
                      encode = "multipart",
                      body = list(tipo_hmm = "0",
                                  upfile = file_up,
                                  Display = "Display")
    )
    res <- httr::content(res,
                         as = "parsed")
    res <- as.character(
      xml2::xml_children(
        xml2::xml_find_all(res,
                           xpath = ".//table")[2])
    )
    res <- lapply(res, function(x){
      z <- gsub("<.*?>", "", x)
      z <- unlist(strsplit(z, "\n"))
      z <- z[z != ""]
      z <- trimws(z)
      z
    }
    )
    res <- do.call(rbind, res)
    res <- as.data.frame(res,
                         stringsAsFactors = FALSE)

    colnames(res) <- c("id",
                       "omega_site",
                       "specificity",
                       "seq")
    row1 <- grep("omega-site position",
                 res[,2])
    res <- res[(row1+1):nrow(res),]
    res <- res[,-4]
    collected_res[[i]] <- res
    unlink(file_list[i])
    Sys.sleep(1)
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  collected_res <- do.call(rbind,
                           collected_res)
  collected_res <- merge(data.frame(id),
                         collected_res,
                         all.x = TRUE,
                         sort = FALSE)
  collected_res
  collected_res$id <- as.character(collected_res$id)
  collected_res$omega_site <- as.integer(collected_res$omega_site)
  collected_res$specificity <- sub("%$", "", collected_res$specificity)
  collected_res$specificity <- as.numeric(collected_res$specificity)/100
  collected_res$is.gpi <- collected_res$specificity >= spec
  collected_res
}


#' @rdname get_pred_gpi
#' @method get_pred_gpi data.frame
#' @export

get_pred_gpi.data.frame <- function(data,
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
  res <- get_pred_gpi.character(file_name, ...)
  return(res)
}

#' @rdname get_pred_gpi
#' @method get_pred_gpi list
#' @export

get_pred_gpi.list <- function(data, ...){
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
  res <- get_pred_gpi.character(data = file_name, ...)
  return(res)
}

#' @rdname get_pred_gpi
#' @method get_pred_gpi default
#' @export

get_pred_gpi.default <- function(sequence, id, ...){
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
  res <- get_pred_gpi.character(data = file_name, ...)
  return(res)
}


