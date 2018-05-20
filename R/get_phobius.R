#' Scraping Phobius web server.
#'
#' Phobius web server is a combined transmembrane topology and signal peptide (N-sp) predictor. Currently only "normal prediction" of signal peptides is supported by the function.
#'
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class "SeqFastaAA" resulting from seqinr::read.fasta call.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#'
#' @return  A data frame with columns:
#' \describe{
#' \item{Name}{Character, name of the submitted sequence.}
#' \item{tm}{Integer, the number of predicted transmembrane segments.}
#' \item{sp}{Character, Y/0 indicator if a signal peptide was predicted or not.}
#' \item{prediction}{Character string, predicted topology of the protein.}
#' \item{cut_site}{Integer, first amino acid after removal of the signal peptide}
#' \item{is.phobius}{Logical, did Phobius predict the presence of a signal peptide}
#'}
#'
#' @details
#' The topology (prediction column of the result) is given as the position of the transmembrane helices separated by 'i' if the loop is on the cytoplasmic or 'o' if it is on the non-cytoplasmic side. A signal peptide is given by the position of its h-region separated by a n and a c, and the position of the last amino acid in the signal peptide and the first of the mature protein separated by a /.
#'
#' @source \url{http://phobius.sbc.su.se/}
#'
#' @references Kall O. Krogh A. Sonnhammer E. L. L. (2004) A Combined Transmembrane Topology and Signal Peptide Prediction Method. Journal of Molecular Biology 338(5): 1027-1036
#' @seealso \code{\link[ragp]{get_phobius_file}}
#'
#' @examples
#' library(ragp)
#' data(at_nsp)
#' 
#' phobius_pred <- get_phobius(at_nsp[1:20,],
#'                             sequence,
#'                             Transcript.id)
#' @export

get_phobius <- function(data = NULL, sequence, id){
  tmr <- paste("temp_",
               gsub("^X",
                    "",
                    make.names(Sys.time())),
               ".fasta",
               sep = "")
  if(missing(data)){
    if (missing(sequence)){
      stop("protein sequence must be provided to obtain predictions")
    }
    if (missing(id)){
      stop("protein id must be provided to obtain predictions")
    }
    id <- as.character(id)
    sequence <- toupper(as.character(sequence))
    if (length(sequence) != length(id)){
      stop("id and sequence vectors are not of same length")
    }
    sequence <- sub("\\*$", "", sequence)
    file_name <- tmr
    seqinr::write.fasta(sequence = strsplit(sequence, ""),
                        name = id, file = file_name)
  }
  if(class(data[[1]]) ==  "SeqFastaAA"){
    dat <- lapply(data, paste0, collapse ="")
    id <- names(dat)
    sequence <- toupper(as.character(unlist(dat)))
    sequence <- sub("\\*$", "", sequence)
    file_name <- tmr
    seqinr::write.fasta(sequence = strsplit(sequence, ""),
                        name = id, file = file_name)
  }
  if(class(data) == "data.frame"){
    if(missing(sequence)){
      stop("the column name with the sequences must be specified")
    }
    if(missing(id)){
      stop("the column name with the sequence id's must be specified")
    }
    id <- as.character(substitute(id))
    sequence <- as.character(substitute(sequence))
    if (length(id) != 1L){
      stop("only one column name for 'id' must be specifed")
    }
    if (length(sequence) != 1L){
      stop("only one column name for 'sequence' must be specifed")
    }
    id <- if(id %in% colnames(data)){
      data[[id]]
    } else {
      stop("specified 'id' not found in data")
    }
    id <- as.character(id)  
    sequence  <- if(sequence %in% colnames(data)){
      data[[sequence]]
    } else {
      stop("specified 'sequence' not found in data")
    }
    sequence <- toupper(as.character(sequence))
    sequence <- sub("\\*$", "", sequence)
    file_name <- tmr
    seqinr::write.fasta(sequence = strsplit(sequence, ""),
                        name = id, file = file_name)
  }
  if(class(data) == "character"){
    if (file.exists(data)){
      file_name <- data
    } else {
      stop("cannot find file in the specified path")
    }
  }
  file_list = ragp::split_fasta(path_in = file_name,
                                path_out = "temp_phob_",
                                num_seq = 500)
  len = length(file_list)
  url <- "http://phobius.binf.ku.dk/cgi-bin/predict.pl"
  collected_res = vector("list", len)
  for (i in 1 : len){
    file_up <- httr::upload_file(file_list[i])
    res <- httr::POST(url = url,
                      encode = "multipart",
                      body = list(`protfile` = file_up ,
                                  `format` = "short"))
    res <- httr::content(res, as = "parsed")
    res <- xml2::xml_text(res, "//pre")
    res <- unlist(strsplit(as.character(res), "\n"))
    res <- res[(grep("SEQENCE", res)+1) : (grep("_uacct", res)-1)]
    res <- strsplit(res, " +")
    res <- do.call(rbind, res)
    res <- as.data.frame(res, stringsAsFactors = FALSE)
    colnames(res) <- c("Name", "tm", "sp", "prediction")
    collected_res[[i]] <- res
    unlink(file_list[i])
  }
  if(class(data) != "character"){
    if(file_name == tmr){
      unlink(file_name)
    }
  }
  collected_res <- do.call(rbind, collected_res)
  collected_res$cut_site <- stringr::str_extract(collected_res$prediction,
                                                 "(?<=/)\\d+")
  collected_res$is.phobius <- collected_res$sp == "Y"
  return(collected_res)
}
