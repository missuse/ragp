#' Scraping PredGPI web server.
#'
#' PredGPI web server is a predictor of GPI modification sites.
#'
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
#' @export


get_pred_gpi <- function(data = NULL,
                         sequence,
                         id,
                         spec = 0.99){
  tmr <- paste("temp_",
               gsub("^X",
                    "",
                    make.names(Sys.time())),
               ".fasta",
               sep = "")
  if (!is.numeric(as.numeric(spec))){
    spec <- 0.99
    warning("spec could not be converted to numeric, setting to default: spec = 0.99")
  }
  if (is.na(spec)){
    spec <- 0.99
    warning("spec was set to NA, setting to default: spec = 0.99")
  }
  if(length(spec) > 1){
    spec <- spec[1]
    warning("spec has more than one element, using spec[1]")
  }
  if (as.numeric(spec) > 1) {
    spec <- 0.99
    warning("spec must take values in the range 0 - 1,
            it was set to the default: spec = 0.99")
  }
  if (as.numeric(spec) < 0) {
    spec <- 0.99
    warning("spec must take values in the range 0 - 1,
            it was set to the default: spec = 0.99")
  }    
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
                                path_out = "temp_predgpi_",
                                num_seq = 500)
  len = length(file_list)
  if(class(data) != "character"){
    if(file_name == tmr){
      unlink(file_name)
    }
  }
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




