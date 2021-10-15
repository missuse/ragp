#' Query the Conserved Domain Database.
#'
#' Conserved Domain Database is a protein annotation resource that consists of a collection of well-annotated multiple sequence alignment models for ancient domains and full-length proteins.
#'
#' @aliases get_cdd get_cdd.default get_cdd.character get_cdd.data.frame get_cdd.list get_cdd.AAStringSet
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class \code{\link[seqinr]{SeqFastaAA}} resulting from \code{\link[seqinr]{read.fasta}} call. Alternatively an \code{\link[Biostrings]{AAStringSet}} object. Should be left blank if vectors are provided to sequence and id arguments.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param output  One of "short", "standard" or "full". Default is "short". "Short" returns the highest scoring hit, for each region of the query sequence; "standard" returns the best-scoring hit from each source database, for each region of the query sequence; "full" returns the complete set of hits.
#' @param evalue Numeric, all sequences with E-value lower or equal to this value will be retained in the function output. Used to filter out low similarity matches. If set some queried sequences might be discarded from the output. Suggested values: 1e-2 - 1e-5.
#' @param bitscore Numeric, all sequences with bitscore greater or equal to this value will be retained in the function output. Used to filter out low similarity. If set some queried sequences might be discarded from the output. Suggested values: 10 - 20.
#' @param splitter An integer indicating the number of sequences to be in each .fasta file that is to be sent to the server. Default is 100. Change only in case of a server side error. Accepted values are in range of 1 to 4000.
#' @param progress Boolean, whether to show the progress bar, at default set to FALSE.
#' @param ... currently no additional arguments are accepted apart the ones documented bellow.
#'  
#' @return A data frame with columns: 
#' \describe{
#'   \item{id}{Character, as from input}
#'   \item{hit.type}{Character, CD-Search results can include hit types that represent various confidence levels (specific hits, non-specific hits) and domain model scope (superfamilies, multi-domains).}
#'   \item{pssm.id}{Character, unique identifier for a domain models position-specific scoring matrix (PSSM).}
#'   \item{align_start}{Integer, start of the alignement in the query protein sequence.}
#'   \item{align_end}{Integer, end of the alignement in the query protein sequence.}
#'   \item{evalue}{Numeric, the expect value, or E-value, indicates the statistical significance of the hit as the likelihood the hit was found by chance.}
#'   \item{bitscore}{Numeric, this values is derived from the raw alignment score in which the statistical properties of the scoring system used have been taken into account.}
#'   \item{acc}{Character, the accession number of the hit, which can either be a domain model or a superfamily cluster.}
#'   \item{desc}{Character, the short name of a conserved domain, which concisely defines the domain.}  
#'   \item{incomplete}{Character, if the hit to a conserved domain is partial (i.e., if the alignment found by RPS-BLAST omitted more than 20 percent of the CDs extent at either the n- or c-terminus or both), this column will be populated with one of the following values: N - incomplete at the N-terminus, C - incomplete at the C-terminus, NC - incomplete at both the N-terminus and C-terminus. If the hit to a conserved domain is complete, then this column will be populated with a dash (-).} 
#'   \item{superfamily}{Character, the short name of a conserved domain, which concisely defines the domain.}
#'   \item{definition}{Character, this column is populated only for domain models that are specific or non-specific hits, and it lists the accession number of the superfamily to which the domain model belongs.}   
#'   }
#'   
#'
#' @note This function creates temporary files in the working directory.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi}
#' @references Lu S, Wang J, Chitsaz F, Derbyshire MK, Geer RC, Gonzales NR, Gwadz M, Hurwitz DI, Marchler GH, Song JS, Thanki N, Yamashita RA, Yang M, Zhang D, Zheng C, Lanczycki CJ, Marchler-Bauer A (2020) CDD/SPARCLE: the conserved domain database in 2020, Nucleic Acids Res. 48(D1):D265-D268. doi: https://doi.org/10.1093/nar/gkz991
# 
#' @seealso \code{\link[ragp]{get_hmm}} \code{\link[ragp]{plot_prot}}
#'
#' @examples
#' 
#' library(ragp)
#' cdd_pred <- get_cdd(data = at_nsp[1:10,],
#'                     sequence,
#'                     Transcript.id)
#'                     
#' cdd_pred[,1:11]
#' 
#' @import seqinr
#' @import httr
#' @import xml2
#' @export

get_cdd <- function(data, ...){
  if (missing(data) || is.null(data)) get_cdd.default(...)
  else UseMethod("get_cdd")
}

#' @rdname get_cdd
#' @method get_cdd character
#' @export

get_cdd.character <- function(data,
                              output = c("short", "standard", "full"),
                              splitter = 100L,
                              progress = FALSE,
                              evalue = NULL,
                              bitscore = NULL,
                              ...){
  if (missing(splitter)) {
    splitter <- 100L
  }
  if (length(splitter) > 1){
    splitter <- 100L
    warning("splitter should be of length 1, setting to default: splitter = 100",
            call. = FALSE)
  }
  if (!is.numeric(splitter)){
    splitter <- as.numeric(splitter)
    warning("splitter is not numeric, converting using 'as.numeric'",
            call. = FALSE)
  }
  if (is.na(splitter)){
    splitter <- 100L
    warning("splitter was set to NA, setting to default: splitter = 100",
            call. = FALSE)
  }
  if (is.numeric(splitter)) {
    splitter <- floor(splitter)
  }
  if (!(splitter %in% 1:4000)) {
    splitter <- 100L
    warning("Illegal splitter input, splitter will be set to 100",
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
  if (!missing(evalue)) {
    if (length(evalue) > 1){
      stop("evalue should be of length 1",
           call. = FALSE)
    }
    if (!is.numeric(evalue)){
      stop("evalue should be numeric of length 1",
           call. = FALSE)
    }
    if (is.nan(evalue)){
      stop("evalue should be numeric of length 1",
           call. = FALSE)
    }
    if (is.na(evalue)){
      stop("evalue should be numeric of length 1",
           call. = FALSE)
    }
  }
  if (!missing(bitscore)) {
    if (length(bitscore) > 1){
      stop("bitscore should be of length 1",
           call. = FALSE)
    }
    if (!is.numeric(bitscore)){
      stop("bitscore should be numeric of length 1",
           call. = FALSE)
    }
    if (is.nan(bitscore)){
      stop("bitscore should be numeric of length 1",
           call. = FALSE)
    }
    if (is.na(bitscore)){
      stop("bitscore should be numeric of length 1",
           call. = FALSE)
    }
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
  
  output <- match.arg(output)

  dmode <- c("rep", "std", "full")
  names(dmode) <- c("short", "standard", "full")
  dmode <- dmode[names(dmode) == output]

  file_list <- ragp::split_fasta(path_in = file_name,
                                 path_out = "tmp_get_cdd_",
                                 num_seq = splitter)
  if(grepl("temp_", file_name)){
    unlink(file_name)
  }
  for_pb <- length(file_list)
  if(progress){
    pb <- utils::txtProgressBar(min = 0,
                                max = for_pb,
                                style = 3)
  }
  url <- "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?"
  
  out <- vector("list", length(file_list))
  
  for(k in seq_along(file_list)){
    x <- file_list[k]
    file_up <- httr::upload_file(x)
    
    res <- httr::POST(
      url = url,
      encode = "multipart",
      body = list(db = "cdd",
                  smode = "auto",
                  compbasedadj = "1",
                  filter = "false",
                  fqueries = file_up))
    
    xml2::xml_text(
      xml2::xml_find_all(
        httr::content(res,
                      as = "parsed",
                      encoding = "UTF-8"),
        ".//td[@id='id_cdsid']")
    ) -> cdd_id
    
    res2 <- httr::POST(url = url,
                       body = list(tdata = "hits",
                                   cdsid = cdd_id,
                                   dmode = dmode,
                                   cddefl = "true"))
    
    tb <- read.delim(text = httr::content(res2,
                                          as = "parsed",
                                          encoding = "UTF-8"),
                     header = FALSE,
                     stringsAsFactors = FALSE)
    
    cdsid <- tb$V2[grepl("cdsid", tb$V1)]
    
    status <- "3"
    while(status == "3"){
      res2 <- httr::POST(url = url,
                         body = list(tdata = "hits",
                                     cdsid = cdsid,
                                     dmode = dmode,
                                     cddefl = "true"))
      status <- gsub("^.*status\\t(\\d).*$", "\\1",
                     httr::content(res2,
                                   as = "text",
                                   encoding = "UTF-8"))
      Sys.sleep(3)
    }
    
    if(status != "0"){
      stop("something went wrong on the server side")
    }
    
    res2 <- httr::POST(url = url,
                       body = list(tdata = "hits",
                                   cdsid = cdsid,
                                   dmode = dmode,
                                   cddefl = "true"))
    
    res2 <-  httr::content(res2,
                           as = "text",
                           encoding = "UTF-8")
    
    res2 <- sub("^.*#status\\tsuccess", "", res2)
    
    res2 <- read.delim(text = res2,
                       stringsAsFactors = FALSE)

    out[[k]] <- res2
    unlink(x)
  }
  if(progress){
    utils::setTxtProgressBar(pb,
                             for_pb)
    close(pb)
  }
  
  out <- do.call(rbind,
                 out)
  
  colnames(out)[1] <- "id"
  colnames(out)[colnames(out) == "Accession"] <- "acc"
  colnames(out)[colnames(out) == "Short.name"] <- "desc"
  colnames(out)[colnames(out) == "From"] <- "align_start"
  colnames(out)[colnames(out) == "To"] <- "align_end"
  colnames(out)[colnames(out) == "E.Value"] <- "evalue"
  colnames(out)[colnames(out) == "Bitscore"] <- "bitscore"
  
  colnames(out) <- tolower(colnames(out))
  
  out$id <- sub("^Q#\\d+ - >", "", out$id)
  
  out$incomplete <- trimws(out$incomplete)
  out$superfamily <- trimws(out$superfamily)
  
  out$superfamily <- ifelse(out$superfamily == "-",
                            NA_character_,
                            out$superfamily)
  
  out$align_start <- as.integer(trimws(out$align_start))
  out$align_end <- as.integer(trimws(out$align_end))
  out$evalue <- as.numeric(trimws(out$evalue))
  out$bitscore <- as.numeric(trimws(out$bitscore))
  if (!missing(evalue)) {
    out <- out[!is.na(out$evalue) & out$evalue <= evalue,]
  }
  if (!missing(bitscore)) {
    out <- out[!is.na(out$bitscore) & out$bitscore >= bitscore,]
  }
  
  return(out)
}

#' @rdname get_cdd
#' @method get_cdd data.frame
#' @export

get_cdd.data.frame <- function(data,
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
  res <- get_cdd.character(data = file_name, ...)
  return(res)
}

#' @rdname get_cdd
#' @method get_cdd list
#' @export


get_cdd.list <- function(data,
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
  res <- get_cdd.character(data = file_name, ...)
  return(res)
}

#' @rdname get_cdd
#' @method get_cdd default
#' @export

get_cdd.default <- function(data = NULL,
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
  res <- get_cdd.character(data = file_name, ...)
  return(res)
}

#' @rdname get_cdd 
#' @method get_cdd AAStringSet
#' @export

get_cdd.AAStringSet <-  function(data,
                                 ...){
  sequence <- as.character(data)
  id <- names(sequence)
  sequence <- unname(sequence)
  sequence <- toupper(sequence)
  sequence <- sub("\\*$",
                  "",
                  sequence)
  
  res <- get_cdd.default(sequence = sequence,
                         id = id,
                         ...)
  return(res)
}

