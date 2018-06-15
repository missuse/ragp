#' Scraping hmmer web server.
#'
#' hmmer web server offers biosequence analysis using profile hidden Markov Models. This function allows searching
#' of a protein sequence vs a profile-HMM database (Pfam-A).
#'
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class "SeqFastaAA" resulting from seqinr::read.fasta call.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param verbose Bolean, whether to print out the output for each sequence, defaults to TRUE.
#' @param sleep Numeric indicating the pause in seconds between server calls, at default set to 1.
#' @param attempts Integer, number of attempts if server unresponsive, at default set to 2.
#' @param timeout Numeric, time in seconds to wait for server response.
#' @param progress Bolean, whether to show the progress bar, at default set to TRUE
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{id}{Character, as supplied in the function call}
#'   \item{name}{Character, PFAM family name}
#'   \item{acc}{Character, PFAM family accession}
#'   \item{desc}{Character, PFAM family description}
#'   \item{clan}{Character, PFAM clan}
#'   \item{align_start}{Numeric, start of domain alignment in query sequence}
#'   \item{align_end}{Numeric, end of domain alignment in query sequence}
#'   \item{model_start}{Numeric, start of alignment in domain model}
#'   \item{model_end}{Numeric, end of alignment in domain model}
#'   \item{ievalue}{Numeric, the "independent E-value", the E-value that the sequence/profile comparison would have received if this were the only domain envelope found in it, excluding any others. This is a stringent measure of how reliable this particular domain may be. The independent E-value uses the total number of targets in the target database.}
#'   \item{cevalue}{Numeric, the "conditional E-value", a permissive measure of how reliable this particular domain may be.}
#'   \item{bitscore}{Numeric, the domain bit score.}
#'   \item{reported}{Logical, is the result reported on the hmmer site. The hmmer web server outputs more hmm profile matches than it presents to the user. Results below a certain threshold are not reported (hidden) on the site.}
#' }
#'
#' @note hmmscan does not handle sequences longer than 1000 amino acids. get_hmm splits these sequences into shorter substrings which overlap by 300 amino acids and queries hmmscan. Some results might be redundant or partially overlapping in this case. When this is an issue it is advisable to provide a subsequence of appropriate length as get_hmm input.
#'
#' @source \url{https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan}
#'
#' @seealso \code{\link[ragp]{pfam2go}}
#'
#' @examples
#'
#' pfam_pred <- get_hmm(data = at_nsp[1:5,],
#'                     sequence = sequence,
#'                     id = Transcript.id,
#'                     verbose = FALSE,
#'                     sleep = 0)
#'
#' @export


get_hmm <- function(data = NULL,
                    sequence,
                    id,
                    verbose = TRUE,
                    sleep = 1,
                    attempts = 2L,
                    timeout = 10,
                    progress = TRUE){
  if (missing(verbose)) {
    verbose <- TRUE
  }
  if (length(verbose) > 1){
    verbose <- TRUE
    warning("verbose should be of length 1, setting to default: verbose = TRUE")
  }
  if (!is.logical(verbose)){
    verbose <- as.logical(verbose)
    warning("verbose is not logical, converting using 'as.logical'")
  }
  if (is.na(verbose)){
    verbose <- TRUE
    warning("verbose was set to NA, setting to default: verbose = TRUE")
  }
  if (missing(progress)) {
    progress <- TRUE
  }
  if (length(progress) > 1){
    progress <- TRUE
    warning("progress should be of length 1, setting to default: progress = TRUE")
  }
  if (!is.logical(progress)){
    progress <- as.logical(progress)
    warning("progress is not logical, converting using 'as.logical'")
  }
  if (is.na(progress)){
    progress <- TRUE
    warning("progress was set to NA, setting to default: progress = TRUE")
  }
  if (missing(sleep)) {
    sleep <- 1
  }
  if (length(sleep) > 1){
    sleep <- 1
    warning("sleep should be of length 1, setting to default: sleep = 1")
  }
  if (!is.numeric(sleep)){
    sleep <- as.numeric(sleep)
    warning("sleep is not numeric, converting using 'as.numeric'")
  }
  if (is.na(sleep)){
    sleep <- 1
    warning("sleep was set to NA, setting to default: sleep = 1")
  }
  if (length(attempts) > 1){
    attempts <- 2L
    warning("attempts should be of length 1, setting to default: attempts = 2")
  }
  if (!is.numeric(attempts)){
    attempts<- as.numeric(attempts)
    warning("attempts is not numeric, converting using 'as.numeric'")
  }
  if (is.na(attempts)){
    attempts <- 2L
    warning("attempts was set to NA, setting to default: attempts = 2")
  }
  if (is.numeric(attempts)) {
    attempts <- floor(attempts)
  }
  if (attempts < 1) {
    attempts <- 2L
    warning("attempts was set to less then 1, setting to default: attempts = 2")
  }
  if (missing(timeout)) {
    timeout <- 10
  }
  if (length(timeout) > 1){
    timeout <- 10
    warning("timeout should be of length 1, setting to default: timeout = 10")
  }
  if (!is.numeric(timeout)){
    timeout <- as.numeric(timeout)
    warning("timeout is not numeric, converting using 'as.numeric'")
  }
  if (is.na(timeout)){
    timeout <- 10
    warning("timeout was set to NA, setting to default: timeout = 10")
  }
  if (timeout < 2){
    timeout <- 10
    warning("timeout was set to less then 2, setting to default: timeout = 10")
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
  }
  if(class(data[[1]]) ==  "SeqFastaAA"){
    dat <- lapply(data, paste0, collapse ="")
    id <- names(dat)
    sequence <- toupper(as.character(unlist(dat)))
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
  }
  if(class(data) == "character"){
    if (file.exists(data)){
      dat <- seqinr::read.fasta(file = data,
                                seqtype = "AA",
                                as.string = FALSE)
      dat <- lapply(dat, paste0, collapse ="")
      id <- names(dat)
      sequence <- toupper(as.character(unlist(dat)))
    } else {
      stop("cannot find file in the specified path")
    }
  }
  sequence <- sub("\\*$", "", sequence)
  n <- length(sequence)
  url <- "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan"
  pfam <- list()
  if(progress){
    pb <- utils::txtProgressBar(min = 0,
                                max = n,
                                style = 3)
  }
  for (i in 1:n) {
    seqi <- sequence[i]
    seqi_id <- id[i]
    nc <- nchar(seqi)
    if(nc > 1000){
      message(cat("\n",
                  "sequence",
                  seqi_id,
                  "is longer than the permitted hmmscan length",
                  "\n",
                  "processing in chunks"))
      star <- c(1, seq(from = 701,
                       to = nc,
                       by = 700))
      end <- unique(c(seq(from = 1000,
                          to = nc,
                          by = 700),
                      nc))
      if(length(star) > length(end)){
        star <- star[-length(star)]
      }
      seqas <- data.frame(seqas = substring(seqi,
                                            first = star,
                                            last = end),
                          star = star,
                          end = end,
                          id = paste(seqi_id,
                                     star,
                                     sep = "_xx_"))
      resis <- ragp::get_hmm(sequence = seqas$seqas,
                             id = seqas$id,
                             progress = FALSE,
                             verbose = FALSE,
                             attempts = attempts,
                             timeout = timeout,
                             sleep = sleep)
      st <- unlist(lapply(strsplit(resis$id,
                                   "_xx_"),
                          function(x)
                            as.numeric(x[2])))-1
      resis$align_start <- resis$align_start + st
      resis$align_end <- resis$align_end + st
      resis$id <- seqi_id
      resis <- resis[!duplicated(resis[,-11]),]
      if (verbose == T) {
        print(i, quote = FALSE)
        print(resis[, 1:11])
        utils::flush.console()
      }
      pfam[[i]] <- resis
    } else {
      res <- NULL
      attempt <- 0
      while (is.null(res) && attempt <= attempts) {
        attempt <- attempt + 1
        try(res <- httr::POST(url = url,
                              encode = "form", 
                              body = list(isajax = "1",
                                          seq = seqi,
                                          hmmdb = "pfam", 
                                          alt = "Paste",
                                          threshold = "cut_ga"),
                              httr::timeout(timeout)), 
            silent = TRUE)
        if (!is.null(res)) {
          break
        }
      }
      if (i == 1) {
        if (attempt > attempts) {
          stop("server was unresponsive, consdier increasing timeout and attempts arguments")
        }
      }
      else {
        if (attempt > attempts) {
          pfam <- suppressWarnings(do.call(rbind, pfam))
          pfam <- as.data.frame(pfam)
          rownames(pfam) <- 1:nrow(pfam)
          pfam$id <- as.character(pfam$id)
          pfam$name <- as.character(pfam$name)
          pfam$acc <- as.character(pfam$acc)
          pfam$desc <- as.character(pfam$desc)
          pfam$clan <- as.character(pfam$clan)
          pfam$align_start <- as.numeric(as.character(pfam$align_start))
          pfam$align_end <- as.numeric(as.character(pfam$align_end))
          pfam$model_start <- as.numeric(as.character(pfam$model_start))
          pfam$model_end <- as.numeric(as.character(pfam$model_end))
          pfam$ievalue <- as.numeric(as.character(pfam$ievalue))
          pfam$cevalue <- as.numeric(as.character(pfam$cevalue))
          pfam$bitscore <- as.numeric(as.character(pfam$bitscore))
          pfam$reported <- as.logical(as.integer(pfam$reported))
          warning(paste("maximum attempts reached at", 
                        id[i], "returning finished queries"))
          return(pfam)
        }
      }
      res <- httr::content(res, as = "parsed")
      hits <- res$results$hits
      if (length(hits) != 0) {
        hitx <- lapply(hits, function(y) {
          dom <- y$domains
          doms <- lapply(dom, function(x) {
            resi <- data.frame(id = seqi_id,
                               name = x$alihmmname,
                               acc = x$alihmmacc, 
                               desc = x$alihmmdesc,
                               clan = ifelse(is.null(x$clan),
                                             NA,
                                             as.character(x$clan)),
                               align_start = x$alisqfrom, 
                               align_end = x$alisqto,
                               model_start = x$alihmmfrom, 
                               model_end = x$alihmmto,
                               ievalue = x$ievalue, 
                               cevalue = x$cevalue,
                               bitscore = x$bitscore, 
                               reported = as.logical(
                                 as.integer(
                                   x$is_included)),
                               stringsAsFactors = F)
            return(resi)
          })
          domx <- data.frame(do.call(rbind, doms))
          return(domx)
        })
        hitx <- do.call(rbind, hitx)
        hitx <- as.data.frame(hitx,
                              stringsAsFactors = F)
        pfam[[i]] <- hitx
        if (verbose == T) {
          print(i, quote = FALSE)
          print(hitx[, 1:11])
          utils::flush.console()
        }
      }
      else {
        hitx <- data.frame(id = seqi_id,
                           name = NA,
                           acc = NA,
                           desc = NA, 
                           clan = NA,
                           align_start = NA,
                           align_end = NA, 
                           model_start = NA,
                           model_end = NA,
                           ievalue = NA, 
                           cevalue = NA,
                           bitscore = NA,
                           reported = FALSE)
        pfam[[i]] <- hitx
        if (verbose == T) {
          print(hitx[, 1:11])
          utils::flush.console()
        }
      }
    }
    if(progress){
      utils::setTxtProgressBar(pb, i)
    }
    Sys.sleep(sleep)
  }
  if(progress){
    close(pb)
  }
  
  pfam <- suppressWarnings(do.call(rbind, pfam))
  pfam <- as.data.frame(pfam)
  rownames(pfam) <- 1:nrow(pfam)
  pfam$id <- as.character(pfam$id)
  pfam$name <- as.character(pfam$name)
  pfam$acc <- as.character(pfam$acc)
  pfam$desc <- as.character(pfam$desc)
  pfam$clan <- as.character(pfam$clan)
  pfam$align_start <- as.numeric(as.character(pfam$align_start))
  pfam$align_end <- as.numeric(as.character(pfam$align_end))
  pfam$model_start <- as.numeric(as.character(pfam$model_start))
  pfam$model_end <- as.numeric(as.character(pfam$model_end))
  pfam$ievalue <- as.numeric(as.character(pfam$ievalue))
  pfam$cevalue <- as.numeric(as.character(pfam$cevalue))
  pfam$bitscore <- as.numeric(as.character(pfam$bitscore))
  pfam$reported <- as.logical(as.integer(pfam$reported))
  return(pfam)
}

