#' Query hmmer web server.
#'
#' hmmer web server offers biosequence analysis using profile hidden Markov Models. This function allows searching
#' of a protein sequence vs a profile-HMM database (Pfam-A).
#'
#' @aliases get_hmm get_hmm.default get_hmm.character get_hmm.data.frame get_hmm.list get_hmm.AAStringSet
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class \code{\link[seqinr]{SeqFastaAA}} resulting from \code{\link[seqinr]{read.fasta}} call. Alternatively an \code{\link[Biostrings]{AAStringSet}} object. Should be left blank if vectors are provided to sequence and id arguments.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param verbose Boolean, whether to print out the output for each sequence, defaults to FALSE.
#' @param sleep Numeric indicating the pause in seconds between server calls, at default set to 1.
#' @param attempts Integer, number of attempts if server unresponsive, at default set to 2.
#' @param timeout Numeric, time in seconds to wait for server response.
#' @param ievalue Numeric, all sequences with independent E-value lower or equal to this value will be retained in the function output. Used to filter out low similarity matches. If set some queried sequences might be discarded from the output. Suggested values: 1e-2 - 1e-5.
#' @param bitscore Numeric, all sequences with bitscore greater or equal to this value will be retained in the function output. Used to filter out low similarity. If set some queried sequences might be discarded from the output. Suggested values: 10 - 20.
#' @param progress Boolean, whether to show the progress bar, at default set to FALSE.
#' @param ... currently no additional arguments are accepted apart the ones documented bellow.
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
#' pfam_pred
#'
#' @import seqinr
#' @import httr
#' @import xml2
#' @export

get_hmm <- function (data, ...){
  if (missing(data) || is.null(data)) get_hmm.default(...)
  else UseMethod("get_hmm")
}

#' @rdname get_hmm
#' @method get_hmm character
#' @export

get_hmm.character <- function(data,
                              ...){
  if (file.exists(data)){
    dat <- seqinr::read.fasta(file = data,
                              seqtype = "AA",
                              as.string = FALSE)
    dat <- lapply(dat,
                  paste0,
                  collapse ="")
    id <- names(dat)
    sequence <- toupper(as.character(unlist(dat)))
    sequence <- sub("\\*$",
                    "",
                    sequence)
  } else {
    stop("cannot find file in the specified path",
         call. = FALSE)
  }
  res <- get_hmm.default(sequence = sequence,
                         id = id,
                         ...)
  return(res)
}

#' @rdname get_hmm
#' @method get_hmm data.frame
#' @export

get_hmm.data.frame <- function(data,
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
  res <- get_hmm.default(sequence = sequence,
                         id = id,
                         ...)

  return(res)
}

#' @rdname get_hmm
#' @method get_hmm list
#' @export

get_hmm.list <- function(data,
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
  }
  res <- get_hmm.default(sequence = sequence,
                         id = id,
                         ...)
  return(res)
}

#' @rdname get_hmm
#' @method get_hmm default
#' @export

get_hmm.default <- function(data = NULL,
                            sequence,
                            id,
                            verbose = FALSE,
                            sleep = 1,
                            attempts = 2L,
                            timeout = 10,
                            progress = FALSE,
                            ievalue = NULL,
                            bitscore = NULL,
                            ...){
  if (missing(verbose)) {
    verbose <- FALSE
  }
  if (length(verbose) > 1){
    verbose <- FALSE
    warning("verbose should be of length 1, setting to default: verbose = FALSE",
            call. = FALSE)
  }
  if (!is.logical(verbose)){
    verbose <- as.logical(verbose)
    warning("verbose is not logical, converting using 'as.logical'",
            call. = FALSE)
  }
  if (is.na(verbose)){
    verbose <- FALSE
    warning("verbose was set to NA, setting to default: verbose = FALSE",
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
  if (!missing(ievalue)) {
    if (length(ievalue) > 1){
      stop("ievalue should be of length 1",
              call. = FALSE)
    }
    if (!is.numeric(ievalue)){
      stop("ievalue should be numeric of length 1",
           call. = FALSE)
    }
    if (is.nan(ievalue)){
      stop("ievalue should be numeric of length 1",
           call. = FALSE)
    }
    if (is.na(ievalue)){
      stop("ievalue should be numeric of length 1",
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
  if (length(sleep) > 1){
    sleep <- 1
    warning("sleep should be of length 1, setting to default: sleep = 1",
            call. = FALSE)
  }
  if (!is.numeric(sleep)){
    sleep <- as.numeric(sleep)
    warning("sleep is not numeric, converting using 'as.numeric'",
            call. = FALSE)
  }
  if (is.na(sleep)){
    sleep <- 1
    warning("sleep was set to NA, setting to default: sleep = 1",
            call. = FALSE)
  }
  if (missing(sleep)) {
    sleep <- 1
  }
  if (length(sleep) > 1){
    sleep <- 1
    warning("sleep should be of length 1, setting to default: sleep = 1",
            call. = FALSE)
  }
  if (!is.numeric(sleep)){
    sleep <- as.numeric(sleep)
    warning("sleep is not numeric, converting using 'as.numeric'",
            call. = FALSE)
  }
  if (is.na(sleep)){
    sleep <- 1
    warning("sleep was set to NA, setting to default: sleep = 1",
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
  if (missing(timeout)) {
    timeout <- 10
  }
  if (length(timeout) > 1){
    timeout <- 10
    warning("timeout should be of length 1, setting to default: timeout = 10",
            call. = FALSE)
  }
  if (!is.numeric(timeout)){
    timeout <- as.numeric(timeout)
    warning("timeout is not numeric, converting using 'as.numeric'",
            call. = FALSE)
  }
  if (is.na(timeout)){
    timeout <- 10
    warning("timeout was set to NA, setting to default: timeout = 10",
            call. = FALSE)
  }
  if (timeout < 2){
    timeout <- 10
    warning("timeout was set to less then 2, setting to default: timeout = 10",
            call. = FALSE)
  }
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
      message("\n",
              "sequence ",
              seqi_id,
              " is longer than the permitted hmmscan length",
              "\n",
              "processing in chunks")
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
      if (verbose == TRUE) {
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
        if (!is.null(res)){
          if(res$status_code == 500){ #if(res$status_code != 200){
            res <- NULL
          }
        }
        if (!is.null(res)) {
          break
        }
      }
      if (i == 1) {
        if (attempt > attempts) {
          stop("server was unresponsive, consdier increasing timeout and attempts arguments",
               call. = FALSE)
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
                        id[i], "returning finished queries"),
                  call. = FALSE)
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
  if (!missing(ievalue)) {
    pfam <- pfam[!is.na(pfam$ievalue) & pfam$ievalue <= ievalue,]
  }
  if (!missing(bitscore)) {
    pfam <- pfam[!is.na(pfam$bitscore) & pfam$bitscore >= bitscore,]
  }
  rownames(pfam) <- NULL
  return(pfam)
}

#' @rdname get_hmm 
#' @method get_hmm AAStringSet
#' @export

get_hmm.AAStringSet <-  function(data,
                                 ...){
  sequence <- as.character(data)
  id <- names(sequence)
  sequence <- unname(sequence)
  sequence <- toupper(sequence)
  sequence <- sub("\\*$",
                  "",
                  sequence)
  
  res <- get_hmm.default(sequence = sequence,
                         id = id,
                         ...)
  return(res)
}
