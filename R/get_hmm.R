#' Scraping hmmer web server.
#'
#' hmmer web server offers biosequence analysis using profile hidden Markov Models. This function allows searching
#' of a protein sequence vs a profile-HMM database (Pfam-A).
#'
#' @param sequence A vector of strings representing protein amino acid sequences
#' @param id A vector of strings representing the names of the corresponding sequences
#' @param verbose Bolean wheather to print out the output for each sequence, defaults to T
#' @param sleep Numeric indicating the pause in seconds betwean server calls, at default set to 1
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{id}{Character, as suplied in the function call}
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
#'   \item{reported}{Logical, is the result reported on the hmmer site. The hmmer web server outputs more hmm profile matches than it presents to the user. Results below a certain treshold are not reported (hidden) on the site.}
#' }
#'
#'@source \url{https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan}
#'
#'@seealso \code{\link[ragp]{pfam2go}}
#'
#'@examples
#'
#'pfam_pred <- get_hmm(sequence = at_nsp$sequence[1:20],
#'                   id = at_nsp$Transcript.id[1:20])
#'
#'@export


get_hmm = function (sequence, id, verbose = NULL, sleep = NULL)
{
  if (missing(verbose)) {
    verbose <- T
  }
  if (missing(sleep)) {
    sleep <- 1
  }
  sequence <- toupper(as.character(sequence))
  id <- as.character(id)
  n <- length(sequence)
  url <- "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan"
  pfam <- list()
  for (i in 1:n) {
    seqi <- sequence[i]
    res <- httr::POST(url = url, encode = "form", body = list(isajax = "1",
                                                              seq = seqi,
                                                              hmmdb = "pfam",
                                                              alt = "Paste",
                                                              threshold = "cut_ga"))
    res <- httr::content(res, as = "parsed")
    hits <- res$results$hits
    if (length(hits) != 0) {
      hitx <- lapply(hits, function(y) {
        dom <- y$domains
        doms <- lapply(dom, function(x) {
          resi <- data.frame(name = x$alihmmname, acc = x$alihmmacc,
                             desc = x$alihmmdesc, clan = ifelse(is.null(x$clan),
                                                                NA, as.character(x$clan)), align_start = x$alisqfrom,
                             align_end = x$alisqto, model_start = x$alihmmfrom,
                             model_end = x$alihmmto, ievalue = x$ievalue,
                             cevalue = x$cevalue, bitscore = x$bitscore,
                             reported = x$is_included, stringsAsFactors = F)
          return(resi)
        })
        domx <- data.frame(do.call(rbind, doms))
        return(domx)
      })
      hitx <- do.call(rbind, hitx)
      hitx <- as.data.frame(hitx, stringsAsFactors = F)
      pfam[[i]] <- hitx
      if (verbose == T) {
        print(hitx[, 1:11])
        flush.console()
      }
      Sys.sleep(sleep)
    }
    else {
      hitx <- data.frame(name = NA, acc = NA, desc = NA,
                         clan = NA, align_start = NA, align_end = NA,
                         model_start = NA, model_end = NA, ievalue = NA,
                         cevalue = NA, bitscore = NA, reported = 0)
      pfam[[i]] <- hitx
      if (verbose == T) {
        print(hitx[, 1:11])
        flush.console()
      }
      Sys.sleep(sleep)
    }
  }
  idi <- rep(id, times = unlist(lapply(pfam, function(x) nrow(x))))
  pfam <- suppressWarnings(do.call(rbind, pfam))
  pfam <- data.frame(id = idi, pfam)
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
