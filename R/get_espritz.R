#' Query Espritz web server.
#'
#' Espritz web server predicts disordered regions from primary sequence. It utilizes Bi-directional Recursive Neural Networks and can process proteins on a genomic scale with little effort and state-of-the-art accuracy.
#'
#' @aliases get_espritz get_espritz.default get_espritz.character get_espritz.data.frame get_espritz.list
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class "SeqFastaAA" resulting from \code{\link[seqinr]{read.fasta}} call. Should be left blank if vectors are provided to sequence and id arguments.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param model One of c('X-Ray', 'Disprot', 'NMR'), default is 'X-Ray'. Determines the model to be used for prediction. See details.
#' @param FPR One of c('best Sw', '5"\%" FPR'). default is 'best Sw'. Determines the cutoff probability for prediction. 'best Sw' maximizes a weighted score rewarding correctly disorder prediction more than order prediction.
#' @param simplify A Boolean indicating the type of returned object, defaults to TRUE.
#' @param ... currently no additional arguments are accepted apart the ones documented bellow.
#'
#' @return If simplify == TRUE:
#' A data frame (one row per disordered region) with columns:
#' \describe{
#'   \item{start}{Integer, indicating the sequence position of disordered region start.}
#'   \item{end}{Integer, indicating the sequence position of disordered region end.}
#'   \item{id}{Character, indicating the protein identifier.}
#'   }
#'
#' If simplify == FALSE:
#' A data frame (one row per protein) with columns:
#' \describe{
#'   \item{id}{Character, indicating the protein identifier.}
#'   \item{probability}{List column of numeric vectors, vectors contain probabilities of disorder for each residue.}
#'   \item{prediction}{Character, indicating the prediction: D - disordered, O - ordered for each residue.}
#'   }
#'
#' @details Three models trained on different data sets are available and can be selected via the argument model: X-Ray - based on missing atoms from the Protein Data Bank (PDB) X-ray solved structures. If this option is chosen then the predictors with short disorder options are executed. Disprot -  contains longer disorder segments compared to x-ray. In particular, disprot a manually curetted database which is often based on functional attributes of the disordered region was used for this definition. Disorder residues are defined if the disprot curators consider the residue to be disordered at least once. All other residues are considered ordered. If this option is chosen then the predictors with long disorder options are executed. 'NMR' - based on NMR mobility. NMR flexibility is calculated using the Mobi server optimized to replicate the ordered-disordered NMR definition used in CASP8. These models provide quite different predictions. For further details visit \url{http://protein.bio.unipd.it/espritz/help_pages/help.html} and \url{http://protein.bio.unipd.it/espritz/help_pages/methods.html}.
#'
#' @note The Espritz web server has a limit on the amount of daily queries by ip. The function will inform the user when the limit has been exceeded.
#'
#' @source \url{http://protein.bio.unipd.it/espritz/}
#' @references Walsh I, Martin AJM, Di domenico T, Tosatto SCE (2012) ESpritz: accurate and fast prediction of protein disorder. Bioinformatics 28(4): 503 - 509
#'
#' @examples
#' library(ragp)
#'
#' espritz_test <- get_espritz(at_nsp[1:10,],
#'                             sequence,
#'                             Transcript.id)
#'
#' @import seqinr
#' @import httr
#' @import stringr
#' @import xml2
#' @export get_espritz

get_espritz <- function (data, ...){
  if (missing(data) || is.null(data)) get_espritz.default(...)
  else UseMethod("get_espritz")
}

#' @rdname get_espritz
#' @method get_espritz character
#' @export

get_espritz.character <- function(data,
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
  res <- get_espritz.default(sequence = sequence,
                             id = id,
                             ...)
  return(res)
}

#' @rdname get_espritz
#' @method get_espritz data.frame
#' @export

get_espritz.data.frame <- function(data,
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
  res <- get_espritz.default(sequence = sequence,
                             id = id,
                             ...)

  return(res)
}

#' @rdname get_espritz
#' @method get_espritz list
#' @export

get_espritz.list <- function(data,
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
  res <- get_espritz.default(sequence = sequence,
                             id = id,
                             ...)
  return(res)
}

#' @rdname get_espritz
#' @method get_espritz default
#' @export

get_espritz.default <- function(data = NULL,
                                sequence,
                                id,
                                model = c("X-Ray", "Disprot", "NMR"),
                                FPR = c("best Sw", "5% FPR"),
                                simplify = TRUE,
                                ...){
  if (missing(model)) {
    model <- "X-Ray"
  }
  if (!model %in% c("X-Ray", "Disprot", "NMR")) {
    stop("model should be one of: 'X-Ray', 'Disprot', 'NMR'",
         call. = FALSE)
  }
  if (length(model) > 1){
    stop("model should be one of: 'X-Ray', 'Disprot', 'NMR'",
         call. = FALSE)
  }
  if (missing(FPR)) {
    FPR <- "best Sw"
  }
  if (!FPR %in% c("best Sw", "5% FPR")) {
    stop("FPR should be one of: 'best Sw', '5% FPR'",
         call. = FALSE)
  }
  if (length(FPR) > 1){
    stop("FPR should be one of: 'best Sw', '5% FPR'",
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
  aa_regex <- "[^ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv]"
  if (any(grepl(aa_regex, sequence))){
    warning(paste("sequences: ",
                  paste(id[grepl(aa_regex,
                                 sequence)],
                        collapse = ", "),
                  " contain symbols not corresponding to amino acids",
                  sep = ""),
            call. = FALSE)
  }
  if (length(simplify) > 1){
    simplify <- TRUE
    warning("simplify should be of length 1, setting to default: simplify = TRUE",
            call. = FALSE)
  }
  if (!is.logical(simplify)){
    simplify <- as.logical(simplify)
    warning("simplify is not logical, converting using 'as.logical'",
            call. = FALSE)
  }
  if (is.na(simplify)){
    simplify <- TRUE
    warning("simplify was set to NA, setting to default: simplify = TRUE",
            call. = FALSE)
  }

  id_num <- seq_along(sequence)
  mat <- data.frame(id = id_num,
                    sequence = sequence,
                    id_correct = id)
  splt <- 0:(nrow(mat) - 1) %/% 15000
  mat_split <- split(mat, splt)
  url_query <- "http://protein.bio.unipd.it/espritz/Espritz.jsp"
  url_2 <- "http://protein.bio.unipd.it/espritz/work/pid_"
  tot <- length(mat_split)*5
  pb <- utils::txtProgressBar(min = 0,
                              max = tot,
                              style = 3)
  res <- lapply(seq_along(mat_split), function(z) {
    matz <- mat_split[[z]]
    spltz <- 0:(nrow(matz) - 1) %/% 3000
    mat_splitz <- split(matz, spltz)
    pid <- lapply(seq_along(mat_splitz), function(i) {
      mat <- mat_splitz[[i]]
      up <- paste(">",
                  mat$id,
                  "\n",
                  mat$sequence,
                  "\n",
                  collapse = "",
                  sep = "")

      res <- httr::POST(url = url_query,
                        encode = "multipart",
                        body = list(sequence = up,
                                    model = model,
                                    FPR = FPR
                        )
      )
      res <- httr::content(res,
                           as = "parsed",
                           encoding = "UTF-8")
      really_bad <- "You have exceeded your daily number of runs or number of concurrent runs."
      if(grepl(really_bad, as.character(res))){
        stop("Your job has been blocked.",
             "\nYou have exceeded your daily",
             " number of runs or number of concurrent runs.",
             "\nPlease wait that your running jobs",
             " complete or try again tomorrow.",
             call. = FALSE)
      }

      res <- xml2::xml_find_all(res,
                                ".//*[@id='userinfo']")
      res <- xml2::xml_text(res)
      res <- gsub(".*?(\\d+).*",
                  "\\1",
                  res)
      res
    }
    )
    pid <- unlist(pid)
    res <- lapply(seq_along(pid), function(i){
      pidi <- pid[i]
      repeat{
        url_res <- paste0(url_2, pidi, "/espritz.html")
        res2 <- httr::GET(url_res)
        res2 <- as.character(httr::content(res2,
                                           as = "text",
                                           encoding = "UTF-8"))

        bad <- "Please do not close this window unless bookmarked."
        if(!grepl(bad,
                  res2)){
          break
        }
        Sys.sleep(1)
      }
      mat <- mat_splitz[[i]]
      ids <- mat$id
      id_correct <- mat$id_correct
      collected_res <- vector("list", length(ids))
      for(k in seq_along(ids)){
        idk <- ids[k]
        id_correctk <- id_correct[k]
        url_prob <- paste0(url_2,
                           pidi,
                           "//batch/",
                           idk,
                           ".fasta.espritz")
        res_prob <- utils::read.table(url_prob)
        colnames(res_prob) <- c("prediction",
                                "probability")
        res_fin <- data.frame(id = id_correctk,
                              probability = I(list(res_prob$probability)),
                              prediction = paste(res_prob$prediction,
                                                 collapse = "",
                                                 sep = ""),
                              stringsAsFactors = TRUE)

        collected_res[[k]] <- res_fin
      }
      collected_res <- do.call(rbind,
                               collected_res)
      utils::setTxtProgressBar(pb, (z-1)*5+i)
      collected_res
    }
    )
    res <- do.call(rbind,
                   res)
    if(simplify){
      start <- stringr::str_locate_all(res$prediction,
                                       "(?<=O|^)D")
      start <- lapply(start, function(x){
        if(nrow(x) == 0){
          cbind(start = NA, end = NA)
        } else {
          x
        }
      }
      )
      id_reps <- rep(res$id,
                     unlist(lapply(start,
                                   nrow)))

      start <- do.call(rbind, start)

      end <- stringr::str_locate_all(res$prediction,
                                     "D(?=O|$)")

      end <- lapply(end, function(x){
        if(nrow(x) == 0){
          cbind(start = NA, end = NA)
        } else {
          x
        }
      }
      )
      end <- do.call(rbind,
                     end)

      res <- data.frame(start = start[,1],
                        end = end[,1],
                        id = id_reps,
                        stringsAsFactors = FALSE)
    }

    utils::setTxtProgressBar(pb, tot)
    return(res)
  }
  )
  close(pb)
  res <- do.call(rbind,
                 res)
  return(res)
}
