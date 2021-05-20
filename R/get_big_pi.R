#' Query big-PI Plant Predictor web server.
#'
#' big-PI Plant Predictor is a web server utilizing a scoring algorithm for prediction of GPI modification sites in plants.
#'
#' @aliases get_big_pi get_big_pi.default get_big_pi.character get_big_pi.data.frame get_big_pi.list
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class "SeqFastaAA" resulting from \code{\link[seqinr]{read.fasta}} call. Should be left blank if vectors are provided to sequence and id arguments.
#' @param sequence An appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id An appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param simplify A boolean indicating the type of returned object, defaults to TRUE.
#' @param sleep A numeric indicating the pause in seconds between server calls, at default set to 1.
#' @param progress Boolean, whether to show the progress bar, at default set to FALSE.
#' @param ... currently no additional arguments are accepted apart the ones documented bellow.
#' 
#' @return If simplify == TRUE:
#' 
#' A data frame with columns:
#' \describe{
#'   \item{id}{Character, indicating the protein identifier}
#'   \item{is.gpi}{Logical, did big-Pi predict the presence of a GPI}
#'   \item{Quality}{Character, indicating the quality of the highest scoring omega-site}
#'   \item{omega_site}{Integer, indicating the sequence position of the highest scoring omega-site}
#'   \item{PValue}{Numeric, indicating the p-value for the prediction of the highest scoring omega-site}
#'   }
#'
#' If simplify == FALSE:
#' 
#' A list of predictions, each element named according to the sequence id, containing a two element list:
#' \describe{
#'   \item{prediction}{data frame, resembling the one returned by simplify == TRUE, along with alternative site predictions (if present)}
#'   \item{calculation}{data frame, with profile dependent and profile independent scores}
#'   }
#'   
#' @note If the server is unable to make a prediction due to non-amino acid letters or length of the sequence, the returned prediction is NA (is.bigpi column).
#'
#' @source \url{http://mendel.imp.ac.at/gpi/plant_server.html}
#' @references Eisenhaber B. Wildpaner M. Schultz CJ. Borner GHH. Dupree P. Eisenhaber F. (2003) Glycosylphosphatidylinositol lipid anchoring of plant proteins. Sensitive prediction from sequence- and genome-wide studies for Arabidopsis and rice. Plant Physiology 133(4): 1691-701
#'
#' @examples
#' library(ragp)
#' data(at_nsp)
#'
#' #indexes of some sequences in at_nsp
#' ind <- c(129, 145, 147, 160, 170)
#'
#' big_pi_pred <- get_big_pi(sequence = at_nsp$sequence[ind],
#'                           id = at_nsp$Transcript.id[ind],
#'                           simplify = FALSE)
#'
#' big_pi_pred <- get_big_pi(data = at_nsp[ind,],
#'                           sequence = sequence,
#'                           id = Transcript.id,
#'                           simplify = TRUE)
#' big_pi_pred
#' 
#' @import seqinr
#' @import httr
#' @import xml2
#' @export get_big_pi

get_big_pi <- function (data, ...){
  if (missing(data) || is.null(data)) get_big_pi.default(...)
  else UseMethod("get_big_pi")
}

#' @rdname get_big_pi
#' @method get_big_pi character
#' @export

get_big_pi.character <- function(data,
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
  res <- get_big_pi.default(sequence = sequence,
                            id = id,
                            ...)
  return(res)
}

#' @rdname get_big_pi
#' @method get_big_pi data.frame
#' @export


get_big_pi.data.frame <- function(data,
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
  res <- get_big_pi.default(sequence = sequence,
                            id = id,
                            ...)

  return(res)
}

#' @rdname get_big_pi
#' @method get_big_pi list
#' @export

get_big_pi.list <- function(data,
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
  res <- get_big_pi.default(sequence = sequence,
                            id = id,
                            ...)
  return(res)
}

#' @rdname get_big_pi
#' @method get_big_pi default
#' @export

get_big_pi.default <- function(data = NULL,
                               sequence,
                               id,
                               simplify = TRUE,
                               sleep = 1,
                               progress = FALSE,
                               ...){
  if (missing(simplify)){
    simplify <- TRUE
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
  if (missing(sleep)){
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
  id_sort <- id
  Terms <- c("Total Score",
             "Profile Score",
             "Term  0",
             "Term  1",
             "Term  2",
             "Term  3",
             "Term  4",
             "Term  5",
             "Term  6",
             "Term  7",
             "Term  8",
             "Term  9",
             "Term 10",
             "Term 11",
             "Term 12",
             "Term 13",
             "Term 14",
             "Term 15",
             "Term 16",
             "Term 17",
             "Term 18",
             "Term 19",
             "Term 20",
             "Term 21",
             "Profile independent Score")
  Term_pos <- c(1, 1, 4, 4, 4, 4, 5,
                4, 4, 2, 4, 4, 4, 4,
                rep(3, 9), 4, 1)
  extract_val <- function(x, y){
    res <- regmatches(y[grep(x, y)],
                      gregexpr(num_reg,
                               y[grep(x, y)]))
    res <- unlist(res)
    res <- as.numeric(res)
    return(res)
  }
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
  crop_1 <- ".*Use of the prediction function for VIRIDIPLANTAE"
  site_1 <- "Potential GPI-modification site was found"
  site_2 <- "Potential alternative GPI-modification site was found"
  num_reg <- "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?"

  short_seq <- which(nchar(sequence) <= 55)
  bad_seq <- grep(aa_regex, sequence)

  if (length(bad_seq) != 0){
    Best <- rep(NA, length(Terms))
    Alternative <- rep(NA, length(Terms))
    total_score <- NA
    pred_bad <- data.frame(omega_site = rep("invalid character", 2),
                           Quality = rep("None", 2),
                           PValue =  rep(NA, 2),
                           stringsAsFactors = FALSE)
    calc <- data.frame(Terms,
                       Best,
                       Alternative,
                       stringsAsFactors = FALSE)
    res_bad <- list(prediction = pred_bad,
                    calculation = calc)
    res_bad <- rep(list(res_bad), length(bad_seq))
    names(res_bad) <- id[bad_seq]
  }
  if (length(short_seq) != 0){
    Best <- rep(NA, length(Terms))
    Alternative <- rep(NA, length(Terms))
    total_score <- NA
    pred_short <- data.frame(omega_site = rep("short sequence", 2),
                             Quality = rep("None", 2),
                             PValue =  rep(NA, 2),
                             stringsAsFactors = FALSE)
    calc <- data.frame(Terms,
                       Best,
                       Alternative,
                       stringsAsFactors = FALSE)
    res_short <- list(prediction = pred_short,
                      calculation = calc)
    res_short <- rep(list(res_short), length(short_seq))
    names(res_short) <- id[short_seq]
  }

  if (length(c(short_seq, bad_seq)) != 0) {
    id_else <- id[-sort(c(bad_seq, short_seq))]
    sequence_else <- sequence[-sort(c(bad_seq, short_seq))]
  } else {
    id_else <- id
    sequence_else <- sequence
  }
  mat <- data.frame(id = id_else,
                    sequence =  sequence_else)
  splt <- 0:(nrow(mat)-1) %/% 100
  mat_split <- split(mat, splt)
  sleep <- 2
  url <- "https://mendel.imp.ac.at/gpi/cgi-bin/gpi_pred_plants.cgi"
  tot <- length(mat_split)
  if(progress){
    pb <- utils::txtProgressBar(min = 0,
                                max = tot,
                                style = 3)
  }
  res_good <- lapply(seq_along(mat_split), function(i){
    mat <- mat_split[[i]]
    up <- paste(">",
                mat$id,
                "\n", mat$sequence,
                "\n", collapse = "",
                sep = "")
    res <- httr::POST(url = url,
                      encode = "form",
                      body = list(Sequence = up)
    )
    res <- httr::content(res,
                         as = "parsed",
                         encoding = "UTF-8")
    res <- xml2::xml_find_all(res,
                              xpath = ".//pre")
    resulti <- xml2::xml_text(res)
    out <- lapply(strsplit(resulti, "Query sequence")[[1]], function(x){
      unlist(strsplit(x, "\\\n"))
    })
    Sys.sleep(sleep)
    if(progress){
      utils::setTxtProgressBar(pb, i)
    }
    return(out)
  })

  res_good <- unlist(res_good, recursive = FALSE)

  res_good <- res_good[which(grepl(crop_1,
                                   res_good))]

  res_out <- lapply(res_good, function(resulti){
    impro <- unlist(strsplit(resulti, "\\\n"))
    impro <- impro[grep(crop_1, impro):length(impro)]
    if (length(grep(site_1, impro)) > 0 &
        length(grep(site_2, impro)) > 0){
      Best <- unlist(lapply(seq_along(Term_pos), function(x){
        extract_val(x = Terms[x], y = impro)[Term_pos[x]]
      }))
      Alternative <- unlist(lapply(seq_along(Term_pos), function(x){
        extract_val(x = Terms[x], y = impro)[Term_pos[x]+1]
      }))

      position_gpi <- extract_val("Sequence position of the omega-site",
                                  impro)
      total_score <- extract_val("Total Score\\.",
                                 impro)[1]
      best_p <- extract_val("PValue = ", impro)[2]
      alt_p <- extract_val("PValue = ", impro)[4]

      site_q <- unlist(strsplit(impro[grep("Quality of the site", impro)],
                                " {2,}", perl = TRUE))[c(2, 4)]

      pred <- data.frame(omega_site = position_gpi,
                         Quality = site_q ,
                         PValue =  c(best_p, alt_p),
                         stringsAsFactors = FALSE)

      calc <- data.frame(Terms,
                         Best,
                         Alternative,
                         stringsAsFactors = FALSE)
      res <- list(prediction = pred,
                  calculation = calc)
    } else {
      if (length(grep(site_1, impro)) > 0 &
          length(grep(site_2, impro)) == 0){
        Best <- unlist(lapply(seq_along(Term_pos), function(x){
          extract_val(x = Terms[x], y = impro)[Term_pos[x]]
        }))
        Alternative <- rep(NA, length(Terms))
        position_gpi <- extract_val("Sequence position of the omega-site",
                                    impro)
        total_score <- extract_val("Total Score\\.",
                                   impro)[1]
        best_p <- extract_val("PValue = ", impro)[2]
        site_q <- unlist(strsplit(impro[grep("Quality of the site", impro)],
                                  " {2,}", perl = TRUE))[2]

        pred <- data.frame(omega_site = c(position_gpi[1], NA),
                           Quality = c(site_q[1], NA) ,
                           PValue =  c(best_p, NA),
                           stringsAsFactors = FALSE)
        calc <- data.frame(Terms,
                           Best,
                           Alternative,
                           stringsAsFactors = FALSE)
        res <- list(prediction = pred,
                    calculation = calc)
      } else {
        Best <- unlist(lapply(seq_along(Term_pos), function(x){
          extract_val(x = Terms[x], y = impro)[Term_pos[x]]
        }))
        Alternative <- rep(NA, length(Terms))
        best_p <- extract_val("PValue = ", impro)[2]
        position_gpi <- extract_val("Among all positions checked",
                                    impro)
        total_score <- extract_val("Total Score\\.",
                                   impro)[1]

        pred <- data.frame(omega_site = c(position_gpi[1], NA),
                           Quality = c("None", NA) ,
                           PValue =  c(best_p, NA),
                           stringsAsFactors = FALSE)
        calc <- data.frame(Terms,
                           Best,
                           Alternative,
                           stringsAsFactors = FALSE)
        res <- list(prediction = pred,
                    calculation = calc)

      }
    }
  }
  )
  if(progress){
    close(pb)
  }
  names(res_out) <- id_else

  if (length(short_seq) != 0){
    res_out <- c(res_out,
                 res_short)
  }

  if (length(bad_seq) != 0){
    res_out <- c(res_out,
                 res_bad)
  }

  res_out <- res_out[id_sort]
  if(simplify){
    res_out <- lapply(res_out, function(x){
      return(x$pred[1,])
    })
    res_out <- do.call(rbind, res_out)
    res_out$id <- id
    res_out$is.gpi <- res_out$Quality != "None"
    res_out$is.gpi <- ifelse(is.na(res_out$PValue),
                             NA, 
                             res_out$is.gpi)
    res_out$omega_site <- as.integer(res_out$omega_site)
    rownames(res_out) <- NULL
    res_out <- res_out[,c(4,5,2,1,3)]
  }
  return(res_out)
}
