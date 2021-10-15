#' Find AG glycomodules in protein sequences
#'
#' AG glycomodules are amino acid dipeptides: OA, OS, OT, AO, SO and TO (and probably OG, OV, GO and VO) which are in close proximity to each other (Tan et al., 2003).
#' Where: O - hydroxyproline, A - alanine, S - serine, T - threonine, G - glycine and V - valine. This function attempts to find the mentioned dipeptides according to user specified rules. Since the positions of hydroxyprolines are usually unknown, all prolines are considered instead. If any sequence from the supplied contains "O" the function will consider only true AG glycomodules.
#'
#' @aliases scan_ag scan_ag.default scan_ag.character scan_ag.data.frame scan_ag.list scan_ag.AAStringSet
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class \code{\link[seqinr]{SeqFastaAA}} resulting from \code{\link[seqinr]{read.fasta}} call. Alternatively an \code{\link[Biostrings]{AAStringSet}} object. Should be left blank if vectors are provided to sequence and id arguments.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param dim An integer defining the minimum number of close dipeptides to be considered, at default set to 3.
#' @param div An integer defining the maximum number of amino acids that can separate the dipeptides for them to be considered, at default to 10
#' @param type One of c("conservative", "extended"), if conservative only A, S and T will be considered as possible P|O partners in dipeptides, if extended dipeptides involving P|O with A, S, T, G and V will be considered. At default set to "extended".
#' @param exclude_ext One of c("no", "yes", "all"), should extensin (SPPP+) regions be excluded from the search: "no" - do not exclude SPPP+; "yes" - exclude all SPPP+; "all" - exclude all PPP+
#' @param simplify Boolean, should the function return a data frame or a list additional values.
#' @param tidy Boolean, should the function return a tidy data frame instead of a list if simplify = FALSE.
#' @param ... currently no additional arguments are accepted apart the ones documented bellow.
#'
#' @return If simplify = TRUE, a data frame with one row per sequence, containing columns:
#' \describe{
#'   \item{id}{Character, as supplied in the function call.}
#'   \item{sequence}{Character, input sequence with amino acids in dipeptides (which satisfy the user set conditions) in uppercase.}
#'   \item{AG_aa}{Integer, number of matched amino acids in dipeptides.}
#'   \item{total_length}{Integer, total length of the found stretches of dipeptides including the amino acids between dipeptides in a match.}
#'   \item{longest}{Integer, maximum length of the found stretches of dipeptides including the amino acids between dipeptides in a match.}
#' }
#'
#' if simplify = FALSE and tidy = TRUE, a data frame with one row per match, with columns:
#' \describe{
#'   \item{id}{Character, as supplied in the function call.}
#'   \item{sequence}{Character, input sequence with amino acids in dipeptides that satisfy the user set conditions in uppercase}
#'   \item{location.start}{Integer, start of a match.}
#'   \item{location.end}{Integer, end of a match.}
#'   \item{P_pos}{List column, each element is an integer vector with AG-proline positions in each match.}
#'   \item{AG_aa}{Integer, number of amino acids in dipeptides in each match}
#' }
#'
#'If simplify = FALSE and tidy = FALSE, a list with elements:
#'
#' \describe{
#'   \item{id}{Character vector, as supplied in the function call.}
#'   \item{sequence}{Character vector, each element corresponding to one input sequence, with matched letters (amino acids in dipeptides that satisfy the user set conditions) in uppercase}
#'   \item{AG_aa}{Integer vector, each element corresponding to the number of matched letters (amino acids in dipeptides that satisfy the user set conditions) in each input sequence}
#'   \item{AG_locations}{Named (by id) list of Integer vectors, each element corresponding to the locations of found dipeptides}
#'   \item{total_length}{Integer vector, with elements corresponding to the total length of the found stretches of dipeptides (including the amino acids between dipeptides in a match) in each sequence}
#'   \item{longest}{Integer vector, with elements corresponding to the maximum length of the found stretches of dipeptides (including the amino acids between dipeptides in a match) in each sequence}
#'   \item{locations}{Named (by id) list of numeric matrices, each element describing the start and end locations of the found stretches of dipeptides (including the amino acids between dipeptides in a match)}
#'   \item{dim}{Integer, as from input, default dim = 3}
#'   \item{div}{Integer, as from input, default div = 10}
#'   \item{type}{Character, as from input, one of c("conservative", "extended")}
#' }
#'
#' @note The function can be supplied with the sequences resulting from predict_hyp in which case only AG glycomodules containing O instead of P will be considered.
#'
#' @references Tan L, Leykam JF, Kieliszewski MJ. (2003) Glycosylation motifs that direct arabinogalactan addition to arabinogalactan proteins. Plant Physiol 132: 1362-136
#' @seealso \code{\link[ragp]{maab}} \code{\link[ragp]{predict_hyp}}
#'
#' @examples
#' data(at_nsp)
#'
#' # find all stretches of AP, SP, TP, PA, PS and PT dipeptides where there are at least
#' # 3 dipeptides separated by a maximum of 10 amino acids between each two dipeptides
#' at_nsp_ag <- scan_ag(sequence = at_nsp$sequence[1:20],
#'                      id = at_nsp$Transcript.id[1:20],
#'                      dim = 3,
#'                      div = 10,
#'                      type = "conservative")
#'
#' # find all stretches of AP, SP, TP, GP, VP, PA, PS, PT PG, and PV dipeptides where there
#' # are at least 2 dipeptides separated by a maximum of 4 amino acids between them
#' at_nsp_ag <- scan_ag(sequence = at_nsp$sequence[1:20],
#'                      id = at_nsp$Transcript.id[1:20],
#'                      dim = 2,
#'                      div = 4,
#'                      type = "extended")
#'
#' # check how much the results differ when extensin regions are excluded
#' at_sp_ag <- scan_ag(sequence = at_nsp$sequence,
#'                      id = at_nsp$Transcript.id,
#'                      dim = 3,
#'                      div = 6,
#'                      type = "extended")
#'
#'
#' at_sp_ag_ext <- scan_ag(sequence = at_nsp$sequence,
#'                      id = at_nsp$Transcript.id,
#'                      dim = 3,
#'                      div = 6,
#'                      type = "extended", exclude_ext = "yes")
#'
#' at_sp_ag_ext$sequence[at_sp_ag_ext$sequence != at_sp_ag$sequence]
#'
#' @import seqinr
#' @import stringr
#' @export

scan_ag <- function (data, ...){
  if (missing(data) || is.null(data)) scan_ag.default(...)
  else UseMethod("scan_ag")
}

#' @rdname scan_ag
#' @method scan_ag character
#' @export

scan_ag.character <- function(data,
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
  res <- scan_ag.default(sequence = sequence,
                         id = id,
                         ...)
  return(res)
}

#' @rdname scan_ag
#' @method scan_ag data.frame
#' @export

scan_ag.data.frame <- function(data,
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
  res <- scan_ag.default(sequence = sequence,
                         id = id,
                         ...)

  return(res)
}

#' @rdname scan_ag
#' @method scan_ag list
#' @export

scan_ag.list <- function(data,
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
  res <- scan_ag.default(sequence = sequence,
                         id = id,
                         ...)
  return(res)
}

#' @rdname scan_ag
#' @method scan_ag default
#' @export

scan_ag.default <- function(data = NULL,
                            sequence,
                            id,
                            dim = 3L,
                            div = 10L,
                            type = c("conservative", "extended"),
                            exclude_ext = c("no", "yes", "all"),
                            simplify = TRUE,
                            tidy = FALSE,
                            ...){
  if (length(dim) > 1){
    dim <- 3L
    warning("dim should be of length 1, setting to default: dim = 3",
            .call = FALSE)
  }
  if (!is.numeric(dim)){
    dim <- as.numeric(dim)
    warning("dim is not numeric, converting using 'as.numeric'")
  }
  if (is.na(dim)){
    dim <- 3L
    warning("dim was set to NA, setting to default: dim = 3",
            .call = FALSE)
  }
  if (is.numeric(dim)) {
    dim <- floor(dim)
  }
  if (length(div) > 1){
    div <- 10L
    warning("div should be of length 1, setting to default: div = 10",
            .call = FALSE)
  }
  if (!is.numeric(div)){
    div <- as.numeric(div)
    warning("div is not numeric, converting using 'as.numeric'")
  }
  if (is.na(div)){
    div <- 10L
    warning("div was set to NA, setting to default: div = 10",
            .call = FALSE)
  }
  if (is.numeric(div)) {
    div <- floor(div)
  }
  if (length(simplify) > 1){
    simplify <- TRUE
    warning("simplify should be of length 1, setting to default: simplify = TRUE",
            .call = FALSE)
  }
  if (!is.logical(simplify)){
    simplify <- as.logical(simplify)
    warning("simplify is not logical, converting using 'as.logical'")
  }
  if (is.na(simplify)){
    simplify <- TRUE
    warning("simplify was set to NA, setting to default: simplify = TRUE",
            .call = FALSE)
  }
  if (length(tidy) > 1){
    tidy <- FALSE
    warning("tidy should be of length 1, setting to default: tidy = FALSE",
            .call = FALSE)
  }
  if (!is.logical(tidy)){
    tidy <- as.logical(tidy)
    warning("tidy is not logical, converting using 'as.logical'",
            .call = FALSE)
  }
  if (is.na(tidy)){
    tidy <- FALSE
    warning("tidy was set to NA, setting to default: tidy = FALSE",
            .call = FALSE)
  }
  if (missing(type)){
    type <- "extended"
  }
  if (!type %in% c("conservative", "extended")){
    stop ("type should be one of: 'conservative', 'extended'",
          .call = FALSE)
  }
  if (missing(exclude_ext)){
    exclude_ext <- "no"
  }
  if (!exclude_ext %in% c("no", "yes", "all")){
    stop ("exclude_ext should be one of: 'no', 'yes', 'all'",
          .call = FALSE)
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
  aa_regex <- "[^ARNDCQEGHILKMFPSTWYVO]"
  if (any(grepl(aa_regex, sequence))){
    warning(paste("sequences: ",
                  paste(id[grepl(aa_regex,
                                 sequence)],
                        collapse = ", "),
                  " contain symbols not corresponding to amino acids",
                  sep = ""),
            call. = FALSE)
  }
  if (type == "extended"){
    aa <- "[ASTGV]"
    aaa <- "ASTGV"
  } else {
    aa <- "[AST]"
    aaa <- "AST"
  }
  regex <- paste("(P",  aa, "P(?!" , aa, ")|", aa, "P",
                 aa, "(?!P)|", aa, "P|P", aa, ")(.{0,", div, "}(P",
                 aa, "P(?!", aa, ")|", aa, "P", aa, "(?!P)|", aa,
                 "P|P", aa, ")){", dim-1, ",}", sep = "")

  regex2 <- paste("P",  aa, "P(?!" , aa, ")|", aa, "P",
                  aa, "(?!P)|", aa, "P|P", aa, sep = "")

  if (exclude_ext == "yes"){
    sequence <- stringr::str_replace_all(sequence,
                                         "S[PO]{3,}",
                                         tolower)
  }
  if (exclude_ext == "all"){
    sequence <- stringr::str_replace_all(sequence,
                                         "[PO]{3,}",
                                         tolower)
  }
  hyp <- "P"
  hyps <- FALSE
  if (any(grepl("O", sequence))) {
    message("sequence vector contains O, O will be considered instead of P")
    regex <- gsub("P",
                  "O",
                  regex,
                  fixed = TRUE)
    regex2 <- gsub("P",
                   "O",
                   regex2,
                   fixed = TRUE)
    hyp <- "O"
    hyps <- TRUE
  }
  lower <- tolower(sequence)
  locations <- stringr::str_locate_all(sequence,
                                       regex)
  n <- length(sequence)
  upper_PAST <- vector("character", n)
  length_detected <- vector("integer", n)
  longest_detected <- vector("integer", n)
  for (i in 1:n){
    locationi <- locations[[i]]
    seqi <- lower[i]
    if (nrow(locationi) >= 1){
      for (k in 1:nrow(locationi)){
        startk <- locationi[k,1]
        stopk <- locationi[k,2]
        substr(seqi,
               start = startk,
               stop = stopk) <- toupper(substr(seqi,
                                               start = startk,
                                               stop = stopk))
      }
      length_detected[i] <- sum((locationi[,2] + 1) - locationi[,1])
      longest_detected[i] <- max((locationi[,2] + 1) - locationi[,1])
    } else {
      length_detected[i] <- 0
      longest_detected[i] <- 0
    }
    upper_PAST[i] <- seqi
  }
  if (exclude_ext == "yes"){
    upper_PAST <- stringr::str_replace_all(upper_PAST,
                                           "S[POpo]{3,}",
                                           tolower)
  }
  if (exclude_ext == "all"){
    upper_PAST <- stringr::str_replace_all(upper_PAST,
                                           "[POpo]{3,}",
                                           tolower)
  }
  locations2 <- stringr::str_locate_all(upper_PAST,
                                        regex2)
  for (i in 1:n){
    locationi2 <- locations2[[i]]
    seqi2 <- lower[i]
    if (nrow(locationi2) >= 1){
      for (k in 1:nrow(locationi2)){
        startk <- locationi2[k,1]
        stopk <- locationi2[k,2]
        substr(seqi2,
               start = startk,
               stop = stopk) <- toupper(substr(seqi2,
                                               start = startk,
                                               stop = stopk))
      }
      upper_PAST[i] <- seqi2
    }
  }
  char_any <- c(hyp, unlist(strsplit(aaa, "")))
  pastgv <- paste0(char_any,
                   sep = "",
                   collapse = "|")
  AG_sum <- stringr::str_count(upper_PAST,
                               pastgv)
  AG_locations <- stringr::str_locate_all(upper_PAST,
                                          pastgv)
  if (simplify){
    out <- data.frame(id = as.character(id),
                      sequence = upper_PAST,
                      AG_aa = as.integer(as.character(AG_sum)),
                      total_length = as.integer(as.character(length_detected)),
                      longest = as.integer(as.character(longest_detected)),
                      stringsAsFactors = FALSE)
    return(out)
  }
  if (tidy){
    loc <- unlist(lapply(locations, nrow))

    tidy_seq <- rep(upper_PAST,
                    ifelse(loc == 0,
                           1,
                           loc))
    tidy_id <- rep(as.character(id),
                   ifelse(loc == 0,
                          1,
                          loc))
    tidy_loc <- lapply(locations,
                       function(x){
                         if(nrow(x) != 0){
                           x
                         } else {
                           cbind(start = NA,
                                 end = NA)
                         }
                       })
    tidy_location <- do.call(rbind, tidy_loc)

    substri <- substring(tidy_seq,
                         tidy_location[,1],
                         tidy_location[,2])

    if(hyps) {
      P_loc <- gregexpr("O", substri)
    } else {
      P_loc <- gregexpr("P", substri)
    }

    P_loc <- lapply(seq_along(P_loc), function(x){
      as.integer(unname(as.vector(P_loc[[x]] + tidy_location[x,1] - 1)))
    }
    )
    tidy <- data.frame(sequence = as.character(tidy_seq),
                       id = as.character(tidy_id),
                       location = tidy_location,
                       P_pos = I(P_loc),
                       length = tidy_location[,2] - tidy_location[,1] + 1,
                       AG_aa = as.integer(stringr::str_count(substri, pastgv)),
                       stringsAsFactors = FALSE)
    attributes(tidy)
    return(tidy)
  } else {
    AG_locations <- lapply(AG_locations,
                           function (x) as.integer(unique(as.vector(x))))
    names(AG_locations) <- as.character(id)
    names(locations) <- as.character(id)
    lst <- list(id = as.character(id),
                sequence = upper_PAST,
                AG_aa = as.integer(as.character(AG_sum)),
                AG_locations = AG_locations,
                total_length = as.integer(as.character(length_detected)),
                longest = as.integer(as.character(longest_detected)),
                locations = locations,
                dim = as.integer(dim),
                div = as.integer(div),
                type = type)
    return(lst)
  }
}

#' @rdname scan_ag
#' @method scan_ag AAStringSet
#' @export

scan_ag.AAStringSet <-  function(data,
                                 ...){
  sequence <- as.character(data)
  id <- names(sequence)
  sequence <- unname(sequence)
  sequence <- toupper(sequence)
  sequence <- sub("\\*$",
                  "",
                  sequence)
  
  res <- scan_ag.default(sequence = sequence,
                         id = id,
                         ...)
  return(res)
}
