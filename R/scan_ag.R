#' Find AGII glycomodules in protein sequences
#'
#' AGII glycomodules are amino acid dimmers: OA, OS, OT, AO, SO and TO (and probably OG, OV, GO and VO) which are in close proximity to each other (Tan et al., 2003).
#' Where: O - hydroxyproline, A - alanine, S - serine, T - threnonine, G - glycine and V - valine. This function attempts to find the mentioned dimmers according to user specified rules. Since the positions of hydroxyprolines are usually unknown, all prolines are considered instead.
#'
#'@param sequence A vector of strings representing protein amino acid sequences.
#'@param id An optional vector of strings representing the names of the corresponding sequences.
#'@param dim An integer defining the minimum number of close dimmers to be considered, at default set to 3.
#'@param div An integer defining the maximum number of amino acids that can separate the dimmers for them to be considered, at default to 7
#'@param type One of c("conservative", "extended"), if conservative only A, S and T will be considered as possible P partners in dimmers, if extended dimmers involving P with A, S, T, G and V will be consdered. At default set to "extended".
#'@param exclude_ext Logical, should extensin (SPPP+) regions be excluded from the search. At default set to FALSE.
#'
#' @return A named list with components:
#' \describe{
#'   \item{$id}{Character vector, as supplied in the function call or NA if not supplied}
#'   \item{$sequence}{Character vector, each element corresponding to one input sequence, with matched letters (amino acids in dimmers that satisfy the user set conditions) in uppercase}
#'   \item{$AG_aa}{Numeric vector, each element corresponding to the number of matched letters (amino acids in dimmers that satisfy the user set conditions) in each input sequence}
#'   \item{$AG_locations}{List of numeric vectors, each element corresponding to the locations of found dimmers}
#'   \item{$total_length}{Numeric vector, with elements corresponding to the total length of the found stretches of dimmers (including the amino acids betwean dimers in a match) in each sequence}
#'   \item{$longest}{Numeric vector, with elements corresponding to the maximum length of the found stretches of dimmers (including the amino acids betwean dimers in a match) in each sequence}
#'   \item{$locations}{List of matrices, each element describing the start and end locations of the found stretches of dimmers (including the amino acids betwean dimers in a match)}
#'   \item{$dim}{Integer, as from input, default dim = 3}
#'   \item{$div}{Integer, as from input, default div = 7}
#'   \item{$type}{Character, as from input, one of c("conservative", "extended")}
#' }
#'
#'@references Tan L, Leykam JF, Kieliszewski MJ. (2003) Glycosylation motifs that direct arabinogalactan addition to arabinogalactan proteins. Plant Physiol 132: 1362-136
#'@seealso \code{\link[ragp]{calculate_bias}}
#'
#'@examples
#' data(at_nsp)
#'
#' # find all stretches of AP, SP, TP, PA, PS and PT dimmers where there are at least
#' # 3 dimmers separated by a maximim of 10 amino acids betwean each two dimmers
#' at_nsp_ag <- scan_ag(sequence = at_nsp$sequence[1:20],
#'                      id = at_nsp$Transcript.id[1:20],
#'                      dim = 3,
#'                      div = 10,
#'                      type = "conservative")
#'
#' # find all stretches of AP, SP, TP, GP, VP, PA, PS, PT PG, and PV dimmers where there
#' # are at least 2 dimmers separated by a maximim of 4 amino acids betwean them
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
#'                      type = "extended", exclude_ext = T)
#'
#' at_sp_ag_ext$sequence[at_sp_ag_ext$sequence! = at_sp_agt$sequence]
#'
#'
#'@export


scan_ag <- function(sequence, id = NULL, dim = NULL, div = NULL, type = c("conservative", "extended"), exclude_ext = NULL){
  if (is.null(id)) id <- rep(NA, length(sequence))
  id <- as.character(id)
  if (is.null(dim)) dim <- 3
  dim <- as.integer(dim)
  if (length(dim) != 1)
    stop ("dim should be an integer of length 1")
  if (is.null(div)) div <- 7
  div <- as.integer(div)
  if (length(div) != 1)
    stop("div should be an integer of length 1")
  if(missing(type)) type <- "extended"
  if (!type %in% c("conservative", "extended"))
    stop ("type should be one of: conservative or extended")
  if (type == "extended") aa <- "[ASTGV]" else aa <- "[AST]"
  sequence <- toupper(sequence)
  if (missing(exclude_ext)) exclude_ext <- F
  if (exclude_ext == T){
    sequence <- stringr::str_replace_all(sequence, "SP{3,}", tolower)
  }

  regex <- paste("(P",  aa, "P(?!" , aa, ")|", aa, "P",
                 aa, "(?!P)|", aa, "P|P", aa, ")(.{0,", div, "}(P",
                 aa, "P(?!", aa, ")|", aa, "P", aa, "(?!P)|", aa,
                 "P|P", aa, ")){", dim-1, ",}", sep="")
  if (missing(sequence)) stop ("no sequence specified")
  sequence <- as.character(sequence)
  lower <- tolower(sequence)
  locations <- stringr::str_locate_all(sequence, regex)
  n <- length(sequence)
  upper_PAST <- 1:n
  length_detected <- 1:n
  longest_detected <- 1:n
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

    } else {length_detected[i]<- 0
    longest_detected[i] <- 0}
    upper_PAST[i] <- seqi
  }
  if (type == "extended"){
    oldi <- c("CDEFHIKLMNQRWY")
    newi <- tolower(c("CDEFHIKLMNQRWY"))
    upper_PAST <- chartr(old = oldi, new = newi, upper_PAST)
    oldi_p <- paste(oldi, "P", sep="")
    newi_p <- paste(newi, "p", sep="")
    upper_PAST <- gsub(paste("(?<=[", newi_p, oldi_p,
                             "])P{1}(?=[", newi_p, oldi_p, "])",
                             sep=""), "p", upper_PAST, perl = TRUE)
    oldi_ASTGV <- paste(oldi, "ASTGV", sep="")
    newi_ASTGV <- paste(newi, "astgv", sep="")

    chary <- c("A", "S", "T", "G", "V")
    regy <- lapply(chary, function(x){
      paste("(?<=[", newi_ASTGV, oldi_ASTGV,
            "])", x, "{1}(?=[", newi_ASTGV, oldi_ASTGV,
            "])", sep="")
    })
    for (i in 1:length(regy)){
      upper_PAST <- gsub(regy[[i]], tolower(chary[i]),
                         upper_PAST, perl = TRUE)
    }

    char_any <- c("P", "A", "S", "T", "G", "V")
    regy_any <- lapply(char_any, function(x){
      paste("(?<=[a-z])[", x, "]{1}(?=[a-z])", sep = "")
    }
    )
    for (i in 1:length(regy_any)){
      upper_PAST <- gsub(regy_any[[i]], tolower(char_any[i]),
                         upper_PAST, perl = TRUE)
    }

    upper_PAST <- gsub(paste("^P(?=[", oldi_p, newi_p, "])", sep = ""),
                       "p", upper_PAST, perl = TRUE)

    regy_start <- lapply(chary, function(x){
      paste("^", x, "(?=[", newi_ASTGV, oldi_ASTGV, "])", sep="")
    })
    for (i in 1:length(regy_start)){
      upper_PAST <- gsub(regy_start[[i]], tolower(chary[i]),
                         upper_PAST, perl = TRUE)
    }

    upper_PAST <- gsub(paste("(?<=[", newi_p, oldi_p, "])P$", sep=""),
                       "p", upper_PAST, perl = TRUE)

    regy_end <- lapply(chary, function(x){
      paste("(?<=[", newi_ASTGV, oldi_ASTGV, "])", x, "$", sep = "")
    })

    for (i in 1:length(regy_end)){
      upper_PAST <- gsub(regy_end[[i]], tolower(chary[i]),
                         upper_PAST, perl = TRUE)
    }
    AG_sum <- stringr::str_count(upper_PAST, "P|A|S|T|G|V")
    AG_locations <- stringr::str_locate_all(upper_PAST,"P|A|S|T|G|V")
  } else {
    oldi <- c("CDEFHIKLMNQRWYGV")
    newi <- tolower(c("CDEFHIKLMNQRWYGV"))
    upper_PAST <- chartr(old=oldi, new=newi, upper_PAST)
    oldi_p <- paste(oldi, "P", sep="")
    newi_p <- paste(newi, "p", sep="")
    upper_PAST <- gsub(paste("(?<=[", newi_p, oldi_p,
                             "])P{1}(?=[", newi_p, oldi_p, "])",
                             sep=""), "p", upper_PAST, perl = TRUE)
    oldi_AST <- paste(oldi, "AST", sep="")
    newi_AST <- paste(newi, "ast", sep="")

    chary <- c("A", "S", "T")
    regy <- lapply(chary, function(x){
      paste("(?<=[", newi_AST, oldi_AST,
            "])", x, "{1}(?=[", newi_AST, oldi_AST,
            "])", sep="")
    })
    for (i in 1:length(regy)){
      upper_PAST <- gsub(regy[[i]], tolower(chary[i]),
                         upper_PAST, perl = TRUE)
    }
    char_any <- c("P", "A", "S", "T")
    regy_any <- lapply(char_any, function(x){
      paste("(?<=[a-z])[", x, "]{1}(?=[a-z])", sep = "")
    })
    for (i in 1:length(regy_any)){
      upper_PAST <- gsub(regy_any[[i]], tolower(char_any[i]),
                         upper_PAST, perl = TRUE)
    }
    upper_PAST <-   gsub(paste("^P(?=[", oldi_p, newi_p, "])", sep = ""),
                         "p", upper_PAST, perl = TRUE)

    regy_start <- lapply(chary, function(x){
      paste("^", x, "(?=[", newi_AST, oldi_AST, "])", sep="")
    })
    for (i in 1:length(regy_start)){
      upper_PAST <- gsub(regy_start[[i]], tolower(chary[i]),
                         upper_PAST, perl = TRUE)
    }

    upper_PAST <- gsub(paste("(?<=[", newi_p, oldi_p, "])P$", sep=""),
                       "p", upper_PAST, perl = TRUE)

    regy_end <- lapply(chary, function(x){
      paste("(?<=[", newi_AST, oldi_AST, "])", x, "$", sep = "")
    })

    for (i in 1:length(regy_end)){
      upper_PAST <- gsub(regy_end[[i]], tolower(chary[i]),
                         upper_PAST, perl = TRUE)
    }

    AG_sum <- stringr::str_count(upper_PAST, "P|A|S|T")
    AG_locations <- stringr::str_locate_all(upper_PAST,"P|A|S|T")
  }
  list <- list(id = as.character(id),
               sequence = upper_PAST,
               AG_aa = as.numeric(as.character(AG_sum)),
               AG_locations = lapply(AG_locations, function (x) unique(as.vector(x))),
               total_length = as.numeric(as.character(length_detected)),
               longest = as.numeric(as.character(longest_detected)),
               locations = locations,
               dim = dim,
               div = div,
               type = type)
  class(list) <- c("AGII", "scan", "list")
  return(list)
}
