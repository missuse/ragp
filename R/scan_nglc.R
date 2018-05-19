#' Detection of motifs for N-glycosylation on asparagine residues.
#'
#' Detection is based on PROSITE pattern PS00001. Mean local hydrophilicity (Hopp and Woods, 1981) is used to assess if the asparagines are buried.
#'
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class "SeqFastaAA" resulting from seqinr::read.fasta call.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param span An integer specifying how many amino acids around the target asparagine residues is used to calculate hydrophilicity. At default set to 5: asparagine position - 5 to asparagine position +5 residues. Range to consider: 3 - 10. Acceptable values are 0 - 20.
#' @param cutoff An integer specifying the cutoff value for hydrophilicity. Range to consider: -1 - 1. Values lower then -3.4 exclude hydrophilicity as a parameter, while values higher than 3 result in no motifs being found.
#' @param nsp An integer, either of length 1 or length equal the number of sequences, specifying the number of N-terminal amino acids to exclude. 
#' 
#' @return A data frame with columns:
#' \describe{
#'   \item{id}{Character, as supplied in the function call.}
#'   \item{align_start}{Integer, start of motif match.}
#'   \item{motif}{Character, the motif matched.}
#'   \item{hydrophilicity}{the average hydrophilicity.}
#'   \item{is.nglc}{Boolean, is the N-glycosylation likely based on hydrophilicity.}
#'   \item{nsp}{Optional integer column provided when nsp argument is of equal length to the number of input sequences}
#' }
#'
#' @note For N-glycosylation to happen the protein must enter the endoplasmic reticulum. Please check if the proteins are likely to contain a N-terminal signal peptide. The motif Asp-Xaa-Ser/Thr (where Xaa is not Pro) on which N-glycosylation occurs is relatively common, however for N-glycosylation to occur the motif needs to be located on the protein surface. Mean local hydrophilicity (Hopp and Woods, 1981) is used here to evaluate if the asparagines are located in a hydrophilic surrounding which is more likely on the protein surface.
#'
#' @author original R code by Thomas Shafee, modified by Milan Dragićević
#'
#' @references Hopp TP. Woods KR. (1981) Prediction of protein antigenic determinants from amino acid sequences. Proceedings of the National Academy of Sciences of the United States of America, 78(6): 3824-8
#'
#' @source \url{https://prosite.expasy.org/PDOC00001}
#'
#' @examples
#' library(ragp)
#' data(at_nsp)
#'
#' nglc_pred <- scan_nglc(data = at_nsp,
#'                        sequence = sequence,
#'                        id = Transcript.id)
#'
#' @export

scan_nglc <- function(data = NULL, sequence, id, span = 5L, cutoff = 0, nsp = 15L){
  if (length(span) > 1){
    span <- 5
    warning("span should be of length 1, setting to default: span = 5")
  }
  if (!is.numeric(span)){
    span <- as.numeric(span)
    warning("span is not numeric, converting using 'as.numeric'")
  }
  if (is.na(span)){
    span <- 5
    warning("span was set to NA, setting to default: span = 5")
  }
  if (is.numeric(span)) {
    span <- floor(span)
  }
  if (!(span %in% 0:20)) {
    span <- 5
    warning(paste("Illegal span input, span will be set to 5"))
  }
  if (length(cutoff) > 1){
    cutoff <- 0
    warning("cutoff should be of length 1, setting to default: cutoff = 0")
  }
  if (!is.numeric(cutoff)){
    cutoff <- as.numeric(cutoff)
    warning("cutoff is not numeric, converting using 'as.numeric'")
  }
  if (is.na(cutoff)){
    cutoff <- 0
    warning("cutoff was set to NA, setting to default: cutoff = 0")
  }
  if (is.numeric(cutoff)) {
    cutoff <- floor(cutoff)
  }
  if (cutoff < -3.4) {
    warning(paste("cutoff < 3.4 , hydrophilicity will not be considered when determining motifs"))
  }
  if (cutoff > 3) {
    warning(paste("cutoff > 3, hydrophilicity set to high, no motifs can be found"))
  }
  if (length(nsp) > 1){
    if (length(nsp) != length(sequence)){
      nsp <- 15
    warning("nsp should be of length 1, or equal to the number of sequences, setting to default: nsp = 15")
    }
  }
  if (!is.numeric(nsp)){
    stop("nsp is not numeric")
  }
  if (is.numeric(span)) {
    span <- floor(span)
  }
  if (missing(data)) {
    if (missing(sequence)) {
      stop("protein sequence must be provided to obtain predictions")
    }
    if (missing(id)) {
      stop("protein id must be provided to obtain predictions")
    }
    id <- as.character(id)
    sequence <- toupper(as.character(sequence))
    if (length(sequence) != length(id)) {
      stop("id and sequence vectors are not of same length")
    }
  }
  if (class(data[[1]]) == "SeqFastaAA") {
    dat <- lapply(data, paste0, collapse = "")
    id <- names(dat)
    sequence <- toupper(as.character(unlist(dat)))
  }
  if (class(data) == "data.frame") {
    if (missing(sequence)) {
      stop("the column name with the sequences must be specified")
    }
    if (missing(id)) {
      stop("the column name with the sequence id's must be specified")
    }
    id <- as.character(substitute(id))
    sequence <- as.character(substitute(sequence))
    if (length(id) != 1L) {
      stop("only one column name for 'id' must be specifed")
    }
    if (length(sequence) != 1L) {
      stop("only one column name for 'sequence' must be specifed")
    }
    id <- if (id %in% colnames(data)) {
      data[[id]]
    }
    else {
      stop("specified 'id' not found in data")
    }
    id <- as.character(id)
    sequence <- if (sequence %in% colnames(data)) {
      data[[sequence]]
    }
    else {
      stop("specified 'sequence' not found in data")
    }
    sequence <- toupper(as.character(sequence))
  }
  if (class(data) == "character") {
    if (file.exists(data)) {
      dat <- seqinr::read.fasta(file = data, seqtype = "AA", 
                                as.string = FALSE)
      dat <- lapply(dat, paste0, collapse = "")
      id <- names(dat)
      sequence <- toupper(as.character(unlist(dat)))
    }
    else {
      stop("cannot find file in the specified path")
    }
  }
  sequence <- sub("\\*$", "", sequence)
  hoop_woods <- c("A" = -0.5,
                  "R" = 3,
                  "N" = 0.2,
                  "D" = 3,
                  "C" = -1,
                  "Q" = 0.2,
                  "E" = 3,
                  "G" = 0,
                  "H" = -0.5,
                  "I" = -1.8,
                  "L" = -1.8,
                  "K" = 3,
                  "M" = -1.3,
                  "F" = -2.5,
                  "P" = 0,
                  "S" = 0.3,
                  "T" = -0.4,
                  "W" = -3.4,
                  "Y" = -2.3,
                  "V" = -1.5)
  aa_regex <- "[^ARNDCQEGHILKMFPSTWYV]"
  bad_seq <- grepl(aa_regex, sequence)
  if (any(bad_seq)) {
    warning("sequences: ",
            paste(id[bad_seq], collapse = ", "),
            " contain symbols not corresponding to amino acids")
  }

  nglc <- gregexpr('N[^P][ST]',
                   sequence,
                   ignore.case = FALSE)
  countnglc <- unlist(lapply(nglc,length))
  
  motif <- regmatches(sequence,
                      nglc)
  
  motif[sapply(motif,
               function(x) length(x)==0L)] <- NA
  
  nglc2 <- data.frame(id = rep(id, countnglc),
                      sequence = rep(sequence, countnglc),
                      align_start = as.integer(unlist(nglc)),
                      motif = as.character(unlist(motif)))
  

  nglc2$align_start[nglc2$align_start == -1] = NA_integer_

  seq_sub <- substring(nglc2$sequence,
                       first = nglc2$align_start - span,
                       last = nglc2$align_start + span)

  hoop <- sapply(strsplit(seq_sub, ""), function(x){
    mean(hoop_woods[match(x, names(hoop_woods))])
    
  })
  nglc2$hydrophilicity <- hoop
  nglc2 <- nglc2[,-2]
  nglc2$substr <- seq_sub
  nglc2$is.nglc <- nglc2$hydrophilicity >= cutoff
  nglc2$is.nglc[is.na(nglc2$is.nglc)] <- FALSE
  if(length(nsp) > 1) {
    nsp <- rep(nsp, countnglc)
    nsp[is.na(nsp)] <- 0
    nglc2$nsp <- as.integer(nsp)
  }
  nglc2$is.nglc[nglc2$align_start < nsp] <- FALSE
  return(nglc2)
}