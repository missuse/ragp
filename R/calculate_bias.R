#' Compositional bias of a protein sequence
#'
#' Calculate the compositional bias of whole protein sequences, or slide along the sequences and calculate the bias for all subsequences of a defined length (frame). Sequence motifs can also be searched.
#'
#' @param sequence A vector of strings representing protein amino acid sequences
#' @param id An optional vector of strings representing the names of the corresponding sequences
#' @param frame An optional integer defining the frame length for sliding along the sequences
#' @param type One of c("past", "pvyk", "psky", "user"), if set to "past" P, A, S and T will be consdered (AGP bias), if type "pvyk", P, V, Y, and K (PRP bias) will be consdered, if type "psky" P, S, K and Y (EXT bias) will be consdered and if type  set to "user" a user specified input must be provided
#' @param user An optional character vector, if type set to "user", vector should contain amino acid one letter symbols eg. c("P", "V", "K", "C", "Y", "T") or strings corresponding to amino acid motifs eg c("PTYK", "PVKC").
#' @param simplify Logical, should the function return a data.frame instead of a list. If simplify set to TRUE only the maximum bias_percent and bias_sum from all frames in a sequence are kept in the output. Defaults to FALSE.
#'
#' @return if simplify == F a named list with components:
#' \describe{
#' \item{$bias_sum}{An integer vector, or list of vectors (if frame specified), each element indicating the number of specified amino acids in each input sequence}
#' \item{$seq_len}{An integer vector, each element indicating the length of each input sequence}
#' \item{$bias_percent}{An integer vector, or list of vectors (if frame specified), each element indicating the percent of specified amino acids in each input sequence or per frame if frame specified}
#' \item{$id}{Character vector, as supplied in the function call or NA if not supplied}
#' \item{$frame}{Integer, as supplied in the function call or NULL if not supplied}
#' \item{$type}{Character, as supplied in the function call, defaults to "past"}
#' \item{$user}{Character vector, as supplied in the function call or NULL if not supplied}
#'}
#'if simplify == TRUE a data frame with column names corresponding to the list components.
#'
#' @seealso \code{\link[ragp]{scan_ag}}
#'
#' @examples
#'
#' library(ragp)
#' data(at_nsp)
#'
#' test_bias <- calculate_bias(sequence = at_nsp$sequence[1:20],
#'                             id = at_nsp$Transcript.id[1:20])
#'
#' #when frame is not specified the output can be converted to a data frame for easier manipulation
#' as.data.frame(do.call(cbind, test_bias))
#'
#' test_bias <- calculate_bias(sequence = at_nsp$sequence[1:20],
#'                             id = at_nsp$Transcript.id[1:20],
#'                             frame = 20,
#'                             user = c("P", "V", "K", "C", "Y", "T"))
#'
#' #for motif search input a string instead of a character vector to user argument
#' test_bias <- calculate_bias(sequence = at_nsp$sequence,
#'                             id = at_nsp$Transcript.id,
#'                             user = "PTYK")
#'
#' #count all sequencs with a "PAST" bias over 45%
#' test_bias <- calculate_bias(sequence = at_nsp$sequence,
#'                             id = at_nsp$Transcript.id,
#'                             simplify = TRUE)
#' nrow(test_bias[test_bias$bias_percent >= 45, ])
#'
#' # count all sequencs with 2 or more SPPP motifs
#' test_bias <- calculate_bias(sequence = at_nsp$sequence,
#'                             id = at_nsp$Transcript.id,
#'                             simplify = TRUE,
#'                             user = "SPPP")
#' nrow(test_bias[test_bias$bias_sum >= 2,])
#'
#' @export

calculate_bias <- function(sequence, id = NULL, frame = NULL, type = c("past", "pvyk", "psky", "user"), user = NULL, simplify = NULL){
  if (is.null(id)) id <- rep(NA, length(sequence))
  id <- as.character(id)
  if (missing(type)) type <- "past"
  if (!missing(user)) type <- "user"
  if (!type %in% c("past", "pvyk", "psky", "user")){
    stop ("type should be one of: past, pvyk, psky or user")
    }
  if (type == "user" & missing(user)){
    stop ("please specify a vector of amino acids. Symbols to be consdered: A R N D C Q E G H I L K M F P S T W Y V")
    }
  amino_acids = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

  if (!missing(user)){
    if (!any(unlist(strsplit(user, "")) %in% amino_acids)){
    stop ("type = user, please specify a vector of amino acids. Symbols to be consdered: A R N D C Q E G H I L K M F P S T W Y V")
    }
  }
  if (is.null(simplify)) simplify <- F
  sequence <- toupper(as.character(sequence))
  if (type == "past") aa <- c("P", "A", "S", "T")
  if (type == "pvyk") aa <- c("P", "V", "Y", "K")
  if (type == "psky") aa <- c("P", "S", "K", "Y")
  if (type == "user") aa <- user
  if (is.null(frame)){
    bias_sum <- stringr::str_count(sequence, paste(aa, sep = "" , collapse = "|"))
    seq_len <- nchar(sequence)
    bias_percent <- bias_sum / seq_len * 100
    id <- as.character(id)
    frame <- NULL
    } else {
      bias_sum <- lapply(sequence, function(x){
        if (frame > nchar(x)) frame <- nchar(x)
        past_count <- vector("numeric", length = (nchar(x) - frame + 1))
        for (i in 1 : nchar(x) - frame + 1){
          past_count[i] <- stringr::str_count(substr(x,
                                                     start = i - 1,
                                                     stop = i - 1 + frame),
                                              paste(aa, sep = "" , collapse = "|"))
          }
        return(past_count)
        }
        )
      seq_len <- nchar(sequence)
      bias_percent <- lapply(bias_sum, function(x) x / frame * 100)
      id <- as.character(id)
      frame <- frame
      if (simplify == TRUE) {
        bias_sum <- unlist(lapply(bias_sum, max))
        bias_percent <- unlist(lapply(bias_percent, max))
      }
    }

  out = list(bias_sum = bias_sum,
             seq_len = seq_len,
             bias_percent = bias_percent,
             id = id,
             frame = frame,
             type = type,
             user = user)
  if (simplify == TRUE) {
    out <- data.frame(bias_sum = out$bias_sum,
                      seq_len = out$seq_len,
                      bias_percent = out$bias_percent,
                      id = out$id)
    out$bias_sum <- as.numeric(as.character(out$bias_sum))
    out$bias_percent <- as.numeric(as.character(out$bias_percent))
    out$seq_len <- as.numeric(as.character(out$seq_len))
  }
  return (out)
}
