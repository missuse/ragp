#' Split a fasta formatted file.
#'
#' The function splits a fasta formatted file to a defined number of smaller .fasta files for further processing.
#'
#' @param path_in A path to the .FASTA formatted file that is to be processed
#' @param path_out A path where the resulting .FASTA formatted files should be stored. The path should also contain the prefix name of the fasta files on which _n (integer from 1 to number of fasta files generated) will be appended along with the extension ".fa"
#' @param num_seq Integer defining the number of sequences to be in each resulting .fasta file. Defaults to 20000.
#' @param trim Logical, should the sequences be trimmed to 4000 amino acids to bypass the CBS server restrictions. Defaults to FALSE.
#'
#' @return  A Character vector of the paths to the resulting .FASTA formatted files
#'
#' @examples
#' \dontrun{
#' library(ragp)
#' #create a fasta file to be processed, not needed if the input file is already present
#' data(at_nsp)
#' library(seqinr)
#' write.fasta(sequence = strsplit(at_nsp$sequence, ""),
#'             name = at_nsp$Transcript.id,
#'             file = "at_nsp.fasta")
#'
#' #assumes input/output file are in working directory:
#' file_paths <- split_fasta(path_in = "at_nsp.fasta",
#'                           path_out = "at_nsp_split",
#'                           num_seq = 500)
#' }
#' 
#' @export

split_fasta <- function(path_in, path_out, num_seq = 20000, trim = FALSE){
  if (!requireNamespace("seqinr", quietly = TRUE)) {
    stop("seqinr needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!file.exists(path_in)){
    stop("cannot find file in the specified path_in")
  }
  if (missing(num_seq)){
    num_seq <- 20000
  }
  if (length(num_seq) > 1){
    num_seq <- 20000
    warning("num_seq should be of length 1, setting to default: num_seq = 20000")
  }
  if (!is.numeric(num_seq)){
    num_seq <- as.numeric(num_seq)
    warning("num_seq is not numeric, converting using 'as.numeric'")
  }
  if (is.na(num_seq)){
    num_seq <- 20000
    warning("num_seq was set to NA, setting to default: num_seq = 20000")
  }
  if (is.numeric(num_seq)){
    num_seq <- floor(num_seq)
  }
  if (missing(trim)){
    trim <- FALSE
  }
  if (length(trim) > 1){
    trim <- FALSE
    warning("trim should be of length 1, setting to default: trim = FALSE")
  }
  if (!is.logical(trim)){
    trim <- as.logical(trim)
    warning("trim is not logical, converting using 'as.logical'")
  }
  if (is.na(trim)){
    trim <- FALSE
    warning("trim was set to NA, setting to default: trim = FALSE")
  }
  temp_file <- seqinr::read.fasta(file = path_in, seqtype = "AA")
  if (trim == TRUE){
    temp_file = lapply(temp_file, function(x){
      len = length(x)
      if (len > 4000){
        out = x[1:4000]
      } else {
        out = x
      }
      return(out)
    }
    )
  }
  len <- length(temp_file)
  splt <- num_seq
  pam <- ((seq(len)-1) %/% splt)+1
  m_split <- split(temp_file, pam)
  file_list <- vector("character", length(m_split))
  for (i in 1 : length(m_split)){
    seqinr::write.fasta(sequences = m_split[[i]],
                        names = names(m_split[[i]]),
                        file.out = paste(path_out, i, ".fa", sep=""))
    file_list[i] <- paste(path_out, i, ".fa", sep="")
  }
  return(file_list)
}
