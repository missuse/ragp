#' Scraping Phobius web server.
#'
#' Phobius web server is a combined transmembrane topology and signal peptide (N-sp) predictor. Currently only "normal prediction" of signal peptides is supported by the function.
#'
#' @param file path to the fasta formated file to be analyzed
#'
#' @return  A data frame with columns:
#' \describe{
#' \item{Name}{Character, name of the submitted sequence.}
#' \item{tm}{Integer, the number of predicted transmembrane segments.}
#' \item{sp}{Character, Y/0 indicator if a signal peptide was predicted or not.}
#' \item{prediction}{Character string, predicted topology of the protein.}
#' \item{cut_site}{Integer, first amino acid after removal of the signal peptide}
#'}
#'
#' @details
#' The topology (prediction column of the result) is given as the position of the transmembrane helices separated by 'i' if the loop is on the cytoplasmic or 'o' if it is on the non cytoplasmic side. A signal peptide is given by the position of its h-region separated by a n and a c, and and the position of the last amino acid in the signal peptide and the first of the mature protein separated by a /.
#'
#' @source \url{http://phobius.sbc.su.se/}
#'
#' @references Kall O. Krogh A. Sonnhammer E. L. L. (2004) A Combined Transmembrane Topology and Signal Peptide Prediction Method. Journal of Molecular Biology 338(5): 1027-1036
#' @seealso \code{\link[ragp]{get_signalp_file}}
#'
#' @examples
#'
#' data(at_nsp)
#'
#' #create a fasta file:
#' library(seqinr)
#' write.fasta(sequence = strsplit(at_nsp$sequence, ""), name = at_nsp$Transcript.id, file = "at_nsp.fasta")
#'
#' #get predictions
#' phobius_pred <- get_phobius_file("at_nsp.fasta")
#' @export


get_phobius_file = function(file){
  file_list = ragp::split_fasta(path_in = file, path_out = "temp_phob_", num_seq = 500)
  len = length(file_list)
  url <- "http://phobius.binf.ku.dk/cgi-bin/predict.pl"
  collected_res = vector("list", len)
  for (i in 1 : len){
    file_up <- httr::upload_file(file_list[i])
    res <- httr::POST(url = url,
                      encode = "multipart",
                      body = list(`protfile` = file_up ,
                                  `format` = "short"))
    res <- httr::content(res, as = "parsed")
    res <- xml2::xml_text(res, "//pre")
    res <- unlist(strsplit(as.character(res), "\n"))
    res <- res[(grep("SEQENCE", res)+1) : (grep("_uacct", res)-1)]
    res <- strsplit(res, " +")
    res <- do.call(rbind, res)
    res <- as.data.frame(res, stringsAsFactors = FALSE)
    colnames(res) <- c("Name", "tm", "sp", "prediction")
    collected_res[[i]] <- res
    unlink(file_list[i])
  }
  collected_res <- do.call(rbind, collected_res)
  collected_res$cut_site <- stringr::str_extract(collected_res$prediction,
                                                 "(?<=/)\\d+")
  return(collected_res)
}

