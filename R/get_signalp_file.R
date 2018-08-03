#' Scraping SignalP web server.
#'
#' SignalP 4.1 server predicts the presence and location of signal peptide cleavage sites in amino acid sequences from different organisms: Gram-positive prokaryotes, Gram-negative prokaryotes, and eukaryotes. The method incorporates a prediction of cleavage sites and a signal peptide/non-signal peptide prediction based on a combination of several artificial neural networks.
#'
#' @param file A path to the FASTA formatted file that is to be processed
#' @param org_type One of c("euk", "gram-", "gram+"), defaults to "euk". Which model should be used for prediction.
#' @param Dcut_type One of c("default", "sensitive", "user"), defaults to "default". The default cutoff values for SignalP 4 are chosen to optimize the performance measured as Matthews Correlation Coefficient (MCC). This results in a lower sensitivity (true positive rate) than SignalP 3.0 had. Setting this argument to "sensitive" will yield the same sensitivity as SignalP 3.0. This will make the false positive rate slightly higher, but still better than that of SignalP 3.0.
#' @param Dcut_noTM A numeric value, with range 0 - 1, defaults to 0.45. For experimenting with cutoff values.
#' @param Dcut_TM A numeric value, with range 0 - 1, defaults to 0.5. For experimenting with cutoff values.
#' @param method One of c("best", "notm"), defaults to "best". Signalp 4.1 contains two types of neural networks. SignalP-TM has been trained with sequences containing transmembrane segments in the data set, while SignalP-noTM has been trained without those sequences. Per default, SignalP 4.1 uses SignalP-TM as a preprocessor to determine whether to use SignalP-TM or SignalP-noTM in the final prediction (if 4 or more positions are predicted to be in a transmembrane state, SignalP-TM is used, otherwise SignalP-noTM). An exception is Gram-positive bacteria, where SignalP-TM is used always. If you are confident that there are no transmembrane segments in your data, you can get a slightly better performance by choosing "Input sequences do not include TM regions", which will tell SignalP 4.1 to use SignalP-noTM always.
#' @param minlen An integer value corresponding to the minimal predicted signal peptide length, at default set to 10. SignalP 4.0 could, in rare cases, erroneously predict signal peptides shorter than 10 residues. These errors have in SignalP 4.1 been eliminated by imposing a lower limit on the cleavage site position (signal peptide length). The minimum length is by default 10, but you can adjust it. Signal peptides shorter than 15 residues are very rare. If you want to disable this length restriction completely, enter 0 (zero).
#' @param trunc An integer value corresponding to the N-terminal truncation of input sequence, at default set to 70. By default, the predictor truncates each sequence to max. 70 residues before submitting it to the neural networks. If you want to predict extremely long signal peptides, you can try a higher value, or disable truncation completely by entering 0 (zero).
#' @param splitter An integer indicating the number of sequences to be in each .fasta file that is to be sent to the server. Defaults to 500. Change only in case of a server side error. Accepted values are in range of 1 to 2000.
#' @param sleep A numeric indicating the pause in seconds between POST and GET server calls, at default set to 1s. Decreasing is not recommended.
#' @return  A data frame with columns:
#' \describe{
#'   \item{Cmax}{Numeric, C-score (raw cleavage site score). The output from the CS networks, which are trained to distinguish signal peptide cleavage sites from everything else. Note the position numbering of the cleavage site: the C-score is trained to be high at the position immediately after the cleavage site (the first residue in the mature protein).}
#'   \item{Cmax.pos}{Integer, position of Cmax. position immediately after the cleavage site (the first residue in the mature protein).}
#'   \item{Smax}{Numeric, S-score (signal peptide score). The output from the SP networks, which are trained to distinguish positions within signal peptides from positions in the mature part of the proteins and from proteins without signal peptides.}
#'   \item{Smax.pos}{Integer, position of Smax}
#'   \item{Ymax}{Numeric, Y-score (combined cleavage site score), A combination (geometric average) of the C-score and the slope of the S-score, resulting in a better cleavage site prediction than the raw C-score alone. This is due to the fact that multiple high-peaking C-scores can be found in one sequence, where only one is the true cleavage site. The Y-score distinguishes between C-score peaks by choosing the one where the slope of the S-score is steep.}
#'   \item{Ymax.pos}{Integer, position of Ymax}
#'   \item{Smean}{Numeric, The average S-score of the possible signal peptide (from position 1 to the position immediately before the maximal Y-score)}
#'   \item{Dmean}{Numeric, D-score (discrimination score). A weighted average of the mean S and the max. Y scores. This is the score that is used to discriminate signal peptides from non-signal peptides.}
#'   \item{is.sp}{Character, does the sequence contain a N-sp}
#'   \item{Dmaxcut}{Numeric, as from input, Dcut_noTM if SignalP-noTM network used and Dcut_TM if SignalP-TM network used}
#'   \item{Networks.used}{Character, which network was used for the prediction: SignalP-noTM or SignalP-TM}
#'   \item{id}{Character, as from input}
#'   \item{is.signalp}{Logical, did SignalP predict the presence of a signal peptide}
#'   }
#'
#' @source \url{http://www.cbs.dtu.dk/services/SignalP-4.1/}
#' @references Petersen TN. Brunak S. Heijne G. Nielsen H. (2011) SignalP 4.0: discriminating signal peptides from transmembrane regions. Nature Methods 8: 785-786
#'
#' @seealso \code{\link[ragp]{get_targetp_file}}
#'
#' @examples
#' \dontrun{
#' library(seqinr)
#' write.fasta(sequence = strsplit(at_nsp$sequence, ""),
#'             name = at_nsp$Transcript.id, file = "at_nsp.fasta")
#' signalp_pred <- get_signalp_file("at_nsp.fasta")
#' #should provide predictions for all 2700 sequences in under 1 min
#' }
#' @export


get_signalp_file = function (file, org_type = c("euk", "gram-", "gram+"), Dcut_type = c("default",
                                                                                        "sensitive", "user"), Dcut_noTM = 0.45, Dcut_TM = 0.5, method = c("best",
                                                                                                                                                           "notm"), minlen = NULL, trunc = NULL, splitter = 500, sleep = 3)
{
  message("This function has been deprecated, please use get_signalp")
  if (!requireNamespace("seqinr", quietly = TRUE)) {
    stop("seqinr needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (missing(splitter)) {
    splitter <- 500
  }
  if (length(splitter) > 1){
    splitter <- 500
    warning("splitter should be of length 1, setting to default: splitter = 500")
  }
  if (!is.numeric(splitter)){
    splitter <- as.numeric(splitter)
    warning("splitter is not numeric, converting using 'as.numeric'")
  }
  if (is.na(splitter)){
    splitter <- 500
    warning("splitter was set to NA, setting to default: splitter = 500")
  }
  if (is.numeric(splitter)) {
    splitter <- floor(splitter)
  }
  if (!(splitter %in% 1:2000)) {
    splitter <- 500
    warning(paste("Illegal splitter input, splitter will be set to 500"))
  }
  if (missing(sleep)) {
    sleep <- 3
  }
  if (length(sleep) > 1){
    sleep <- 3
    warning("sleep should be of length 1, setting to default: sleep = 3")
  }
  if (!is.numeric(sleep)){
    sleep <- as.numeric(sleep)
    warning("sleep is not numeric, converting using 'as.numeric'")
  }
  if (is.na(sleep)){
    sleep <- 3
    warning("sleep was set to NA, setting to default: sleep = 3")
  }
  if (sleep < 2){
    warning("setting sleep to less than 2s can cause problems when fetching results from the server")
  }
  if (missing(org_type)) {
    org_type <- "euk"
  }
  if (!org_type %in% c("euk", "gram-", "gram+")) {
    stop("org_type should be one of: 'euk', 'gram-', 'gram+'")
  }
  if (length(org_type) > 1){
    stop("org_type should be one of: 'euk', 'gram-', 'gram+'")
  }
  if (missing(Dcut_type)) {
    Dcut_type <- "default"
  }
  if (!Dcut_type %in% c("default", "sensitive", "user")) {
    stop("Dcut_type should be one of: 'default', 'sensitive', 'user'")
  }
  if (length(Dcut_type) > 1){
    stop("Dcut_type should be one of: 'default', 'sensitive', 'user'")
  }
  if (missing(Dcut_noTM)) {
    Dcut_noTM <- "0.45"
  }  else {
    Dcut_noTM <- as.character(Dcut_noTM)[1]
  }
  if (!is.numeric(as.numeric(Dcut_noTM))){
    Dcut_noTM <- "0.45"
    warning("Dcut_noTM could not be converted to numeric, setting to default: Dcut_noTM = '0.45'")
  }
  if (is.na(Dcut_noTM)) {
    Dcut_noTM <- "0.45"
    warning("Dcut_noTM was set to NA, setting to default: Dcut_noTM = '0.45'")
  }
  if (as.numeric(Dcut_noTM[1]) > 1) {
    Dcut_noTM <- "0.45"
    warning("Dcut_noTM must take values in the range 0 - 1,
            it was set to the default: Dcut_noTM = '0.45'")
  }
  if (as.numeric(Dcut_noTM[1]) < 0) {
    Dcut_noTM <- "0.45"
    warning("Dcut_noTM must take values in the range 0 - 1,
            it was set to the default: Dcut_noTM = '0.45'")
  }    
  if (missing(Dcut_TM)) {
    Dcut_TM <- "0.5"
  } else {
    Dcut_TM <- as.character(Dcut_TM)[1]
  }
  if (!is.numeric(as.numeric(Dcut_TM))){
    Dcut_TM <- "0.5"
    warning("Dcut_TM could not be converted to numeric, setting to default: Dcut_TM = '0.5'")
  }
  if (is.na(Dcut_TM)) {
    Dcut_TM <- "0.5"
    warning("Dcut_noTM was set to NA, setting to default: Dcut_TM = '0.5'")
  }
  if (as.numeric(Dcut_TM[1]) > 1) {
    Dcut_TM <- "0.5"
    warning("Dcut_TM must take values in the range 0 - 1,
            it was set to the default: Dcut_TM = '0.5'")
  }
  if (as.numeric(Dcut_TM[1]) < 0) {
    Dcut_TM <- "0.5"
    warning("Dcut_TM must take values in the range 0 - 1,
            it was set to the default: Dcut_TM = '0.5'")
  } 
  if (missing(method)) {
    method <- "best"
  }
  if (!method %in% c("best", "notm")){
    stop("method should be one of: 'best', 'notm'")
  }
  if (length(method) > 1){
    stop("method should be one of: 'best', 'notm'")
  }
  if (missing(minlen)) {
    minlen <- ""
  }  else {
    minlen <- as.character(minlen)[1]
  }
  if (missing(trunc)) {
    trunc <- ""
  }  else {
    trunc <- as.character(trunc)[1]
  }
  if (!file.exists(file)){
    stop("cannot find file in the specified adress")
  }
  url <- "http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi"
  file_list <- ragp::split_fasta(path_in = file,
                                 path_out = "temp_signalp_",
                                 num_seq = splitter)
  jobid <- vector("character", length(file_list))
  for (i in 1:length(file_list)) {
    file_up <- httr::upload_file(file_list[i])
    res <- httr::POST(url = url, encode = "multipart", body = list(configfile = "/usr/opt/www/pub/CBS/services/SignalP-4.1/SignalP.cf",
                                                                   SEQSUB = file_up, orgtype = org_type, `Dcut-type` = Dcut_type,
                                                                   `Dcut-noTM` = Dcut_noTM, `Dcut-TM` = Dcut_TM, graphmode = NULL,
                                                                   format = "short", minlen = minlen, method = method,
                                                                   trunc = trunc))
    res <- httr::content(res, as = "parsed")
    res <- xml2::xml_find_all(res,
                               ".//input[@name='jobid']")
    jobid[i] <- xml2::xml_attr(res,
                               "value")
    Sys.sleep(3)
  }
  for (i in 1:length(file_list)) {
    unlink(file_list[i])
  }
  collected_res = vector("list", length(jobid))
  for (i in 1:length(jobid)) {
    repeat {
      res2 <- httr::GET(url = url,
                        query = list(jobid = jobid[i],
                                     wait = "20"))
      bad = xml2::xml_text(xml2::xml_find_all(httr::content(res2,
                                                            as = "parsed"), "//head"))
      if (grepl("Illegal", bad)) {
        prt = xml2::xml_text(xml2::xml_find_all(httr::content(res2,
                                                              as = "parsed"), "//li"))
        stop(paste0(prt, ". Problem in file: ", "temp_",
                    i, ".fa"))
      }
      res2 <- as.character(
        xml2::xml_find_all(
          httr::content(res2,
                        as = "parsed"),
          ".//pre")
      )
      res2_split <- unlist(strsplit(res2, "\n"))
      if (any(grepl("Cmax", res2_split))) {
        break
      }
    }
    res2_split = res2_split[(which(grepl("name", res2_split))[1] +
                               1):(which(grepl("/pre", res2_split ))[1] - 1)]
    if(any(grepl("hr", res2_split))){
      res2_split = res2_split[1:(which(grepl("<hr>", res2_split))[1] - 1)]
    }
    res2_split <- strsplit(res2_split, " +")
    res2_split <- do.call(rbind, res2_split)
    res2_split <- as.data.frame(res2_split, stringsAsFactors = F)
    colnames(res2_split) <- c("id", "Cmax", "Cmax.pos", "Ymax",
                              "Ymax.pos", "Smax", "Smax.pos", "Smean", "Dmean",
                              "is.sp", "Dmaxcut", "Networks.used")
    collected_res[[i]] <- res2_split
  }
  collected_res <- do.call(rbind, collected_res)
  collected_res$is.signalp <- collected_res$is.sp == "Y"
  return(collected_res)
}
