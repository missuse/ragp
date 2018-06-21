#' Plotting SignalP prediction.
#' 
#' Plots the SignalP prediction for one protein sequence using base graphics. SignalP 4.1 server predicts the presence and location of signal peptide cleavage sites in amino acid sequences from different organisms: Gram-positive prokaryotes, Gram-negative prokaryotes, and eukaryotes. The method incorporates a prediction of cleavage sites and a signal peptide/non-signal peptide prediction based on a combination of several artificial neural networks.
#'
#' @param sequence String representing a protein amino acid sequence.
#' @param id String representing a protein identifier.
#' @param org_type One of c("euk", "gram-", "gram+"), defaults to "euk". Which model should be used for prediction.
#' @param Dcut_type One of c("default", "sensitive", "user"), defaults to "default". The default cutoff values for SignalP 4 are chosen to optimize the performance measured as Matthews Correlation Coefficient (MCC). This results in a lower sensitivity (true positive rate) than SignalP 3.0 had. Setting this argument to "sensitive" will yield the same sensitivity as SignalP 3.0. This will make the false positive rate slightly higher, but still better than that of SignalP 3.0.
#' @param Dcut_noTM A numeric value, with range 0 - 1, defaults to 0.45. For experimenting with cutoff values.
#' @param Dcut_TM A numeric value, with range 0 - 1, defaults to 0.5. For experimenting with cutoff values.
#' @param method One of c("best", "notm"), defaults to "best". Signalp 4.1 contains two types of neural networks. SignalP-TM has been trained with sequences containing transmembrane segments in the data set, while SignalP-noTM has been trained without those sequences. Per default, SignalP 4.1 uses SignalP-TM as a preprocessor to determine whether to use SignalP-TM or SignalP-noTM in the final prediction (if 4 or more positions are predicted to be in a transmembrane state, SignalP-TM is used, otherwise SignalP-noTM). An exception is Gram-positive bacteria, where SignalP-TM is used always. If you are confident that there are no transmembrane segments in your data, you can get a slightly better performance by choosing "Input sequences do not include TM regions", which will tell SignalP 4.1 to use SignalP-noTM always.
#' @param c.score.col Plotting color of the C-score line. At default set to: '#ff0000'.
#' @param s.score.col Plotting color of the S-score line. At default set to: '#728fcc'.
#' @param y.score.col Plotting color of the Y-score line. At default set to: '#728fcc'.
#' @param t.col Plotting color of the threshold line. At default set to: '#551a8b'.
#' @param main Title of the plot.
#' @param sleep A numeric indicating the pause in seconds between POST and GET server calls, at default set to 5s. Decreasing is not recommended.
#' @return A list with two elements:
#' \describe{
#'   \item{prediction}{Data frame with the prediction results.}
#'   \item{plot}{Data frame with values used for plotting.}
#' }
#'
#' @source \url{http://www.cbs.dtu.dk/services/SignalP-4.1/}
#' @references Petersen TN. Brunak S. Heijne G. Nielsen H. (2011) SignalP 4.0: discriminating signal peptides from transmembrane regions. Nature Methods 8: 785-786
#'
#' @seealso \code{\link[ragp]{get_signalp}}
#'
#' @examples
#' library(ragp)
#' pred <- plot_signalp(sequence = at_nsp$sequence[5],
#'                      id = at_nsp$Transcript.id[5])
#' 
#' @export

plot_signalp <- function(sequence,
                         id,
                         org_type = c("euk", "gram-", "gram+"),
                         Dcut_type = c("default", "sensitive", "user"),
                         Dcut_noTM = 0.45,
                         Dcut_TM = 0.5,
                         method = c("best", "notm"),
                         c.score.col = "#ff0000",
                         s.score.col = "#59a454",
                         y.score.col = "#728fcc",
                         t.col = "#551a8b",
                         main = NULL,
                         sleep = 5L){
  if (missing(sleep)) {
    sleep <- 5
  }
  if (length(sleep) > 1){
    sleep <- 5
    warning("sleep should be of length 1, setting to default: sleep = 5")
  }
  if (!is.numeric(sleep)){
    sleep <- as.numeric(sleep)
    warning("sleep is not numeric, converting using 'as.numeric'")
  }
  if (is.na(sleep)){
    sleep <- 5
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
  if (missing(sequence)){
    stop("protein sequence must be provided to obtain predictions")
  }
  if (missing(id)){
    stop("protein id must be provided to obtain predictions")
  }
  if (length(sequence) != 1){
    stop("one string representing a protein sequence should be provided")
  }
  if (length(id) != 1){
    stop("one string representing a protein id should be provided")
  }
  areColors <- function(x) {
    sapply(x, function(X) {
      tryCatch(is.matrix(grDevices::col2rgb(X)), 
               error = function(e) FALSE)
    })
  }
  if (length(c.score.col) > 1){
    c.score.col <- '#ff0000'
    warning("One color should be provided for c.score.col. Using default: '#ff0000'")
  }
  if (!areColors(c.score.col)){
    c.score.col <- '#ff0000'
    warning("c.score.col provided is not a valid color, default will be used: '#ff0000'")
  }
  if (length(s.score.col) > 1){
    s.score.col <- '#59a454'
    warning("One color should be provided for s.score.col. Using default: '#59a454'")
  }
  if (!areColors(s.score.col)){
    s.score.col <- '#59a454'
    warning("s.score.col provided is not a valid color, default will be used: '#59a454'")
  }
  if (length(y.score.col) > 1){
    y.score.col <- '#728fcc'
    warning("One color should be provided for y.score.col. Using default: '#728fcc'")
  }
  if (!areColors(y.score.col)){
    y.score.col <- '#728fcc'
    warning("y.score.col provided is not a valid color, default will be used: '#728fcc'")
  }
  if (length(t.col) > 1){
    t.col <- '#551a8b'
    warning("One color should be provided for t.col. Using default: '#551a8b'")
  }
  if (!areColors(t.col)){
    t.col <- '#551a8b'
    warning("t.col provided is not a valid color, default will be used: '#551a8b'")
  }
  sequence <- sub("\\*$", "", sequence)
  sequence <- toupper(as.character(sequence))
  id <- as.character(id)
  res <- httr::POST(
    url = "http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi",
    encode = "form",
    body = list(
      `configfile` = "/usr/opt/www/pub/CBS/services/SignalP-4.1/SignalP.cf",
      `SEQPASTE` =  sequence,
      `orgtype` = org_type,
      `Dcut-type` = Dcut_type,
      `Dcut-noTM` = Dcut_noTM,
      `Dcut-TM` = Dcut_TM,
      `graphmode` = NULL,
      `format` = "long",
      `minlen` = "",
      `method` = method,
      `trunc` = ""
    ))
  res <- httr::content(res,
                       as = "parsed")
  res <- xml2::xml_find_all(res,
                            ".//input[@name='jobid']")
  jobid <- xml2::xml_attr(res,
                             "value")
  Sys.sleep(sleep)
  repeat {
    res2 <- httr::GET(
      url = "http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi",
      query = list(
        jobid = jobid,
        wait = "20"
        ))
    
    res2 <- as.character(
      xml2::xml_find_all(
        httr::content(res2,
                      as = "parsed"),
        ".//pre")
    )
    
    res2 <- unlist(strsplit(res2,
                            "\n"))
    if (any(grepl("Length", res2))) {
      break
    }
  }

  tit <- res2[grep("SignalP-4.1",
                   res2)[1]]
  tit <- gsub("#", "", tit)
  tit <- trimws(tit)
  isnoTM <- grepl("noTM",
                  res2[grep("pos",
                            res2)[1]-1])
  if(isnoTM) {
    subtit <- paste("Networks: SignalP-noTM")
  } else {
    subtit <- paste("Networks: SignalP-TM")
  }
  res2_tab <- res2[(grep("pos",
                         res2)[1]+1):(grep("Measure",
                                           res2)[1]-1)]
  res2_tab <- utils::read.table(text = res2_tab,
                                stringsAsFactors = FALSE)
  colnames(res2_tab) <- c("pos",
                          "aa",
                          "C",
                          "S",
                          "Y")
  res2_out <- res2[(grep("Measure",
                         res2)[1]+1):(grep("Measure",
                                           res2)[1]+4)]
  
  res2_out <- utils::read.table(text = res2_out,
                                stringsAsFactors = FALSE)
  last_line <- unlist(strsplit(res2[(grep("Measure",
                                          res2)[1]+5)], " +"))
  res2_out <- rbind(res2_out[,2:4],
                    last_line[2:4])
  
  colnames(res2_out) <- c("Measure",
                          "Position",
                          "value")
  res2_out$Measure <- c("max.C",
                        "max.Y",
                        "max.S",
                        "mean.S",
                        "D")
  res2_out$value <- as.numeric(as.character(res2_out$value))
  res2_out$Cutoff <- c(rep("", 4),
                       rev(last_line)[2])
  out <- list(prediction = res2_out,
              plot = res2_tab)
  graphics::plot(res2_tab$S,
                 ylim = c(-0.1, 1),
                 type = "l",
                 col = s.score.col,
                 xlab = "Position",
                 ylab = "Score",
                 yaxt = "n")
  
  if (missing(main)){
    maint <- paste(tit,
                 "\n",
                 subtit,
                 "\n",
                 "id: ",
                 id[1],
                 sep = "")
  } else {
    maint <- main
  }
  graphics::title(main = maint,
                  adj = 0,
                  cex.main = 1)
  graphics::axis(2, yaxp = c(0, 1, 5), las = 2)
  if(isnoTM) {
    graphics::abline(h = Dcut_noTM,
                     col = t.col,
                     lty = 2)
  } else {
    graphics::abline(h = Dcut_TM,
                     col = t.col,
                     lty = 2)
  }
  graphics::segments(x0 = res2_tab$pos,
                     y0 = res2_tab$C,
                     y1 = 0,
                     col = c.score.col)
  graphics::lines(res2_tab$Y,
                  col = y.score.col)
  graphics::legend("topright",
                   col = c(c.score.col,
                           s.score.col,
                           y.score.col),
                   legend = c("C-score",
                              "S-score",
                              "Y-score"),
                   lty = 1)
  max_pos <- max(as.numeric(res2_tab$pos))
  graphics::text(seq(1,
                     max_pos,
                     by = 1),
                 -0.02 ,
                 labels = res2_tab$aa,
                 pos = 1,
                 xpd = TRUE,
                 cex = 0.7)
  invisible(out)
}
