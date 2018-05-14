#' Predict hydroxyproline positions in plant proteins based on primary structure
#'
#' predict_hyp is a hydroxyproline site prediction algorithm for plant proteins, based on the xgboost distributed gradient boosting library.
#' It was trained on plant sequences with experimentally determined 4-hydroxyprolines from uniprot data base. Prediction is not possible for prolines which are within 10 N-terminal and 6 C-terminal amino acids, they will be excluded from output.
#'
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class "SeqFastaAA" resulting from seqinr::read.fasta call.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param tprob A numeric value indicating the threshold for prediction. Acceptable values are in 0 - 1 range. At default set to 0.3 offering a tradeoff between sensitivity and specificity.
#' @param split A numeric value determining the ratio of vectorized and sequential computation. Should be left at default, lower to 0 - 1 range if low memory errors occur. Increase at your own risk.
#' @return A list with two elements:
#' \describe{
#'   \item{prediction}{data frame with columns:
#'   id - character, indicating the inputted protein id's;
#'   substr - character, indicating the sequence substring which was used for predictions;
#'   P_pos - integer, position of proline in the sequence;
#'   prob - numeric, predicted probability of being hydroxyproline;
#'   HYP - character, is the site predicted as a hydroxyproline}
#'   \item{sequence}{data frame with columns:
#'   sequence - sequences with prolines - P substituted with hydroxyprolines - O according to the prediction;
#'   id - corresponding id's}
#' }
#'
#' 
#' @examples
#' library(ragp)
#' data(at_nsp)
#'
#' #ramdom indexes
#' ind <- c(129, 145, 147, 160, 170,
#'     180, 189, 203, 205, 214, 217, 224)
#'
#' hyp_pred <- predict_hyp(sequence = at_nsp$sequence[ind],
#'                         id = at_nsp$Transcript.id[ind])
#' @export



predict_hyp <- function (data = NULL, sequence, id, tprob = 0.3, split = 1){
  if (missing(tprob)) {
    tprob <- 0.3
  }
  if (length(tprob) > 1){
    tprob <- 0.3
    warning("tprob should be of length 1, setting to default: tprob = 0.3")
  }
  if (!is.numeric(tprob)) {
    tprob <- as.numeric(tprob)
    warning("tprob is not numeric, converting using 'as.numeric'")
  }
  if (is.na(tprob)){
    tprob <- 0.3
    warning("tprob was set to NA, setting to default: tprob = 0.3")
  }
  if (tprob < 0) {
    tprob <- 0.3
    warning(paste("treshold probability for prediction should be in 0 - 1 range,", 
                  "tprob was set to the default 0.3"))
  }
  if (tprob > 1) {
    tprob <- 0.3
    warning(paste("treshold probability for prediction should be in 0 - 1 range,", 
                  "tprob was set to the default 0.3"))
  }
  if (missing(split)) {
    split <- 1
  }
  if (length(split) > 1){
    split <- 1
    warning("split should be of length 1, setting to default: split = 1")
  }
  if (!is.numeric(split)) {
    split <- as.numeric(split)
    warning("split is not numeric, converting using 'as.numeric'")
  }
  if (is.na(split)){
    split <- 1
    warning("split was set to NA, setting to default: split = 1")
  }
  if (split > 2){
    warning("setting split > 2 can cause memory problems")
  }
  if (split < 0.2){
    warning("setting split < 0.2 can increase computation time")
  }
  if(missing(data)){
    if (missing(sequence)){
      stop("protein sequence must be provided to obtain predictions")
    }
    if (missing(id)){
      stop("protein id must be provided to obtain predictions")
    }
    id <- as.character(id)
    sequence <- toupper(as.character(sequence))
    if (length(sequence) != length(id)){
      stop("id and sequence vectors are not of same length")
    }
  }
  if(class(data[[1]]) ==  "SeqFastaAA"){
    dat <- lapply(data, paste0, collapse ="")
    id <- names(dat)
    sequence <- toupper(as.character(unlist(dat)))
  }
  if(class(data) == "data.frame"){
    if(missing(sequence)){
      stop("the column name with the sequences must be specified")
    }
    if(missing(id)){
      stop("the column name with the sequence id's must be specified")
    }
    id <- as.character(substitute(id))
    sequence <- as.character(substitute(sequence))
    if (length(id) != 1L){
      stop("only one column name for 'id' must be specifed")
    }
    if (length(sequence) != 1L){
      stop("only one column name for 'sequence' must be specifed")
    }
    id <- if(id %in% colnames(data)){
      data[[id]]
    } else {
      stop("specified 'id' not found in data")
    }
    id <- as.character(id)  
    sequence  <- if(sequence %in% colnames(data)){
      data[[sequence]]
    } else {
      stop("specified 'sequence' not found in data")
    }
    sequence <- toupper(as.character(sequence))
  }
  if(class(data) == "character"){
    if (file.exists(data)){
      dat <- seqinr::read.fasta(file = data,
                                seqtype = "AA",
                                as.string = FALSE)
      dat <- lapply(dat, paste0, collapse ="")
      id <- names(dat)
      sequence <- toupper(as.character(unlist(dat)))
    } else {
      stop("cannot find file in the specified path")
    }
  }
  sequence <- sub("\\*$", "", sequence)
  if (sum(grepl("P", sequence)) == 0){
    warning("no prolines in the target sequences")
    return(NULL)
  }
  
  splt <- split * 10000
  pam <- ((seq(length(sequence)) - 1)%/%splt) + 1
  m_split <- split(data.frame(sequence, id), pam)
  predict <- utils::getFromNamespace("predict.xgb.Booster", "xgboost")
  result <- lapply(m_split, function(x) {
    sequence <- x[, 1]
    id <- x[, 2]
    seq_kmer <- vector("list", length(sequence))
    for (i in 1:length(sequence)) {
      seq_kmer[[i]] <- getKmer(sequence = sequence[i], 
                               id = id[i], kmer = 10)
    }
    seq_kmer <- do.call(rbind, seq_kmer)
    seq_kmer[, 2] <- as.character(seq_kmer[, 2])
    
    seq_kmer <- cbind(seq_kmer, order = 1:nrow(seq_kmer))
    non_aa <- grep(paste("[^",
                         paste(AADict,
                               sep = "", collapse = ""), 
                         "]",
                         sep = "",
                         collapse = ""),
                   seq_kmer[,2])
    if (sum((grep(paste("[^", paste(AADict, sep = "", collapse = ""), 
                        "]", sep = "", collapse = ""), seq_kmer[,2]))) != 
        0) {
      warning("characters other than single letter code for amino acids are present")
    }
    if(length(non_aa) > 0){
      seq_non_aa <- seq_kmer[non_aa,, drop = FALSE]
      seq_kmer <- seq_kmer[-non_aa,]
    } else {
      seq_non_aa <- NULL
    }
    
    kmer21 <- seq_kmer[nchar(seq_kmer[,2]) == 21,]
    kmer13 <- seq_kmer[nchar(seq_kmer[,2]) < 21 &
                         nchar(seq_kmer[,2]) >= 17 &
                         as.numeric(seq_kmer[,3]) > 10, ,
                       drop = FALSE]
    
    if(length(kmer13) != 0){
      kmer13[,2] <- substring(kmer13[,2], first = 5, last = 17)
      
      MBI_13 <- getAAindex(kmer13[,2],  aaidx)
      MoreauBroto_lag6_13 <- extractMBdesc(kmer13[,2])
      QSO_lag12_13 <- QSOlevel(kmer13[,2])
      ATC_13 <- getAAindex(kmer13[,2], Atchley)
      CTDC_13 <- do.call(rbind, lapply(kmer13[,2], function(x) CTDC(x)))
      
      dtest13 <- data.matrix(cbind(ATC_13,
                                   MoreauBroto_lag6_13,
                                   QSO_lag12_13 ,
                                   MBI_13,
                                   CTDC_13))
      
      
      prob_13 <- predict(model_13, dtest13)
      HYP_13 <- ifelse(prob_13 >= tprob, "Yes", "No")
      prediction13 <- cbind(id = as.character(kmer13[,1]),
                            substr = as.character(kmer13[,2]), 
                            P_pos = as.character(kmer13[,3]),
                            prob = as.character(prob_13), 
                            HYP = as.character(HYP_13),
                            order = kmer13[,5])
    } else {prediction13 <- NULL }
    
    MBI_21 <- getAAindex(kmer21[,2], aaidx)
    MoreauBroto_lag6_21 <- extractMBdesc(kmer21[,2])
    QSO_lag12_21 <- QSOlevel(kmer21[,2])
    ATC_21 <- getAAindex(kmer21[,2], Atchley)
    CTDC_21 <- do.call(rbind, lapply(kmer21[,2], function(x) CTDC(x)))
    
    dtest21 <- data.matrix(cbind(ATC_21,
                                 MoreauBroto_lag6_21,
                                 QSO_lag12_21 ,
                                 MBI_21,
                                 CTDC_21))
    
    
    prob_21 <- predict(model_21, dtest21)
    
    HYP_21 <- ifelse(prob_21 >= tprob, "Yes", "No")
    
    prediction21 <- cbind(id = as.character(kmer21[,1]),
                          substr = as.character(kmer21[,2]), 
                          P_pos = as.character(kmer21[,3]),
                          prob = as.character(prob_21), 
                          HYP = as.character(HYP_21),
                          order = kmer21[,5])
    
    res <- as.data.frame(rbind(prediction21,
                               prediction13),
                         stringsAsFactors = FALSE)
    res[,3] <- as.numeric(res[,3])
    res[,4] <- as.numeric(res[,4])
    if (length(non_aa) > 0){
      prediction_none <- cbind(id = as.character(seq_non_aa[,1]),
                               substr = as.character(seq_non_aa[,2]), 
                               P_pos = as.character(seq_non_aa[,3]),
                               prob = NA, 
                               HYP = NA,
                               order = seq_non_aa[,5])
      res <- rbind(res, prediction_none)
      res <- res[order(as.numeric(res[,6])),]
      res <- res[,-6]
    } else {
      res <- res[order(as.numeric(res[,6])),]
      res <- res[,-6]
    }
    return(res)
  })
  prediction <- as.data.frame(do.call(rbind, result), stringsAsFactors = FALSE)
  prediction$prob <- as.numeric(prediction$prob)
  prediction$P_pos <- as.integer(prediction$P_pos)
  row.names(prediction) <- 1:nrow(prediction)
  seq <- sub_hyp(sequence, id, prediction)
  seq <- data.frame(sequence = seq, id = id, stringsAsFactors = FALSE)
  result <- list(prediction = prediction, sequence = seq)
  return(result)
}