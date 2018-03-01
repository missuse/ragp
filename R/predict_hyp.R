#' Predict hydroxyproline positions in plant proteins based on primary structure
#'
#' predict_hyp is a hydroxyproline site prediction algorithm for plant proteins, based on the xgboost distributed gradient boosting library.
#' It was trained on plant sequences with experimentally determined 4-hydroxyprolines from uniprot data base. N- and C- terminal prolines surrounded by less than 10 amino acid residues will be excluded.
#'
#' @param sequence A vector of strings representing protein amino acid sequences
#' @param id A vector of strings representing protein id's
#' @param tprob A numeric value indicating the treshold for prediction. Acceptable values are in 0 - 1 range. At default set to 0.32 offering a tradeoff between sensitivity and specificity.
#' @param split A numeric value determining the ratio of vectorized and parallelized computation. Should be left at default, lower to 0 - 1 range if low memory errors occur. Increase at your own risk.
#' @return  A list with two elements:
#' \describe{
#'   \item{prediction}{data frame with columns:
#'   id - character, indicating the inputed protein id's;
#'   substr - character, indicating the sequence substring which was used for predictions;
#'   P_pos - integer, position of proline in the sequence;
#'   prob - numeric, predicted probability of beeing hydroxyproline}
#'   \item{sequence}{sequences with prolines - P substituted with hydroxyprolines - O according to the prediction}
#' }
#'
#' @examples
#' library(ragp)
#' data(at_nsp)
#'
#' #ramdom indexes
#' ind <- c(129, 145, 147, 160, 170,
#'     180, 189, 203, 205, 214, 217, 224)
#'
#' hyp_pred <- predict_hyp(sequence = at_nsp$sequence[ind], id = at_nsp$Transcript.id[ind])
#' @export


predict_hyp <- function (sequence, id, tprob = 0.32, split = 1) 
{
  if (missing(sequence)) {
    stop("protein sequence must be provided to obtain predictions")
  }
  if (missing(id)) {
    stop("protein id must be provided to obtain predictions")
  }
  if (missing(tprob)) {
    tprob <- 0.32
  }
  if (tprob < 0) {
    tprob <- 0.32
    warning(paste("treshold probability for prediction should be in 0 - 1 range,", 
                  "tprob was set to the default 0.32"))
  }
  if (tprob > 1) {
    tprob <- 0.32
    warning(paste("treshold probability for prediction should be in 0 - 1 range,", 
                  "tprob was set to the default 0.32"))
  }
  if (sum(grepl("P", sequence)) == 0){
    warning("no prolines in the target sequences")
    return(NULL)
  }
  
  splt <- split * 10000
  sequence <- as.character(sequence)
  id <- as.character(id)
  if (length(sequence) != length(id)) 
    stop("id and sequence vectors are not of same length")
  
  pam <- ((seq(length(sequence)) - 1)%/%splt) + 1
  m_split <- split(data.frame(sequence, id), pam)
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
    
    if (sum((grepl(paste("[^", paste(AADict, sep = "", collapse = ""), 
                         "]", sep = "", collapse = ""), seq_kmer[,2]))) != 
        0) {
      stop("characters other than single letter code for amino acids are present")
    }
    
    kmer21 <- seq_kmer[nchar(seq_kmer[,2]) == 21,]
    kmer13 <- seq_kmer[nchar(seq_kmer[,2]) < 21 & nchar(seq_kmer[,2]) >= 17 & as.numeric(seq_kmer[,3]) > 10, , drop = FALSE]
    
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
      
      
      prob_13 <- xgboost:::predict.xgb.Booster(model_13, dtest13)
      HYP_13 <- ifelse(prob_13 >= tprob, "Yes", "No")
      prediction13 <- cbind(id = as.character(kmer13[,1]),
                            substr = as.character(kmer13[,2]), 
                            P_pos = as.character(kmer13[,3]),
                            prob = as.character(prob_13), 
                            HYP = as.character(HYP_13))
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
    
    
    prob_21 <- xgboost:::predict.xgb.Booster(model_21, dtest21)
    
    HYP_21 <- ifelse(prob_21 >= tprob, "Yes", "No")
    
    prediction21 <- cbind(id = as.character(kmer21[,1]),
                          substr = as.character(kmer21[,2]), 
                          P_pos = as.character(kmer21[,3]),
                          prob = as.character(prob_21), 
                          HYP = as.character(HYP_21))
    
    res <- as.data.frame(rbind(prediction21, prediction13), stringsAsFactors = FALSE)
    res[,3] <- as.numeric(res[,3])
    res[,4] <- as.numeric(res[,4])
    res <- res[order(res[,1], res[,3]),]
    res <- merge(data.frame(id = factor(id, levels = id)), res)
    res <- res[order(res$id, res$P_pos),]
    res <- as.matrix(res)
    return(res)
  })
  prediction <- as.data.frame(do.call(rbind, result), stringsAsFactors = FALSE)
  prediction$prob <- as.numeric(prediction$prob)
  prediction$P_pos <- as.integer(prediction$P_pos)
  row.names(prediction) <- 1:nrow(prediction)
  seq <- sub_hyp(sequence, id, prediction)
  result <- list(prediction = prediction, sequence = seq)
  return(result)
}