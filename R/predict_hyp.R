#' Predict hydroxyproline positions in plant proteins based on primary structure
#'
#' predict_hyp is a hydroxyproline site prediction algorithm for plant proteins, based on the xgboost distributed gradient boosting library.
#' It was trained on plant sequences with experimentally determined 4-hydroxyprolines from swissprot data base. Prediction is not possible for prolines which are within 10 N-terminal and 6 C-terminal amino acids (V1 model version) and 10 N-terminal and 7 C-terminal amino acids (V2 model version), they will be excluded from output.
#' 
#' @aliases predict_hyp predict_hyp.default predict_hyp.character predict_hyp.data.frame predict_hyp.list
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class "SeqFastaAA" resulting from \code{\link[seqinr]{read.fasta}} call. Should be left blank if vectors are provided to sequence and id arguments.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank. Ids should be unique.
#' @param tprob A numeric value indicating the threshold for prediction. Acceptable values are in 0 - 1 range. At default set to 0.3 for "V1" model and 0.224 for "V2" model, offering a tradeoff between sensitivity and specificity.
#' @param version A string indicating which model version to use: the first version "V1", or the second version "V2". Default is "V2".
#' @param split A numeric value determining the ratio of vectorized and sequential computation. Should be left at default, lower to 0 - 1 range if low memory errors occur. Increase at your own risk.
#' @param ... currently no additional arguments are accepted apart the ones documented bellow.
#'
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
#' @details Previously trained xgboost models were re-saved using xgboost 1.1.1.1 to increase compatibility. While using the mentioned xgboost version the returned predictions are equal to previous. However, using earlier xgboost versions with the new models will result in slightly different predicted probabilities.
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
#'
#' @import seqinr
#' @import xgboost
#' @import utils
#' @export

predict_hyp <- function (data, ...){
  if (missing(data) || is.null(data)) predict_hyp.default(...)
  else UseMethod("predict_hyp")
}


#' @rdname predict_hyp
#' @method predict_hyp character
#' @export

predict_hyp.character <- function(data,
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
  res <- predict_hyp.default(sequence = sequence,
                             id = id,
                             ...)
  return(res)
}

#' @rdname predict_hyp
#' @method predict_hyp data.frame
#' @export

predict_hyp.data.frame <- function(data,
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
  res <- predict_hyp.default(sequence = sequence,
                             id = id,
                             ...)

  return(res)
}

#' @rdname predict_hyp
#' @method predict_hyp list
#' @export

predict_hyp.list <- function(data,
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
  res <- predict_hyp.default(sequence = sequence,
                             id = id,
                             ...)
  return(res)
}

#' @rdname predict_hyp
#' @method predict_hyp default
#' @export

predict_hyp.default <- function (data = NULL,
                                 sequence,
                                 id,
                                 version = "V2",
                                 tprob = ifelse(version == "V1", 0.3, 0.224),
                                 split = 1,
                                 ...){
  if (missing(version)) {
    version = "V2"
  }
  if (!version %in% c("V1", "V2")){
    stop("incorect version specification, use 'V1' or 'V2'")
  }
  if (missing(tprob)) {
    tprob <- ifelse(version == "V1", 0.3, 0.224)
  }
  if (length(tprob) > 1){
    stop("tprob should be of length 1",
         call. = FALSE)
  }
  if (!is.numeric(tprob)) {
    tprob <- as.numeric(tprob)
    warning("tprob is not numeric, converting using 'as.numeric'",
            call. = FALSE)
  }
  if (is.na(tprob)){
    stop("tprob was set to NA",
         call. = FALSE)
  }
  if (tprob < 0) {
    stop("treshold probability for prediction should be in 0 - 1 range, ",
         call. = FALSE)
  }
  if (tprob > 1) {
    stop("treshold probability for prediction should be in 0 - 1 range, ",
         call. = FALSE)
  }
  if (missing(split)) {
    split <- 1
  }
  if (length(split) > 1){
    split <- 1
    warning("split should be of length 1, setting to default: split = 1",
            call. = FALSE)
  }
  if (!is.numeric(split)) {
    split <- as.numeric(split)
    warning("split is not numeric, converting using 'as.numeric'",
            call. = FALSE)
  }
  if (is.na(split)){
    split <- 1
    warning("split was set to NA, setting to default: split = 1",
            call. = FALSE)
  }
  if (split > 2){
    warning("setting split > 2 can cause memory problems",
            call. = FALSE)
  }
  if (split < 0.2){
    warning("setting split < 0.2 can increase computation time",
            call. = FALSE)
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
  if (length(unique(id)) != length(id)){
    stop("id's should be unique",
         call. = FALSE)
  }
  sequence <- sub("\\*$",
                  "",
                  sequence)
  if (sum(grepl("P", sequence)) == 0){
    warning("no prolines in the target sequences")
    return(NULL)
  }
  splt <- split * 10000
  pam <- ((seq(length(sequence)) - 1) %/% splt) + 1
  m_split <- split(data.frame(sequence,
                              id),
                   pam)
  predict <- utils::getFromNamespace("predict.xgb.Booster",
                                     "xgboost")
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
      warning("characters other than single letter code for amino acids are present",
              call. = FALSE)
    }
    if(length(non_aa) > 0){
      seq_non_aa <- seq_kmer[non_aa,, drop = FALSE]
      seq_kmer <- seq_kmer[-non_aa,]
    } else {
      seq_non_aa <- NULL
    }
    if(version == "V1"){
      kmer21 <- seq_kmer[nchar(seq_kmer[,2]) == 21,]
      kmer13 <- seq_kmer[nchar(seq_kmer[,2]) < 21 &
                           nchar(seq_kmer[,2]) >= 17 &
                           as.numeric(seq_kmer[,3]) > 10, ,
                         drop = FALSE]
      
      if(length(kmer13) != 0){
        kmer13[,2] <- substring(kmer13[,2], first = 5, last = 17)
        
        MBI_13 <- getAAindex(kmer13[,2],  aaidx)
        MoreauBroto_lag6_13 <- extractMBdesc(kmer13[,2], nlag = 6, aaidx)
        QSO_lag12_13 <- QSOlevel(kmer13[,2])
        ATC_13 <- getAAindex(kmer13[,2], Atchley)
        CTDC_13 <- do.call(rbind, lapply(kmer13[,2], function(x) CTDC(x)))
 
        dtest13 <- xgboost::xgb.DMatrix(data = cbind(ATC_13,
                                                     MoreauBroto_lag6_13,
                                                     QSO_lag12_13,
                                                     MBI_13,
                                                     CTDC_13))
        model_13 <- xgboost::xgb.load(system.file("extdata",
                                                  "model_kmer13.model",
                                                  package = "ragp"))
        
        prob_13 <- predict(model_13, dtest13)
        HYP_13 <- ifelse(prob_13 >= 0.3, "Yes", "No")
        prediction13 <- cbind(id = as.character(kmer13[,1]),
                              substr = as.character(kmer13[,2]),
                              P_pos = as.character(kmer13[,3]),
                              prob = as.character(prob_13),
                              HYP = as.character(HYP_13),
                              order = kmer13[,5])
        } else {
          prediction13 <- NULL
          }
      MBI_21 <- getAAindex(kmer21[,2], aaidx)
      MoreauBroto_lag6_21 <- extractMBdesc(kmer21[,2], nlag = 6, aaidx)
      QSO_lag12_21 <- QSOlevel(kmer21[,2])
      ATC_21 <- getAAindex(kmer21[,2], Atchley)
      CTDC_21 <- do.call(rbind, lapply(kmer21[,2],
                                       function(x) CTDC(x)))
      dtest21 <- xgboost::xgb.DMatrix(data = cbind(ATC_21,
                                                   MoreauBroto_lag6_21,
                                                   QSO_lag12_21,
                                                   MBI_21,
                                                   CTDC_21))
      model_21 <- xgboost::xgb.load(system.file("extdata",
                                                "model_kmer21.model",
                                                package = "ragp"))
        
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
    }
    if( version == "V2"){
      kmer21 <- seq_kmer[nchar(seq_kmer[,2]) == 21,]
      kmer15 <- seq_kmer[nchar(seq_kmer[,2]) < 21 &
                           nchar(seq_kmer[,2]) >= 18 &
                           as.numeric(seq_kmer[,3]) > 10, ,
                         drop = FALSE]
      if(length(kmer15) != 0){
        kmer15[,2] <- substring(kmer15[,2], first = 4, last = 18)
        
        MBI_15 <- getAAindex(kmer15[,2],  aaidx[c(2,4,5,3,1,6),])
        MoreauBroto_lag12_15_mbi <- extractMBdesc(kmer15[,2], nlag = 12, aaidx)
        MoreauBroto_lag12_15_atc <- extractMBdesc(kmer15[,2], nlag = 12, Atchley[,-21])
        QSO_lag12_15 <- QSOlevel(kmer15[,2])
        
        dtest15 <- xgboost::xgb.DMatrix(data = cbind(MBI_15,
                                                     QSO_lag12_15,
                                                     MoreauBroto_lag12_15_mbi,
                                                     MoreauBroto_lag12_15_atc))
        
        model15_v2 <- xgboost::xgb.load(system.file("extdata",
                                                    "model15_v2.model",
                                                    package = "ragp"))
        
        prob_15 <- predict(model15_v2, dtest15)
        HYP_15 <- ifelse(prob_15 >= 0.22, "Yes", "No")
        prediction15 <- cbind(id = as.character(kmer15[,1]),
                              substr = as.character(kmer15[,2]),
                              P_pos = as.character(kmer15[,3]),
                              prob = as.character(prob_15),
                              HYP = as.character(HYP_15),
                              order = kmer15[,5])
      } else {
        prediction15 <- NULL
      }
      MBI_21 <- getAAindex(kmer21[,2], aaidx[c(2,4,5,3,1,6),])
      MoreauBroto_lag12_21_mbi <- extractMBdesc(kmer21[,2], nlag = 12, aaidx)
      MoreauBroto_lag12_21_atc <- extractMBdesc(kmer21[,2], nlag = 12, Atchley[,-21])
      QSO_lag12_21 <- QSOlevel(kmer21[,2])
      dtest21 <- xgboost::xgb.DMatrix(data = cbind(MBI_21,
                                                   QSO_lag12_21,
                                                   MoreauBroto_lag12_21_mbi,
                                                   MoreauBroto_lag12_21_atc))
      model21_v2 <- xgboost::xgb.load(system.file("extdata",
                                                  "model21_v2.model",
                                                  package = "ragp"))
      prob_21 <- predict(model21_v2, dtest21)
      HYP_21 <- ifelse(prob_21 >= tprob, "Yes", "No")
      prediction21 <- cbind(id = as.character(kmer21[,1]),
                            substr = as.character(kmer21[,2]),
                            P_pos = as.character(kmer21[,3]),
                            prob = as.character(prob_21),
                            HYP = as.character(HYP_21),
                            order = kmer21[,5])
      res <- as.data.frame(rbind(prediction21,
                                 prediction15),
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
    }
  }
  )
  prediction <- as.data.frame(do.call(rbind, result),
                              stringsAsFactors = FALSE)
  prediction$prob <- as.numeric(prediction$prob)
  prediction$P_pos <- as.integer(prediction$P_pos)
  row.names(prediction) <- 1:nrow(prediction)
  seq <- sub_hyp(sequence,
                 id,
                 prediction)
  seq <- data.frame(sequence = seq,
                    id = id,
                    stringsAsFactors = FALSE)
  result <- list(prediction = prediction,
                 sequence = seq)
  return(result)
}
