#' Predict hydroxyproline positions in plant proteins based on primary structure
#'
#' predict_hyp is a hydroxyproline site prediction algorithm for plant proteins, based on the xgboost distributed gradient boosting library.
#' It was trained on plant sequences with experimentally determined 4-hydroxyprolines from uniprot data base. N- and C- terminal prolines surrounded by less than 10 amino acid residues will be excluded.
#'
#' @param sequence A vector of strings representing protein amino acid sequences
#' @param id A vector of strings representing protein id's
#' @param tprob A numeric value indicating the treshold for prediction. Acceptable values are in 0 - 1 range. At default set to 0.33 which corresponds to 0.9674 sensitivity and 0.9568 specificity in 5-fold crossvalidation.
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


predict_hyp <- function(sequence, id, tprob = 0.33, split = 1){
  if(missing(sequence)){
    stop("protein sequence must be provided to obtain predictions")
  }
  
  if(missing(id)){
    stop("protein id must be provided to obtain predictions")
  }
  
  if(missing(tprob)){
    tprob <- 0.33
  }
  
  if(tprob < 0){
    tprob <- 0.33
    warning(paste("treshold probability for prediction should be in 0 - 1 range,",
                  "tprob was set to the default 0.33"))
  }
  
  if(tprob > 1){
    tprob <- 0.33
    warning(paste("treshold probability for prediction should be in 0 - 1 range,",
                  "tprob was set to the default 0.33"))
  }
  
  
  splt <- split * 10000
  
  sequence <- as.character(sequence)
  id <- as.character(id)
  if(length(sequence) != length(id)) stop("id and sequence vectors are not of same length")
  
  extractMBdesc <- function(x){
    xSplitted <- strsplit(as.character(x), split = "")
    P <- lapply(xSplitted,
                FUN = function(z) sapply(z, function(x) Pr[,colnames(Pr) == x]))
    P <- do.call(rbind, P)
    nlag <- 6
    N <- 21
    MB <- vector("list", 12)
    for (j in 1:nlag){
      MB[[j]] <- rowSums(P[,1:(N -j)]* P[,1:(N - j) + j])/(N - j)
    }
    MB <- do.call(rbind, MB)
    MB <- matrix(as.vector(MB), ncol = nlag*6, byrow=T)
    colnames(MB) = as.vector(t(outer(props, paste(".lag", 1:nlag,
                                                  sep = ""), paste, sep = "")))
    return(MB)
  }
  
  QSOlevel <- function(m){
    QSObabe <- function(x){
      w <- 0.1
      nlag <- 12
      N <- 21
      xSplitted <- strsplit(as.character(x), split = "")
      xSplitted <- do.call(rbind, xSplitted)
      tau1 <- list()
      for (d in 1:nlag) {
        for (i in 1:(N - d)){
          tau1[[length(tau1)+1]] <- diag(as.matrix(DistMat1)[xSplitted[,i],
                                                             xSplitted[,i + d]])^2
        }
      }
      tau1 <- t(do.call(rbind, tau1))
      tau1 <- cbind(rowSums(tau1[,1:20]),
                    rowSums(tau1[,21:39]),
                    rowSums(tau1[,40:57]),
                    rowSums(tau1[,58:74]),
                    rowSums(tau1[,75:90]),
                    rowSums(tau1[,91:105]),
                    rowSums(tau1[,106:119]),
                    rowSums(tau1[,120:132]),
                    rowSums(tau1[,133:144]),
                    rowSums(tau1[,145:155]),
                    rowSums(tau1[,156:165]),
                    rowSums(tau1[,166:174]))
      fr <- sapply(AADict,
                   function(y) sapply(gregexpr(y,
                                               as.character(x)),
                                      function(x) sum(x>0)))
      Xr1 <- fr/(1 + (w * rowSums(tau1)))
      colnames(Xr1)  <-  paste("Schneider.Xr.", colnames(Xr1), sep = "")
      Xd1 <- (w * tau1)/(1 + (w * rowSums(tau1)))
      colnames(Xd1) <- paste("Schneider.Xd.", 1:nlag, sep = "")
      tau2 <- list()
      for (d in 1:nlag) {
        for (i in 1:(N - d)){
          tau2[[length(tau2)+1]] <- diag(as.matrix(DistMat2)[xSplitted[,i],
                                                             xSplitted[,i + d]])^2
        }
      }
      tau2 <- t(do.call(rbind, tau2))
      tau2 <- cbind(rowSums(tau2[,1:20]),
                    rowSums(tau2[,21:39]),
                    rowSums(tau2[,40:57]),
                    rowSums(tau2[,58:74]),
                    rowSums(tau2[,75:90]),
                    rowSums(tau2[,91:105]),
                    rowSums(tau2[,106:119]),
                    rowSums(tau2[,120:132]),
                    rowSums(tau2[,133:144]),
                    rowSums(tau2[,145:155]),
                    rowSums(tau2[,156:165]),
                    rowSums(tau2[,166:174]))
      Xr2 <- fr/(1 + (w * rowSums(tau2)))
      colnames(Xr2) <- paste("Grantham.Xr.", colnames(Xr2), sep = "")
      Xd2 <- (w * tau2)/(1 + (w * rowSums(tau2)))
      colnames(Xd2) <- paste("Grantham.Xd.", 1:nlag, sep = "")
      QSO <- cbind(Xr1, Xr2, Xd1, Xd2)
      return(QSO)
    }
    splt <- 180
    k <- length(m)/splt
    if(k < 1){
      QSO_all <- QSObabe(m)
    }else{
      pam <- ((seq(length(m))-1) %/% splt) + 1
      m_split <- split(m, pam)
      QSO_all <- lapply(m_split, function(x) QSObabe(x))
      QSO_all <- do.call(rbind, QSO_all)
    }
    return(QSO_all)
  }
  
  getAAindex <- function(test_fragment){
    test_fragment <- as.character(test_fragment)
    test_fragment <- toupper(test_fragment)
    stopifnot(nchar(test_fragment) == 21)
    frag_vector <- unlist(strsplit(test_fragment, ""))
    aaidx <- aaidx[c(2,4,5,3,1,6),]
    frag_var <- sapply(frag_vector,
                       function(x) aaidx[,colnames(aaidx) == x])
    frag_var <- frag_var[,-11]
    frag_var <- as.numeric(frag_var)
    var_names <- expand.grid(as.character(1:nrow(aaidx)),
                             as.character((-10:10)[-11]))
    names(frag_var) <- paste(var_names$Var1, var_names$Var2, sep="_")
    return(frag_var)
  }
  
  getKmer <- function(sequence, id, kmer){
    sequence <- as.character(sequence)
    sequence <- toupper(sequence)
    n_char <- nchar(sequence)
    P_pos <- unlist(gregexpr("P", sequence))
    P_mer <- sapply(P_pos, function(x) substr(sequence,
                                              start = x - kmer,
                                              stop = x + kmer))
    P_mer <- cbind(id = as.character(rep(id, length(P_mer))),
                   substr = as.character(P_mer),
                   pos = as.character(P_pos),
                   nchar = as.character(rep(n_char, length(P_pos))))
    return(P_mer)
  }
  
  sub_hyp <- function(sequence, id, hyp){
    sequence <- unlist(lapply(id, function(x){
      hyp <- hyp[hyp$id == x,]
      p_pos <- as.numeric(as.character(hyp$P_pos[hyp$HYP == "Yes"]))
      sequencei <- as.character(sequence[id == x])
      for(i in p_pos){
        substr(sequencei, start = i, stop = i) <- "O"
      }
      sequencei
    }
    )
    )
  }
  
  k <- length(sequence)/splt
  
  pam <- ((seq(length(sequence))-1) %/% splt) + 1
  m_split <- split(data.frame(sequence, id), pam)
  
  result <- lapply(m_split, function(x){
    sequence <- x[,1]
    id <- x[,2]
    seq_kmer <- vector("list",length(sequence))
    for(i in 1:length(sequence)){
      seq_kmer[[i]] <- getKmer(sequence = sequence[i],
                               id = id[i],
                               kmer = 10)
    }
    
    seq_kmer <- do.call(rbind, seq_kmer)
    seq_kmer <- seq_kmer[nchar(as.character(seq_kmer[,2])) == 21,]
    seq_kmer[,2] <- as.character(seq_kmer[,2])
    seq_kmer <- as.data.frame(seq_kmer)
    
    if(sum((grepl(paste("[^", paste(AADict, sep = "",
                                    collapse = ""),"]",
                        sep = "", collapse = ""),
                  seq_kmer$substr))) != 0){
      stop("characters other than single letter code for amino acids are present")}
    
    MBI_var <- t(sapply(seq_kmer$substr,  FUN = getAAindex))
    row.names(MBI_var) <- 1:length(seq_kmer$substr)
    
    QSO_lag12_var <- QSOlevel(seq_kmer$substr)
    MoreauBroto_lag6_var <- extractMBdesc(seq_kmer$substr)
    
    dtest <- data.matrix(cbind(MBI_var,
                               MoreauBroto_lag6_var,
                               QSO_lag12_var))
    
    prob <- xgboost:::predict.xgb.Booster(model, dtest)
    
    HYP <- ifelse(prob >= tprob, "Yes", "No")
    
    prediction <- cbind(id = as.character(seq_kmer$id),
                        substr = as.character(seq_kmer$substr),
                        P_pos = as.character(seq_kmer$pos),
                        prob = as.character(prob),
                        HYP = as.character(HYP))
    return(prediction)
  }
  )
  prediction <- as.data.frame(do.call(rbind, result), stringsAsFactors = FALSE)
  prediction$prob <- as.numeric(prediction$prob )
  prediction$P_pos <- as.integer(prediction$P_pos)
  
  seq <- sub_hyp(sequence, id, prediction)
  result <- list(prediction = prediction,
                 sequence = seq)
  return(result)
}


