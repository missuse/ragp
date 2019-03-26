#' Normalized Moreau-Broto Autocorrelation Descriptor
#'
#' This function calculates the normalized Moreau-Broto autocorrelation descriptor with lag 6.
#'  
#' @param x A vector of k-mers
#' @param mat A matrix with amino acid atributes. Columns are amino acids, rows are atributes
#' @param nlag The lag parameter
#'  
#' @return A matrix with the Normalized Moreau-Broto Autocorrelation Descriptor (dim: nlag * nrow(mat))
#' 
#' @details Normalized Moreau-Broto Autocorrelation Descriptor is calculated based on the following amino acid attributes:
#' CIDH920105 - Normalized average hydrophobicity scales, BHAR880101 - Average flexibility indices,
#' CHAM820102 - Free energy of solution in water, kcal/mole, BIGC670101 - Residue volume,
#' CHAM810101 - Polarizability parameter and DAYM780201 - Relative mutability 
#' 
#' @note This is an internal function used inside predict_hyp
#' 
#' @author original R code by Nan Xiao, modified by Milan Dragićević
#' 
#' @seealso \code{\link[ragp]{predict_hyp} \link[protr]{extractMoreauBroto}}
#' 
#' @import stats


extractMBdesc <- function(x, mat, nlag){
  N <- nchar(x[1])
  mat_names <- rownames(mat)
  mat_rows <- nrow(mat)
  pmean <- rowMeans(mat)
  psd <- apply(mat, 1, stats::sd) * sqrt((20 - 1)/20)
  mat <- (mat - pmean)/psd
  m1 <- t(mat)[strsplit(paste(x,
                              collapse = ""),
                        "")[[1]], ]
  P <- do.call(cbind,
               lapply(1:N,
                      function(i) m1[seq(i,
                                         nrow(m1),
                                         by = N), ]))
  P <- do.call(cbind, lapply(1:mat_rows,
                             function(x) P[,seq(x,
                                                ncol(P),
                                                by = mat_rows)]))
  P <- matrix(t(P), byrow = T, ncol = N)
  MB <- vector("list", nlag)
  for (j in 1:nlag) {
    MB[[j]] <- rowSums(P[, 1:(N - j)] * P[, 1:(N - j) + j])/(N - j)
  }
  MB <- do.call(cbind, MB)
  MB <- matrix(t(MB),
               byrow = T,
               ncol = nlag * mat_rows)
  colnames(MB) <- as.vector(t(outer(mat_names,
                                    paste(".lag",
                                          1:nlag,
                                          sep = ""),
                                    paste, sep = "")))
  return(MB)
}


#' Quasi-Sequence-Order (QSO) Descriptor
#'
#' This function calculates the Quasi-Sequence-Order (QSO) descriptor with lag 12
#' 
#' @param m A vector of k-mers
#' 
#' @return  A matrix with the Quasi-Sequence-Order Descriptor (dim: 64)
#' 
#' @details QSO descriptor is calculated based on two physicochemical distance matrices: Schneider-Wrede (Schneider and Wrede, 1994) and Grantham physicochemical distance matrix (Grantham, 1974)
#' 
#' @note This is an internal function used inside predict_hyp
#' 
#' @author original R code by Nan Xiao, modified by Milan Dragićević
#' 
#' @references Kuo-Chen Chou. Prediction of Protein Subcellar Locations by Incorporating Quasi-Sequence-Order Effect. Biochemical and Biophysical Research Communications, 2000, 278, 477-483.
#' Gisbert Schneider and Paul Wrede. The Rational Design of Amino Acid Sequences by Artifical Neural Networks and Simulated Molecular Evolution: Do Novo Design of an Idealized Leader Cleavage Site. Biophys Journal, 1994, 66, 335-344.
#' 
#' @seealso \code{\link[ragp]{predict_hyp} \link[protr]{extractQSO}}


QSOlevel <- function (m){
  QSObabe <- function(x) {
    if (length(x) == 1) {
      x <- c(x, x)
    }
    w <- 0.1
    nlag <- 12
    N <- nchar(x[1])
    xSplitted <- strsplit(as.character(x),
                          split = "")
    xSplitted <- do.call(rbind,
                         xSplitted)
    tau1 <- list()
    for (d in 1:nlag) {
      for (i in 1:(N - d)) {
        tau1[[length(tau1) + 1]] <- diag(as.matrix(DistMat1)[xSplitted[, i],
                                                             xSplitted[, i + d]])^2
      }
    }
    tau1 <- t(do.call(rbind, tau1))
    ends <- cumsum(N - 1:nlag)
    starts <- c(1, ends[-length(ends)] + 1)
    tau1 <- do.call(cbind,
                    lapply(mapply(seq,
                                  starts,
                                  ends),
                           function(s) if(length(s) > 1){
                             rowSums(tau1[,s])
                           } else {
                             tau1[,s]
                           }
                    )
    )
    
    fr <- sapply(AADict, function(y) sapply(gregexpr(y,
                                                     as.character(x)),
                                            function(x) sum(x > 0)))
    Xr1 <- fr/(1 + (w * rowSums(tau1)))
    colnames(Xr1) <- paste("Schneider.Xr.",
                           colnames(Xr1), 
                           sep = "")
    Xd1 <- (w * tau1)/(1 + (w * rowSums(tau1)))
    colnames(Xd1) <- paste("Schneider.Xd.",
                           1:nlag,
                           sep = "")
    tau2 <- list()
    for (d in 1:nlag) {
      for (i in 1:(N - d)) {
        tau2[[length(tau2) + 1]] <- diag(as.matrix(DistMat2)[xSplitted[, i],
                                                             xSplitted[, i + d]])^2
      }
    }
    tau2 <- t(do.call(rbind, tau2))
    tau2 <- do.call(cbind,
                    lapply(mapply(seq,
                                  starts,
                                  ends),
                           function(s) if(length(s) > 1){
                             rowSums(tau2[,s])
                           } else {
                             tau2[,s]
                           }
                    )
    )
    Xr2 <- fr/(1 + (w * rowSums(tau2)))
    colnames(Xr2) <- paste("Grantham.Xr.",
                           colnames(Xr2), 
                           sep = "")
    Xd2 <- (w * tau2)/(1 + (w * rowSums(tau2)))
    colnames(Xd2) <- paste("Grantham.Xd.",
                           1:nlag,
                           sep = "")
    QSO <- cbind(Xr1, Xr2, Xd1, Xd2)
    if (length(unique(x)) == 1) {
      QSO <- QSO[1, , drop = FALSE]
    }
    return(QSO)
  }
  splt <- 180
  k <- length(m)/splt
  if (k < 1) {
    QSO_all <- QSObabe(m)
  }
  else {
    pam <- ((seq(length(m)) - 1)%/%splt) + 1
    m_split <- split(m, pam)
    QSO_all <- lapply(m_split, function(x) QSObabe(x))
    QSO_all <- do.call(rbind, QSO_all)
  }
  return(QSO_all)
}



#' Composition descriptor (CTDC) 
#'
#' This function calculates the Composition descriptor
#' 
#' @param x A vector of k-mers (character vector)
#' 
#' @return  A matrix with the Composition descriptor (dim: 21)
#' 
#' @note This is an internal function used inside predict_hyp
#' 
#' @author original R code by Nan Xiao
#' 
#' @references Inna Dubchak, Ilya Muchink, Stephen R. Holbrook and Sung-Hou Kim. Prediction of protein folding class using global description of amino acid sequence. Proceedings of the National Academy of Sciences. USA, 1995, 92, 8700-8704.
#' 
#' @seealso \code{\link[ragp]{predict_hyp} \link[protr]{extractCTDC}}

CTDC <- function (x){
  group1 <- list(hydrophobicity = c("R", "K", "E", "D", "Q", "N"),
                 normwaalsvolume = c("G", "A", "S", "T", "P", "D", "C"),
                 polarity = c("L", "I", "F", "W", "C", "M", "V", "Y"),
                 polarizability = c("G", "A", "S", "D", "T"),
                 charge = c("K", "R"),
                 secondarystruct = c("E", "A", "L", "M", "Q", "K", "R", "H"),
                 solventaccess = c("A", "L", "F", "C", "G", "I", "V", "W"))
  group2 <- list(hydrophobicity = c("G", "A", "S", "T", "P", "H", "Y"),
                 normwaalsvolume = c("N", "V", "E", "Q", "I", "L"),
                 polarity = c("P", "A", "T", "G", "S"),
                 polarizability = c("C", "P", "N", "V", "E", "Q", "I", "L"),
                 charge = c("A", "N", "C", "Q", "G", "H", "I", "L", "M", "F", "P", "S", "T", "W", "Y", "V"),
                 secondarystruct = c("V", "I", "Y", "C", "W", "F", "T"),
                 solventaccess = c("R", "K", "Q", "E", "N", "D"))
  group3 <- list(hydrophobicity = c("C", "L", "V", "I", "M", "F", "W"),
                 normwaalsvolume = c("M", "H", "K", "F", "R", "Y", "W"),
                 polarity = c("H", "Q", "R", "K", "N", "E", "D"),
                 polarizability = c("K", "M", "H", "F", "R", "Y", "W"),
                 charge = c("D", "E"),
                 secondarystruct = c("G", "N", "P", "S", "D"),
                 solventaccess = c("M", "S", "P", "T", "H", "Y"))
  xSplitted <- strsplit(x, split = "")[[1]]
  n <- nchar(x)
  g1 <- lapply(group1,
               function(g) length(which(xSplitted %in% g)))
  names(g1) <- paste(names(g1),
                     "Group1",
                     sep = ".")
  g2 <- lapply(group2,
               function(g) length(which(xSplitted %in% g)))
  names(g2) <- paste(names(g2),
                     "Group2",
                     sep = ".")
  g3 <- lapply(group3,
               function(g) length(which(xSplitted %in% g)))
  names(g3) <- paste(names(g3),
                     "Group3",
                     sep = ".")
  CTDC <- unlist(c(g1, g2, g3))/n
  ids <- unlist(lapply(1:7, function(x) x + c(0, 7, 14)))
  CTDC[ids]
}

#' Extract amino acid attributes
#'
#' This function extracts amino acid attributes
#' 
#' @param x A vector of k-mers (character vector)
#' @param AA a matrix of amino acid attributes with attributes as rownames and amino acids as column names
#' 
#' @return A matrix with specified amino acid attributes (dim: nchar(x) * ncol(AA)
#'  
#' @note This is an internal function used inside predict_hyp
#' 
#' @seealso \code{\link[ragp]{predict_hyp}}

getAAindex <- function(x, AA){
  x <- as.character(x)
  m1 <- t(AA)[strsplit(paste(x, collapse = ""), "")[[1]], ]
  p <- nchar(x[1])
  if(length(x) == 1){
    res <- rbind(unlist(lapply(1:p, function(i) m1[seq(i, nrow(m1), by = p), ])))
  } else {
    res <- do.call(cbind, lapply(1:p, function(i) m1[seq(i, nrow(m1), by = p), ]))
  }
  z <- floor(p/2)
  var_names <- expand.grid(rownames(AA), -z:z)
  var_names <- paste(var_names$Var1, var_names$Var2, sep="_")
  colnames(res) <- var_names
  res <- res[, -grep("_0$", var_names), drop = FALSE] 
  row.names(res) <-   NULL
  class(res) <- "numeric"
  res
}

#' Extract simetric k-mers around prolines
#'
#' This function extracts simetric k-mers around prolines
#' 
#' @param sequence protein sequence in the form of a string
#' @param id protein id in the form of a string
#' @param kmer number of amino acids to retain on each side
#' 
#' @return A matrix with extracted kmers
#'  
#' @note This is an internal function used inside predict_hyp
#' 
#' @seealso \code{\link[ragp]{predict_hyp}}

getKmer <- function(sequence, id, kmer) {
  sequence <- as.character(sequence)
  sequence <- toupper(sequence)
  n_char <- nchar(sequence)
  P_pos <- unlist(gregexpr("P", sequence))
  P_mer <- sapply(P_pos, function(x) substr(sequence, start = x - 
                                              kmer, stop = x + kmer))
  P_mer <- cbind(id = as.character(rep(id, length(P_mer))), 
                 substr = as.character(P_mer), pos = as.character(P_pos), 
                 nchar = as.character(rep(n_char, length(P_pos))))
  return(P_mer)
}

#' Substitute prolines to hydroxyprolines 
#'
#' This function substitutes P with O based on predict_hyp prediction
#' 
#' @param sequence protein sequences in the form of a string vector
#' @param id protein ids in the form of a string vector
#' @param hyp the prediction object returned by predict_hyp
#' 
#' @return A vector of strings where P are changed to O based on predict_hyp prediction
#'  
#' @note This is an internal function used inside predict_hyp
#' 
#' @seealso \code{\link[ragp]{predict_hyp}}

sub_hyp <- function (sequence, id, hyp){
  sequence <- unlist(lapply(id, function(x) {
    hyp <- hyp[hyp$id == x, ]
    hyp <- hyp[!is.na(hyp$HYP),]
    sequencei <- as.character(sequence[id == x])
    if(nrow(hyp) == 0){
      return(sequencei)
    }
    p_pos <- as.numeric(
      as.character(
        hyp$P_pos[hyp$HYP == "Yes"]))
    for (i in p_pos) {
      substr(sequencei, start = i, stop = i) <- "O"
    }
    return(sequencei)
  }))
}




