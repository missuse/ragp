% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_hyp.R
\name{predict_hyp}
\alias{predict_hyp}
\alias{predict_hyp.default}
\alias{predict_hyp.character}
\alias{predict_hyp.data.frame}
\alias{predict_hyp.list}
\alias{predict_hyp.AAStringSet}
\title{Predict hydroxyproline positions in plant proteins based on primary structure}
\usage{
predict_hyp(data, ...)

\method{predict_hyp}{character}(data, ...)

\method{predict_hyp}{data.frame}(data, sequence, id, ...)

\method{predict_hyp}{list}(data, ...)

\method{predict_hyp}{default}(
  data = NULL,
  sequence,
  id,
  version = "V2",
  tprob = ifelse(version == "V1", 0.3, 0.224),
  split = 1,
  ...
)

\method{predict_hyp}{AAStringSet}(data, ...)
}
\arguments{
\item{data}{A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class \code{\link[seqinr]{SeqFastaAA}} resulting from \code{\link[seqinr]{read.fasta}} call. Alternatively an \code{\link[Biostrings]{AAStringSet}} object. Should be left blank if vectors are provided to sequence and id arguments.}

\item{...}{currently no additional arguments are accepted apart the ones documented bellow.}

\item{sequence}{A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.}

\item{id}{A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank. Ids should be unique.}

\item{version}{A string indicating which model version to use: the first version "V1", or the second version "V2". Default is "V2".}

\item{tprob}{A numeric value indicating the threshold for prediction. Acceptable values are in 0 - 1 range. At default set to 0.3 for "V1" model and 0.224 for "V2" model, offering a tradeoff between sensitivity and specificity.}

\item{split}{A numeric value determining the ratio of vectorized and sequential computation. Should be left at default, lower to 0 - 1 range if low memory errors occur. Increase at your own risk.}
}
\value{
A list with two elements:
\describe{
  \item{prediction}{data frame with columns:
  id - character, indicating the inputted protein id's;
  substr - character, indicating the sequence substring which was used for predictions;
  P_pos - integer, position of proline in the sequence;
  prob - numeric, predicted probability of being hydroxyproline;
  HYP - character, is the site predicted as a hydroxyproline}
  \item{sequence}{data frame with columns:
  sequence - sequences with prolines - P substituted with hydroxyprolines - O according to the prediction;
  id - corresponding id's}
}
}
\description{
predict_hyp is a hydroxyproline site prediction algorithm for plant proteins, based on the xgboost distributed gradient boosting library.
It was trained on plant sequences with experimentally determined 4-hydroxyprolines from swissprot data base. Prediction is not possible for prolines which are within 10 N-terminal and 6 C-terminal amino acids (V1 model version) and 10 N-terminal and 7 C-terminal amino acids (V2 model version), they will be excluded from output.
}
\details{
Previously trained xgboost models were re-saved using xgboost 1.1.1.1 to increase compatibility. While using the mentioned xgboost version the returned predictions are equal to previous. However, using earlier xgboost versions with the new models will result in slightly different predicted probabilities.
}
\examples{
library(ragp)
data(at_nsp)

#ramdom indexes
ind <- c(129, 145, 147, 160, 170,
    180, 189, 203, 205, 214, 217, 224)

hyp_pred <- predict_hyp(sequence = at_nsp$sequence[ind],
                        id = at_nsp$Transcript.id[ind])

}
