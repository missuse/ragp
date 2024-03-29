% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_big_pi.R
\name{get_big_pi}
\alias{get_big_pi}
\alias{get_big_pi.default}
\alias{get_big_pi.character}
\alias{get_big_pi.data.frame}
\alias{get_big_pi.list}
\alias{get_big_pi.AAStringSet}
\title{Query big-PI Plant Predictor web server.}
\source{
\url{http://mendel.imp.ac.at/gpi/plant_server.html}
}
\usage{
get_big_pi(data, ...)

\method{get_big_pi}{character}(data, ...)

\method{get_big_pi}{data.frame}(data, sequence, id, ...)

\method{get_big_pi}{list}(data, ...)

\method{get_big_pi}{default}(
  data = NULL,
  sequence,
  id,
  simplify = TRUE,
  sleep = 1,
  progress = FALSE,
  ...
)

\method{get_big_pi}{AAStringSet}(data, ...)
}
\arguments{
\item{data}{A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class \code{\link[seqinr]{SeqFastaAA}} resulting from \code{\link[seqinr]{read.fasta}} call. Alternatively an \code{\link[Biostrings]{AAStringSet}} object. Should be left blank if vectors are provided to sequence and id arguments.}

\item{...}{currently no additional arguments are accepted apart the ones documented bellow.}

\item{sequence}{An appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.}

\item{id}{An appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.}

\item{simplify}{A boolean indicating the type of returned object, defaults to TRUE.}

\item{sleep}{A numeric indicating the pause in seconds between server calls, at default set to 1.}

\item{progress}{Boolean, whether to show the progress bar, at default set to FALSE.}
}
\value{
If simplify == TRUE:

A data frame with columns:
\describe{
  \item{id}{Character, indicating the protein identifier}
  \item{is.gpi}{Logical, did big-Pi predict the presence of a GPI}
  \item{Quality}{Character, indicating the quality of the highest scoring omega-site}
  \item{omega_site}{Integer, indicating the sequence position of the highest scoring omega-site}
  \item{PValue}{Numeric, indicating the p-value for the prediction of the highest scoring omega-site}
  }

If simplify == FALSE:

A list of predictions, each element named according to the sequence id, containing a two element list:
\describe{
  \item{prediction}{data frame, resembling the one returned by simplify == TRUE, along with alternative site predictions (if present)}
  \item{calculation}{data frame, with profile dependent and profile independent scores}
  }
}
\description{
big-PI Plant Predictor is a web server utilizing a scoring algorithm for prediction of GPI modification sites in plants.
}
\note{
If the server is unable to make a prediction due to non-amino acid letters or length of the sequence, the returned prediction is NA (is.bigpi column).
}
\examples{
library(ragp)
data(at_nsp)

#indexes of some sequences in at_nsp
ind <- c(129, 145, 147, 160, 170)

big_pi_pred <- get_big_pi(sequence = at_nsp$sequence[ind],
                          id = at_nsp$Transcript.id[ind],
                          simplify = FALSE)

big_pi_pred <- get_big_pi(data = at_nsp[ind,],
                          sequence = sequence,
                          id = Transcript.id,
                          simplify = TRUE)
big_pi_pred

}
\references{
Eisenhaber B. Wildpaner M. Schultz CJ. Borner GHH. Dupree P. Eisenhaber F. (2003) Glycosylphosphatidylinositol lipid anchoring of plant proteins. Sensitive prediction from sequence- and genome-wide studies for Arabidopsis and rice. Plant Physiology 133(4): 1691-701
}
