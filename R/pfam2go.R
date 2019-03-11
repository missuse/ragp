#' Add GO terms based on pfam accessions
#'
#' The objective of gene ontology (GO) is to provide controlled vocabularies for the description of the biological process, molecular function, and cellular component of gene products. This function maps existing PFAM accessions to corresponding GO terms.
#'
#' @param data_pfam a data frame containing a column with PFAM accessions
#' @param pfam a string defining the column name where the PFAM accessions are stored. Defaults to "acc" as per output of get_hmm function.
#' @return  A merged data frame with columns as provided in data_pfam argument and additional columns:
#' \describe{
#'   \item{Pfam_acc}{Character, PFAM family accession.}
#'   \item{Pfam_name}{Character, PFAM family name.}
#'   \item{GO_name}{Character, GO term name.}
#'   \item{GO_acc}{Character, GO term accession.}
#'   }
#'
#'@source \url{http://geneontology.org/external2go/pfam2go}  
#'
#'        \url{ftp://ftp.geneontology.org/pub/go/external2go/pfam2go}
#'
#'@seealso \code{\link[ragp]{get_hmm}}
#'
#'@examples
#'\dontrun{
#'
#' library(ragp)
#' data(at_nsp)
#'
#' pfam_pred <- get_hmm(sequence = at_nsp$sequence[1],
#'                      id = at_nsp$Transcript.id[1])
#'
#' pfam_pred_go <- pfam2go(data_pfam = pfam_pred, pfam = "acc")
#' }
#' 
#'@export

pfam2go <- function(data_pfam, pfam){
  if (missing(pfam)){
    pfam = "acc"
  }
  if (length(data_pfam[[pfam]]) == 0){
    stop ("please provide the column name of the PFAM accesions, in argument pfam",
          call. = FALSE)
  }
  go_text <- readLines("ftp://ftp.geneontology.org/pub/go/external2go/pfam2go")
  go_text <- go_text[grep("^Pfam:", go_text)]
  go_text <- strsplit(go_text, "{0,} > {0,}")
  go_text <- do.call(rbind, go_text)
  go_1 <- strsplit(go_text[,1], " +")
  go_1 <- do.call(rbind, go_1)
  go_1[,1] <- gsub("Pfam:", "", go_1[,1])
  go_2 <- strsplit(go_text[,2], "{0,} ; {0,}")
  go_2 <- do.call(rbind, go_2)
  go <- data.frame(go_1, go_2, stringsAsFactors = FALSE)
  colnames(go) <- c("Pfam_acc", "Pfam_name", "GO_name", "GO_acc")
  data_pfam[["acc_temp"]] <- substring(data_pfam[[pfam]], first = 1, last = 7)
  data_pfam[["rownames_temp"]] <- 1:nrow(data_pfam)
  out <- merge.data.frame(data_pfam,
                          go,
                          by.x = "acc_temp",
                          by.y = "Pfam_acc",
                          all.x = TRUE, 
                          all.y = FALSE,
                          sort = FALSE)
  out <- out[order(out[["rownames_temp"]]),]
  out <- out[,setdiff(names(out),
                      c("acc_temp",
                        "rownames_temp"))]
  return(out)
  }
