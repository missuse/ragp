#' Motif and amino acid bias of a protein sequence
#'
#' Calculate the compositional bias of whole protein sequences, or slide along the sequences and calculate the bias for all subsequences of a defined length (frame). Sequence motifs can also be searched.
#' 
#' @param sequence A vector of strings representing protein amino acid sequences
#' @param id A vector of strings representing the names of the corresponding sequences
#' @param frame An optional integer defining the frame length for sliding along the sequences
#' @param custom An optional character vector which can contain amino acid one letter symbols eg. c("P", "V", "K", "C", "T"), strings corresponding to amino acid motifs eg c("PTYK", "PVKC"), mixing of types is supported. Common regex operators can also be utilized.
#' 
#' @return by default if, frame is not specified a data frame is returned with columns:
#' \enumerate{
#' \item{past_count:} {summed count of "P", "A", "S" and "T" amino acids}
#' \item{pvyk_count:} {summed count of "P", "V", "Y" and "K" amino acids}
#' \item{psky_count:} {summed count of "P", "S", "K" and "Y" amino acids}
#' \item{p_count:} {count of "P"}
#' \item{ext_sp_count:} {summed count of SPPP+ motifs}
#' \item{ext_fyxy_count:} {summed count of [FY]XY motifs}
#' \item{ext_khy_count:} {summed count of KHY motifs}
#' \item{ext_vyhkde_count:} {summed count of VY[HKDE] motifs}
#' \item{ext_vxy_count:} {summed count of VXY motifs}
#' \item{ext_yy_count:} {summed count of YY motifs}
#' \item{prp_ppvqk_count:} {summed count of PPV[QK] motifs}
#' \item{prp_ppvxkt_count:} {summed count of PPVX[KT] motifs}
#' \item{ppr_kkpcpp:} {summed count of KKPCPP motifs}
#' \item{past_percent:} {summed percent of "P", "A", "S" and "T" amino acids}
#' \item{pvyk_percent:} {summed percent of "P", "V", "Y" and "K" amino acids}
#' \item{psky_percent:} {summed percent of "P", "S", "K" and "Y" amino acids}
#' \item{p_percent:} {percent of "P"}
#' \item{id:} {names of the corresponding sequences, as per input}
#'}
#' 
#' if frame is specified a list with two elements is returned:
#' \enumerate{
#' \item{protein:} {a named list with names corresponding to sequence id's, where each element is a data frame resembling the one returned when frame = NULL, each row of which corresponds to each sequence frame}
#' \item{maab:} {a character vector with search criteria}
#' }
#' 
#' if custom is specified the returned object contains counts for each element specified in custom; additionally percentage is calculated for all elements containing a single symbol.  
#' 
#' @details By default the function provides motif and amino acid bias descriptors used for classification of HRGP's by the MAAB pipeline (Johnson et al. 2017). Arabinogalactan descriptors used in the MAAB pipeline other than the percent of PAST amino acids are not incorporated in this function. Users can specify a sliding frame window in order to calculate partial sequence descriptors or specify custom search criteria.
#' 
#' @references Johnson KL, Cassin AM, Lonsdale A, Bacic A, Doblin MS, Schultz CJ. (2017) Pipeline to Identify Hydroxyproline-Rich Glycoproteins. Plant Physiol 174(2): 886-903.
#'
#' @seealso \code{\link[ragp]{scan_ag}}
#'
#' @examples
#'
#' library(ragp)
#' data(at_nsp)
#'
#' test_bias <- maab(sequence = at_nsp$sequence,
#'                   id = at_nsp$Transcript.id)
#'                   
#' test_bias <- maab(sequence = at_nsp$sequence[1:20],
#'                   id = at_nsp$Transcript.id[1:20],
#'                   frame = 30)
#'                   
#' test_bias <- maab(sequence = at_nsp$sequence[1:20],
#'                   id = at_nsp$Transcript.id[1:20],
#'                   custom = c("KKP", "SD[GD], "P"))                       
#'                   
#' @export                  


maab <- function(sequence, id, frame = NULL, custom = NULL){
  sequence <- toupper(as.character(sequence))
  id <- as.character(id)
  seq_len <- nchar(sequence)
  if (!missing(frame)){
    if (!all.equal(frame, as.integer(frame), check.attributes = FALSE)){
      stop ("frame must be an integer")
    }
  }
  amino_acids = c("A", "R", "N", "D", "C",
                  "Q", "E", "G", "H", "I",
                  "L", "K", "M", "F", "P",
                  "S", "T", "W", "Y", "V")
  if (!missing(custom)){
    symbols <- toupper(gsub("[^\\pL']+", "", custom, perl = T))
    if(any(!unlist(strsplit(symbols, split = "")) %in% amino_acids)){
      warning("custom contains characters other than amino acid symbols")
    }
    }
    if (missing(custom)){
    hrgp_aa <- c( "P|A|S|T",
                  "P|V|Y|K",
                  "P|S|K|Y",
                  "P",
                  "SP{3,}",
                  "[FY].Y",
                  "KHY",
                  "VY[HKDE]",
                  "V.Y",
                  "YY",
                  "PPV[QK]",
                  "PPV.[KT]",
                  "KKPCPP")
    hrgp_name <- c("past_count",
                   "pvyk_count",
                   "psky_count",
                   "p_count",
                   "ext_sp_count",
                   "ext_fyxy_count",
                   "ext_khy_count",
                   "ext_vyhkde_count",
                   "ext_vxy_count",
                   "ext_yy_count",
                   "prp_ppvqk_count",
                   "prp_ppvxkt_count",
                   "ppr_kkpcpp")
    if (missing(frame)){
      out_count <-  lapply(hrgp_aa, function(x){
        stringr::str_count(sequence, x)
      })
      names(out_count) <- hrgp_name
      out_percent <- lapply(1:4, function(x){
        out_count[[x]] / seq_len * 100
      })
      names(out_percent) <- paste0(c("past",
                                     "pvyk",
                                     "psky",
                                     "p"), "_percent")
      out <- c(out_count,
               out_percent)
      out <- do.call(cbind, out)
      out <- as.data.frame(out,
                           stringsAsFactors = FALSE)
      out$id <- id
      return(out)
    } else {
      res <- lapply(sequence, function(x){
        if (frame > nchar(x)) frame <- nchar(x)
        aa_out <- lapply(hrgp_aa, function(k){
          P_count <- vector("numeric", length = (nchar(x) - frame + 1))
          for (i in 1 : (nchar(x) - frame + 1)){
            P_count[i] <- stringr::str_count(substr(x,
                                                    start = i,
                                                    stop = i + frame - 1), k)
            
          }
          P_count 
        }
        )
        names(aa_out) <- hrgp_name 
        out_percent <- lapply(1:4, function(x){
          aa_out[[x]] / frame * 100
        }
        )
        
        names(out_percent) <- paste0(c("past",
                                       "pvyk",
                                       "psky",
                                       "p"), "_percent")
        out <- c(aa_out,
                 out_percent)
        out <- do.call(cbind, out)
        out <- as.data.frame(out,
                             stringsAsFactors = FALSE)
        out$frame_start <- 1 : (nchar(x) - frame + 1)
        out$frame_end <- frame  : nchar(x)
        return(out)
      }
      )
      names(res) <- id
      out <- list(protein = res,
                  maab = hrgp_aa)
      return(out)
    }
  }
  if (!missing(custom)){
    hrgp_aa <- custom
    hrgp_name <- tolower(gsub("[^\\pL']+", "", hrgp_aa, perl = T))
    if (missing(frame)){
      out_count <-  lapply(hrgp_aa, function(x){
        stringr::str_count(sequence, x)
      }
      )
      names(out_count) <- paste0(hrgp_name,
                                 "_count")
      perc <- which(nchar(hrgp_aa) == 1)
      if (length(perc) > 0){
        out_percent <- lapply(perc, function(x){
          out_count[[x]] / seq_len * 100
          }
          )
        names(out_percent) <- paste(hrgp_name[perc],
                                    "percent",
                                    sep = "_")
        out <- c(out_count, out_percent)
        } else {
          out <- out_count
          }
      out <- do.call(cbind, out)
      out <- as.data.frame(out, stringsAsFactors = FALSE)
      out$id <- id
      return(out)
    } else {
      res <- lapply(sequence, function(x){
        if (frame > nchar(x)) frame <- nchar(x)
        aa_out <- lapply(hrgp_aa, function(k){
          P_count <- vector("numeric", length = (nchar(x) - frame + 1))
          for (i in 1 : (nchar(x) - frame + 1)){
            P_count[i] <- stringr::str_count(substr(x,
                                                    start = i,
                                                    stop = i + frame - 1), k)
            
          }
          P_count 
        }
        )
        names(aa_out) <- paste0(hrgp_name, "_count")
        perc <- which(nchar(hrgp_aa) == 1)
        if (length(perc) > 0){
          out_percent <- lapply(perc, function(x){
            aa_out[[x]] / frame * 100
          }
          )
          names(out_percent) <- paste(hrgp_name[perc],
                                      "percent",
                                      sep = "_")
          out <- c(aa_out,
                   out_percent)
        } else {
          out <- aa_out
        }
        out <- do.call(cbind, out)
        out <- as.data.frame(out,
                             stringsAsFactors = FALSE)
        out$frame_start <- 1 : (nchar(x) - frame + 1)
        out$frame_end <- frame  : nchar(x)
        out
        }
        )
      names(res) <- id
      out <- list(protein = res,
                  maab = hrgp_aa)
      return(out)
    }
  }
}

  

  
  
  
    