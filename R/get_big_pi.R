#' Scraping big-PI Plant Predictor web server.
#'
#' big-PI Plant Predictor is a web server utilizing a scoring algorithm for prediction of GPI modification sites in plants.
#'
#' @param data A data frame with protein amino acid sequences as strings in one column and corresponding id's in another. Alternatively a path to a .fasta file with protein sequences. Alternatively a list with elements of class "SeqFastaAA" resulting from seqinr::read.fasta call.
#' @param sequence A vector of strings representing protein amino acid sequences, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param id A vector of strings representing protein identifiers, or the appropriate column name if a data.frame is supplied to data argument. If .fasta file path, or list with elements of class "SeqFastaAA" provided to data, this should be left blank.
#' @param simplify A bolean indicating the type of returned object, defaults to TRUE
#' @param sleep A numeric indicating the pause in seconds between server calls, at default set to 1
#' @param verbose Bolean whether to print out the total score for each sequence, defaults to TRUE
#' @return If simplify == TRUE:
#' A data frame with columns:
#' \describe{
#'   \item{omega_site}{Character, indicating the sequence position of the highest scoring omega-site}
#'   \item{Quality}{Character, indicating the quality of the highest scoring omega-site}
#'   \item{PValue}{Numeric, indicating the p-value for the prediction of the highest scoring omega-site"}
#'   \item{id}{Character, indicating the protein identifier}
#'   \item{is.bigpi}{Logical, did big-Pi predict the presence of a GPI}
#'   }
#'
#' If simplify == FALSE:
#' A list of predictions, each element named according to the sequence id, containing a two element list:
#' \describe{
#'   \item{prediction}{data frame, resembling the one returned by simplify == TRUE, along with alternative site predictions (if present)}
#'   \item{calculation}{data frame, with profile dependent and profile independent scores}
#'   }
#'
#' @source \url{http://mendel.imp.ac.at/gpi/plant_server.html}
#' @references Eisenhaber B. Wildpaner M. Schultz CJ. Borner GHH. Dupree P. Eisenhaber F. (2003) Glycosylphosphatidylinositol lipid anchoring of plant proteins. Sensitive prediction from sequence- and genome-wide studies for Arabidopsis and rice. Plant Physiology 133(4): 1691-701
#'
#' @examples
#' library(ragp)
#' data(at_nsp)
#'
#' #indexes of some sequences in at_nsp
#' ind <- c(129, 145, 147, 160, 170)
#'
#' big_pi_pred <- get_big_pi(sequence = at_nsp$sequence[ind],
#'                           id = at_nsp$Transcript.id[ind],
#'                           simplify = FALSE)
#'                           
#' big_pi_pred <- get_big_pi(data = at_nsp[ind,],
#'                           sequence = sequence,
#'                           id = Transcript.id,
#'                           simplify = TRUE)

#'@export



get_big_pi <- function(data = NULL, sequence, id, simplify = TRUE, sleep = NULL, verbose = TRUE){
  if (missing(simplify)){
    simplify <- TRUE
  }
  if (missing(sleep)){
    sleep <- 1
  }
  if (missing(verbose)){
    verbose <- TRUE
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
    id <- if(deparse(substitute(id)) %in% colnames(data)){
      data[[deparse(substitute(id))]]
    } else if(id %in% colnames(data)){
      data[[id]]
    } else {
      stop("specified id not found in data")
    }
    id <- as.character(id)  
    sequence  <- if(deparse(substitute(sequence)) %in% colnames(data)){
      data[[deparse(substitute(sequence))]]
    } else if(sequence %in% colnames(data)){
      data[[sequence]]
    } else {
      stop("specified id not found in data")
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
  if (length(sequence) != length(id)) 
    stop("id and sequence vectors are not of same length")
  url <- "http://mendel.imp.ac.at/gpi/plant_server.html"
  session <- rvest::html_session(url)
  form <- rvest::html_form(session)[[2]]
  n <- length(sequence)
  sequence <- as.character(sequence)
  aa_regex <- "[^ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv]"
  crop_1 <- ".*Use of the prediction function for VIRIDIPLANTAE"
  site_1 <- "Potential GPI-modification site was found"
  site_2 <- "Potential alternative GPI-modification site was found"
  num_reg <- "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?"
  Terms <- c("Total Score",
             "Profile Score",
             "Term  0",
             "Term  1",
             "Term  2",
             "Term  3",
             "Term  4",
             "Term  5",
             "Term  6",
             "Term  7",
             "Term  8",
             "Term  9",
             "Term 10",
             "Term 11",
             "Term 12",
             "Term 13",
             "Term 14",
             "Term 15",
             "Term 16",
             "Term 17",
             "Term 18",
             "Term 19",
             "Term 20",
             "Term 21",
             "Profile independent Score")
  Term_pos <- c(1, 1, 4, 4, 4, 4, 5,
                4, 4, 2, 4, 4, 4, 4,
                rep(3, 9), 4, 1)
  extract_val <- function(x, y){
    res <- regmatches(y[grep(x, y)], gregexpr(num_reg, y[grep(x, y)]))
    res <- unlist(res)
    res <- as.numeric(res)
    return(res)
  }
  big_pi_list  <-  vector("list", n)
  for (i in 1:n){
    if (grepl(aa_regex, sequence[i])){
      warning(paste("sequence",
                    "[", id[i], "]",
                    " contains symbols not corresponding to amino acids",
                    sep = ""), call. = FALSE)
      
      Best <- rep(NA, length(Terms))
      Alternative <- rep(NA, length(Terms))
      total_score <- NA
      pred <- data.frame(omega_site = rep("invalid character", 2), 
                         Quality = rep("None", 2),
                         PValue =  rep(NA, 2),
                         stringsAsFactors = FALSE)
      
      calc <- data.frame(Terms,
                         Best,
                         Alternative,
                         stringsAsFactors = FALSE)
      res <- list(prediction = pred,
                  calculation = calc)
    } else {
      if (nchar(sequence[i]) <= 55){
        Best <- rep(NA, length(Terms))
        Alternative <- rep(NA, length(Terms))
        total_score <- NA
        pred <- data.frame(omega_site = rep("short sequence", 2), 
                           Quality = rep("None", 2),
                           PValue =  rep(NA, 2),
                           stringsAsFactors = FALSE)
        
        calc <- data.frame(Terms,
                           Best,
                           Alternative,
                           stringsAsFactors = FALSE)
        res <- list(prediction = pred,
                    calculation = calc)
   
      } else {
        form <- rvest::set_values(form, Sequence = sequence[i])
        form_res <- suppressMessages(rvest::submit_form(session, form))
        resulti <- rvest::html_text(rvest::html_nodes(form_res, "pre"))
        impro <- unlist(strsplit(resulti, "\\\n"))
        impro <- impro[grep(crop_1, impro):length(impro)]
        if (length(grep(site_1, impro)) > 0 &
            length(grep(site_2, impro)) > 0){
          Best <- unlist(lapply(seq_along(Term_pos), function(x){
            extract_val(x = Terms[x], y = impro)[Term_pos[x]]
          }))
          Alternative <- unlist(lapply(seq_along(Term_pos), function(x){
            extract_val(x = Terms[x], y = impro)[Term_pos[x]+1]
          }))
          
          position_gpi <- extract_val("Sequence position of the omega-site",
                                      impro)
          total_score <- extract_val("Total Score\\.",
                                     impro)[1]
          best_p <- extract_val("PValue = ", impro)[2]
          alt_p <- extract_val("PValue = ", impro)[4]
          
          site_q <- unlist(strsplit(impro[grep("Quality of the site", impro)],
                                    " {2,}", perl = TRUE))[c(2, 4)]
          
          pred <- data.frame(omega_site = position_gpi, 
                             Quality = site_q ,
                             PValue =  c(best_p, alt_p),
                             stringsAsFactors = FALSE)
          
          calc <- data.frame(Terms,
                             Best,
                             Alternative,
                             stringsAsFactors = FALSE)
          res <- list(prediction = pred,
                      calculation = calc)
        } else {
          if (length(grep(site_1, impro)) > 0 &
              length(grep(site_2, impro)) == 0){
            Best <- unlist(lapply(seq_along(Term_pos), function(x){
              extract_val(x = Terms[x], y = impro)[Term_pos[x]]
            }))
            Alternative <- rep(NA, length(Terms))
            position_gpi <- extract_val("Sequence position of the omega-site",
                                        impro)
            total_score <- extract_val("Total Score\\.",
                                       impro)[1]
            best_p <- extract_val("PValue = ", impro)[2]
            site_q <- unlist(strsplit(impro[grep("Quality of the site", impro)],
                                      " {2,}", perl = TRUE))[2]
            
            pred <- data.frame(omega_site = c(position_gpi[1], NA), 
                               Quality = c(site_q[1], NA) ,
                               PValue =  c(best_p, NA),
                               stringsAsFactors = FALSE)
            calc <- data.frame(Terms,
                               Best,
                               Alternative,
                               stringsAsFactors = FALSE)
            res <- list(prediction = pred,
                        calculation = calc)
          } else {
            Best <- unlist(lapply(seq_along(Term_pos), function(x){
              extract_val(x = Terms[x], y = impro)[Term_pos[x]]
            }))
            Alternative <- rep(NA, length(Terms))
            best_p <- extract_val("PValue = ", impro)[2]
            position_gpi <- extract_val("Among all positions checked",
                                        impro)
            total_score <- extract_val("Total Score\\.",
                                       impro)[1]
            
            pred <- data.frame(omega_site = c(position_gpi[1], NA), 
                               Quality = c("None", NA) ,
                               PValue =  c(best_p, NA),
                               stringsAsFactors = FALSE)
            calc <- data.frame(Terms,
                               Best,
                               Alternative,
                               stringsAsFactors = FALSE)
            res <- list(prediction = pred,
                        calculation = calc)

          }
        }
      }
    }
    big_pi_list[[i]] <- res
    if (verbose == TRUE){
      cat(paste(paste0(id[i], ":"), " Total Score:", total_score, "\n"))
      utils::flush.console()
      }
    Sys.sleep(sleep)
  }
  if(simplify){
    big_pi_list <- lapply(big_pi_list, function(x){
      return(x$pred[1,])
    })
    big_pi_list <- do.call(rbind, big_pi_list)
    big_pi_list$id <- id
    big_pi_list$is.bigpi <- big_pi_list$Quality != "None"
  } else {
    (names(big_pi_list) <- id)
  }
  
  return(big_pi_list)
  }
         
