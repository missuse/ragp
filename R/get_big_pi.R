#' Scraping big-PI Plant Predictor web server.
#'
#' big-PI Plant Predictor is a web server utlizing a scoring algorithm for prediction of GPI modification sites in plants.
#'
#' @param sequence A vector of strings representing protein amino acid sequences
#' @param short.out A bolean indicating the type of returned object, defaults to T
#' @param sleep A numeric indicating the pause in seconds betwean server calls, at default set to 2
#' @param verbose Bolean wheather to print out the total score for each sequence, defaults to T
#' @return If short.out==T:
#' A data frame with columns:
#' \describe{
#'   \item{quality_gpi}{Character, indicating the "Quality of the highest scoring omega-site"}
#'   \item{position_gpi}{Iinteger, indicating the "Sequence position of the highest scoring omega-site"}
#'   \item{total_score}{Numeric, indicating "Total Score of the highest scoring omega-site"}
#'   \item{profile_score}{Numeric, indicating the "Profile Score highest scoring omega-site"}
#'   \item{profile_i_score}{Numeric, indicating "Profile independent Score of the highest scoring omega-site"}
#'   \item{alt_quality_gpi}{Character, indicating the "Quality of the alternative omega-site", if present, else NA}
#'   \item{alt_position_gpi}{Integer, indicating the "Sequence position of the alternative omega-site", if present, else NA}
#'   \item{alt_total_score}{Numeric, indicating "Total Score of the alternative omega-site", if present, else NA}
#'   \item{alt_profile_score}{Numeric, indicating the "Profile Score of the alternative omega-site", if present, else NA}
#'   \item{alt_profile_i_score}{Numeric, indicating "Profile independent Score of the alternative omega-site" if present, else NA}
#'   }
#'
#' If short.out==F:
#' A list where each element correspondes to one queried sequence containing lines of text of the whole big-PI Plant Predictor output.
#'
#' @source \url{http://mendel.imp.ac.at/gpi/plant_server.html}
#' @references Eisenhaber B. Wildpaner M. Schultz CJ. Borner GHH. Dupree P. Eisenhaber F. (2003) Glycosylphosphatidylinositol lipid anchoring of plant proteins. Sensitive prediction from sequence- and genome-wide studies for Arabidopsis and rice. Plant Physiology 133(4): 1691-701
#'
#' @examples
#' library(ragp)
#' data(at_nsp)
#'
#' #indexes of some sequences known to have a positive prediction
#' ind <- c(129, 145, 147, 160, 170,
#'     180, 189, 203, 205, 214, 217, 224)
#'
#' big_pi_pred <- get_big_pi(sequence = at_nsp$sequence[ind])
#'
#' #extract terms not included in the short output:
#'
#' #first define a function for extraction
#'
#' extract_val <- function(x, y){
#'  res <- regmatches(y[grep(x, y)],
#'              gregexpr(num_reg, y[grep(x, y)]))
#'  res <- unlist(res)
#'  res <- as.numeric(res)
#'  return(res)
#' }
#'
#' #define a regex that extracts numerics
#' num_reg <- "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?"
#'
#' #run function with short.out=F
#' big_pi_pred <- get_big_pi(sequence = at_nsp$sequence[ind],
#'                           short.out = F)
#'
#'lapply(big_pi_pred,
#'       function(z) extract_val(x = "Term  3", y = z)[4])
#'
#'lapply(test_big_po,
#'       function(z) extract_val(x = "Total Score", y = z)[1])
#'
#'#interested in the alternative site:
#'lapply(big_pi_pred,
#'       function(x) extract_val(x = "Term  3", y = x)[5])
#'#NAs do not have the alternative site
#'
#'@export


get_big_pi <- function(sequence, short.out = T, sleep = NULL, verbose = NULL){
  short.out <- short.out
  if (missing(short.out)){
    short.out <- T
  }
  if (missing(sleep)){
    sleep <- 2
  }
  if (missing(verbose)){
    verbose <- T
  }
  url <- "http://mendel.imp.ac.at/gpi/plant_server.html"
  session <- rvest::html_session(url)
  form <- rvest::html_form(session)[[2]]
  n <- length(sequence)
  sequence <- as.character(sequence)
  big_pi <- data.frame(quality_gpi = vector("numeric", n),
                       position_gpi = vector("integer", n),
                       total_score = vector("numeric", n),
                       profile_score = vector("numeric", n),
                       profile_i_score = vector("numeric", n),
                       alt_quality_gpi = vector("numeric", n),
                       alt_position_gpi = vector("numeric", n),
                       alt_total_score = vector("numeric", n),
                       alt_profile_score = vector("numeric", n),
                       alt_profile_i_score = vector("numeric", n))
  big_pi_list  <-  vector("list", n)
  aa_regex <- "[^ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv]"
  crop_1 <- ".*Use of the prediction function for VIRIDIPLANTAE"
  site_1 <- "Potential GPI-modification site was found"
  site_2 <- "Potential alternative GPI-modification site was found"
  num_reg <- "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?"
  extract_val <- function(x, y){
    res <- regmatches(y[grep(x, y)], gregexpr(num_reg, y[grep(x, y)]))
    res <- unlist(res)
    res <- as.numeric(res)
    return(res)
  }

  if (short.out){
    for (i in 1:n){
      if (grepl(aa_regex, sequence[i])){
        warning(paste("sequence",
                      "[", i, "]",
                      " contains symbols not corresponding to amino acids",
                      sep = ""), call. = F)
        quality_gpi <- "None"
        position_gpi <- "invalid character"
        total_score <- NA
        profile_score <- NA
        profile_i_score <- NA
        alt_quality_gpi <- NA
        alt_position_gpi <- NA
        alt_total_score <- NA
        alt_profile_score <- NA
        alt_profile_i_score <- NA
      } else {
        if (nchar(sequence[i]) <= 55){
          quality_gpi <- "None"
          position_gpi <- "short sequence"
          total_score <- NA
          profile_score <- NA
          profile_i_score <- NA
          alt_quality_gpi <- NA
          alt_position_gpi <- NA
          alt_total_score <- NA
          alt_profile_score <- NA
          alt_profile_i_score <- NA
      } else {
        form <- rvest::set_values(form, Sequence = sequence[i])
        form_res <- suppressMessages(rvest::submit_form(session, form))
        resulti <- rvest::html_text(rvest::html_nodes(form_res, "pre"))
        impro <- unlist(strsplit(resulti, "\\\n"))
        impro <- impro[grep(crop_1, impro):length(impro)]
        if (length(grep(site_1, impro)) > 0 &
           length(grep(site_2, impro)) == 0){
          quality_gpi <- unlist(
            strsplit(impro[grep("Quality of the site", impro)],
                     " {2,10}", perl = T))[2]
          position_gpi <- extract_val("Sequence position of the omega-site",
                                   impro)
          total_score <- extract_val("Total Score\\.",
                                  impro)
          profile_score <- extract_val("Profile Score\\.",
                                    impro)
          profile_i_score <- extract_val("Profile independent Score",
                                      impro)
          alt_quality_gpi <- NA
          alt_position_gpi <- NA
          alt_total_score <- NA
          alt_profile_score <- NA
          alt_profile_i_score <- NA
          } else {
            if (length(grep(site_1, impro)) > 0 &
                length(grep(site_2, impro)) > 0){
          quality_gpi <- unlist(
            strsplit(head(impro[grep("Quality of the site", impro)], 1),
                     " {2,10}", perl = T))[2]
          position_gpi <- extract_val("Sequence position of the omega-site",
                                   impro)[1]
          total_score <- extract_val("Total Score\\.",
                                  impro)[1]
          profile_score <- extract_val("Profile Score\\.",
                                    impro)[1]
          profile_i_score <- extract_val("Profile independent Score",
                                      impro)[1]
          alt_quality_gpi <- unlist(
            strsplit(tail(impro[grep("Quality of the site", impro)], 1),
                     " {2,10}", perl = T))[2]
          alt_position_gpi <- extract_val("Sequence position of the omega-site",
                                       impro)[2]
          alt_total_score <- extract_val("Total Score\\.",
                                      impro)[2]
          alt_profile_score <- extract_val("Profile Score\\.",
                                        impro)[2]
          alt_profile_i_score <- extract_val("Profile independent Score",
                                          impro)[2]
        } else {
          quality_gpi <- "None"
          position_gpi <- extract_val("Among all positions checked",
                                   impro)
          total_score <- extract_val("Total Score\\.",
                                  impro)[1]
          profile_score <- extract_val("Profile Score\\.",
                                    impro)
          profile_i_score <- extract_val("Profile independent Score",
                                      impro)
          alt_quality_gpi <- NA
          alt_position_gpi <- NA
          alt_total_score <- NA
          alt_profile_score <- NA
          alt_profile_i_score <- NA
        }
        }
      }
      }
      big_pi$quality_gpi[i] <- quality_gpi
      big_pi$position_gpi[i] <- position_gpi
      big_pi$total_score[i] <- total_score
      big_pi$profile_score[i] <- profile_score
      big_pi$profile_i_score[i] <- profile_i_score
      big_pi$alt_quality_gpi[i] <- alt_quality_gpi
      big_pi$alt_position_gpi[i] <- alt_position_gpi
      big_pi$alt_total_score[i] <- alt_total_score
      big_pi$alt_profile_score[i] <- alt_profile_score
      big_pi$alt_profile_i_score[i] <- alt_profile_i_score
      if (verbose == T){
        print(paste(i, "-", " Total Score:", total_score))
        flush.console()
      }
      Sys.sleep(sleep)
    }
    return(big_pi)
  }
  if (!short.out){
    for (i in 1:n){
      if (length(grep(aa_regex, sequence[i])) >= 1){
        warning(paste("sequence", "[", i, "]",
                      " contains symbols not corresponding to amino acids",
                      sep = ""), call. = F)
        impro <- paste("sequence", "[", i, "]",
                    " contains symbols not corresponding to amino acids",
                    sep = "")
      } else {
        if (nchar(sequence[i]) >= 55){
          form <- rvest::set_values(form, Sequence = sequence[i])
          form_res <- suppressMessages(rvest::submit_form(session, form))
          resulti <- rvest::html_text(rvest::html_nodes(form_res, "pre"))
          impro <- unlist(strsplit(resulti, "\\\n"))
          impro <- impro[grep(crop_1, impro):length(impro)]
          score_print <- extract_val("Total Score\\.",
                                  impro)[1]
          if (verbose == T){
            print(paste(i, "-", " Total Score:", score_print))
            flush.console()
          }
        } else {
          impro <- paste("Your sequence length is too short for the",
                      "GPI prediction algorithm! As a minimum,",
                      "55 amino acids are required.")
        }
      }
      big_pi_list[[i]] <- NaN * seq(impro[impro != ""])
      big_pi_list[[i]] <- impro[impro != ""]
      Sys.sleep(sleep)
    }
    return(big_pi_list)
  }
}
