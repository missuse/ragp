ragp 0.3.5
===============

New Features
------------

* added new function `get_signalp5()` which queries SignalP5 web server (http://www.cbs.dtu.dk/services/SignalP)  
* added new function `get_tmhmm()` which queries TMHMM  v. 2.0 web server (http://www.cbs.dtu.dk/services/TMHMM/)  
* `plot_prot()` `nsp` argument can now be `"signalp"`, `"signalp5"` or `"none"`. Default is `"signalp5"`. This argument determines if `get_signalp()` or `get_signalp5()` are used for N-sp prediction. 
* `plot_prot()` `tm` argument can now be `"phobius"`, `"tmhmm"` or `"none"`. Default is `"phobius"`. This argument determines if `get_phobius()` or `get_tmhmm()` are used for TM prediction.  
* all `get_*` and `scan_*` functions, as well as `maab()` now work with `AAStringSet` class objects. #5

Bug Fixes and Improvements
--------------------------

* `get_signalp()` output contains an additional column `sp.length` - integer, length of the predicted signal peptide. This column is a copy of `Ymax.pos`. **Potentially BREAKING CHANGE**  
* `get_phobius()` output contains an additional column `sp.length` - integer, length of the predicted signal peptide. **Potentially BREAKING CHANGE**
* `get_phobius()` the name of the output column "`Name`" has been changed to "`id`". **BREAKING CHANGE**
* `plot_prot()` default value for `gpi` argument has been changed to "`netgpi`". 


ragp 0.3.2.0002
===============

Bug Fixes and Improvements
--------------------------

* `get_netGPI()` now now runs one job at a time. 

* `get_netGPI()` splitter argument default value has been changed to 2500.

* `get_targetp()` has been fixed to work with `org_type = "non_plant"`.  

* `get_espritz()` has been restored since Espritz server is again available at http://old.protein.bio.unipd.it/espritz/.  

* `plot_prot()` has new arguments `gpi_size` - controls the size of the gpi symbol and `gpi_shape`- controls the shape of the gpi symbol.  

* `plot_prot()` has new argument `hyp_scan` - which is a logical and determines if `scan_ag()` (when `ag = TRUE`) should scan for arabinogalactan motifs containing only predicted hydroxyprolines. This argument changes the previous default behavior of `plot_prot()` which was equivalent to `hyp_scan = FALSE`.  

* `get_phobius()` and `get_big_pi()` now use https web server addresses. #6  


ragp 0.3.2.0001
===============

Bug Fixes and Improvements
--------------------------

* `get_espritz()` has been removed since Espritz server is no longer available for obtaining predictions. **BREAKING CHANGE**.
* xgboost models used by ragp were re-saved using xgboost 1.1.1.1 to increase compatibility.
* xgboost models used by ragp are now stored in the inst directory instead of internal sysdata.rda.


ragp 0.3.2.0000
===============

New Features
------------

* `plot_prot()` arguments `nsp`, `domain`, `tm`, `gpi` and `disorder` can also be user supplied data frames obtained by calling `get_signalp()`, `get_hmm()`, `get_phobius()`, `get_big_pi()`, `get_pred_gpi()`, `get_netGPI()` and `get_espritz()` on the same sequences supplied to `plot_prot()`.  


Bug Fixes and Improvements
--------------------------

* `get_big_pi()` output column `is.bigpi` has been renamed to `is.gpi` to bring the output in line with other gpi predicting functions. **BREAKING CHANGE**.  
* `get_big_pi()` output column order has been changed to:  `id`, `is.gpi`, `Quality`, `omega_site` and `PValue`. **BREAKING CHANGE**.


ragp 0.3.1.0000
===============

New Features
------------

*  New `get_netGPI()` queries NetGPI web server (https://services.healthtech.dtu.dk/service.php?NetGPI-1.0) for predictions of GPI-anchored proteins.  

* `maab()` argument `get_gpi` has an additional option `"netgpi"`. When set, the function will query NetGPI web server to resolve ambiguities in maab classes depending on GPI-anchoring predictions. 

* `plot_prot()` argument `gpi` has an additional option `"netgpi"`. When set, the function will query NetGPI web server for prediction of omega sites.


ragp 0.3.0.0002
===============

Bug Fixes and Improvements
--------------------------

* `maab()` now correctly performs when only a single protein sequence is provided as an argument. 

* `get_hmm()` receives additional numeric arguments: `ievalue` and `bitscore`. These arguments are used to filter sequences with lower or equal `ievalue` and higher or equal `bitscore` in the output. This is useful when used from `plot_prot()` to avoid plotting weakly identified domains. 


ragp 0.3.0.0001
===============

Bug Fixes and Improvements
--------------------------

* `get_big_pi()` and `get_pred_gpi()` now return `NA` in respective `is.bigpi` and `is.gpi` columns when the servers are unable to make a prediction due to non-amino acid letters or length of the sequence.

* `maab()` now correctly does not resolve class ambiguities when the logical vector provided as `gpi` argument contains `NA` values. Previously it returned `NA` as `maab_class`. This is useful when `maab()` is called with `get_gpi = 'predgpi'` or `get_gpi = 'bigpi'` arguments, and the corresponding servers are unable to make a prediction for a sequence due to non-amino acid letters or length of the sequence.

ragp 0.3.0.0000
===============

New Features
------------

* `predict_hyp()` internal model is updated to 2nd version ('V2'). Predictions are around 25% faster compared to the first version. The performance in terms of accuracy is similar based on the test set used. 'V2' was created using a more streamlined manner and is the default model. The old model ('V1') is still available using the `version` argument to `predict_hyp()`.

Bug Fixes and Improvements
--------------------------

* `predict_hyp()` now checks if all provided ids are unique. Previously non unique ids caused an error at the end of computation. 

* `predict_hyp()` `sequence` output has changed for sequences containing non amino acid letters. Previously NA was returned for such sequences. At present all "P"" for which the probability is higher then the defined threshold (`tprob` argument) are changed to "O"" and all others are left as "P".

ragp 0.2.1.9100
===============

Bug Fixes and Improvements
--------------------------

* `maab()` now correctly outputs when there are no MAAB classes found and `get_gpi` argument is set to `"predgpi"` or `"bigpi"`. Previously this caused an error due to `stringr::write.fasta()` attempting to write a file with no sequences.

* `pfam2go()` now takes Pfam > GO mappings from ftp://ftp.geneontology.org/pub/go/external2go/pfam2go instead of http://geneontology.org/external2go/pfam2go.


ragp 0.2.1.9000
===============

Bug Fixes and Improvements
--------------------------

* `get_phobius()`, `get_big_pi()`, `get_pred_gpi()`, `maab()` and `plot_prot()`, gain an additional argument `progress`. 
`progress` is a logical value determining whether to show the progress bar, (default `FALSE`).

* `get_targetp()`, `get_signalp()` and `get_hmm()` argument `progress` is now set to `FALSE` at default.

* `get_hmm()` argument `verbose` is now set to `FALSE` at default.

* added `pkgdown` site for `ragp` at: https://missuse.github.io/ragp/.

ragp 0.2.1.0000
===============

Bug Fixes and Improvements
--------------------------

* `get_targetp()` and `get_signalp()` gain additional arguments `progress` and `attempts`. 
`progress` is a logical value determining whether to show the progress bar, (default `TRUE`). 
`attempts` is an integer value determining the number of repeated attempts if server unresponsive (default 2).
These functions now return finished queues if server becomes unresponsive.  
  
  
ragp 0.2.0.0000
===============

New Features
------------

* `get_big_pi()`, `get_espritz()`,`get_hmm()`, `get_phobius()`, `get_pred_gpi()`, `get_signalp()`,
`get_targetp()`, `maab()`, `predict_hyp()`, `scan_ag()` and `scan_nglc()` now use the S3 object
system which will make further extensions of accepted inputs straightforward.

Bug Fixes and Improvements
--------------------------

* Completely rewrote the code for `scan_ag()` which is now simplified and easier to read.
The function is now about 20% slower.

* fixed a bug in `scan_ag()` when argument `exclude_ext = "all"` which resulted
in detection of unwanted amino acids in certain sequence arrangements.
The bug occurred only if AG glycomoduls were detected, it did
not introduce any. If you used this function with `exclude_ext = "all"`
it is advisable to rerun these analyses again. 

* removed deprecated functions: `get_signap_file()`, `get_targetp_file()`,
`get_phobius_file()`.


ragp 0.1.0.0004
===============

New Features
------------

* New `get_espritz()` queries ESpritz web server for predictions on protein
disordered regions.
  
* `plot_prot()` gains an additional argument `disorder`, logical indicating
should the predicted disordered regions be plotted. Defaults to `FALSE`.
  
Bug Fixes and Improvements
--------------------------

* fixed a bug in function `plot_prot()` introduced in 0.1.0.0003 which prevented 
  GPIs to be plotted when `gpi = "bigpi"`.


ragp 0.1.0.0003
===============

New Features
------------

* New `get_pred_gpi()` queries PredGPI web server for predictions on GPI presence
and omega site location.
  
Bug Fixes and Improvements
--------------------------
  
* `maab()` argument `get_gpi` has been changed and now accepts strings
  as input: `get_gpi = c("bigpi", "predgpi", "none")` which indicate whether
  to query Big Pi or PredGPI server or not to resolve class ambiguities.
   
* `plot_prot()` argument `gpi` has been changed and now accepts strings
  as input `gpi = c("bigpi", "predgpi", "none")` which indicate whether
  to query Big Pi or PredGPI server or not to plot GPI positions.
  
* fixed a bug in `plot_prot()` which caused extra-cellular regions to start at 0
  instead of 1.

* fixed a bug in `get_big_pi()` when shorter than 55 amino acid sequences 
  were provided.

  
ragp 0.1.0.0002
===============

Bug Fixes and Improvements
--------------------------

* Removed rvest dependency


ragp 0.1.0.0001
===============

New Features
------------

* Added vignette. 

Bug Fixes and Improvements
--------------------------

* When using `plot_prot()` with argument `dom_sort = "ievalue"` the domains with the
  lowest independent e-value will now be correctly plotted on top.
  
* The ratio of x and y axes in `plot_prot()` output is now calculated based on sequence
  length and should provide more consistent diagrams.
  
* `get_hmm()` can now handle protein sequences of arbitrary length by splitting them
  into several shorter overlapping sequences and querying hmmscan. 
 
 
ragp 0.1.0.0000
===============

New Features
------------

* New `plot_prot()` returns a `ggplot2` diagram of protein structure based on `hmmscan`
domain annotation and several types of predictions.

* New `get_signalp()`, `get_targetp()` and `get_phobius()` replace
deprecated  `get_signalp_file()`, `get_targetp_file()` and `get_phobius_file()`.
These functions accept R objects as well as FASTA files as input (#5).

* New `plot_signalp()` takes one single letter protein sequence and returns a detailed
SignalP prediction along with a plot.

* New `scan_nglc()` detects motifs for N-glycosylation on asparagine residues.
  
* `maab()` gain an additional logical argument `get_gpi` (default `FALSE`);
if `TRUE` `get_big_pi` will be called on all sequences belonging to one of the
HRGP classes thus resolving class ambiguities that depend on GPI knowledge.
  
* Added a `NEWS.md` file to track changes to the package.

Bug Fixes and Improvements
--------------------------

* Significantly improved the speed of `get_big_pi()`. Removed the verbose argument.
  Added a progress bar (#3).

* `get_hmm()` gains new arguments: `timeout` - time in seconds to wait for
  the server response (default is 10s). `attempts` - number of attempts if
  server unresponsive (default is 2 times) (#3).

* `get_signalp()` and `get_targetp()` now perform 10 parallel queries instead
  of arbitrary many. This results in higher stability when running on many
  sequences, and a slight reduction in speed. These functions now have a 
  progress bar.
  
* `get_phobius()` gains a progress bar.

* Changed internal behavior when `trunc` argument is specified for `get_signalp()`.
  Now the truncation is performed prior to sending the files to the server
  resulting in higher efficiency.

* `scan_ag()` gains logical `tidy` argument. The specific output is available when
`simplify = FALSE` and `tidy = TRUE` in function call (#1).