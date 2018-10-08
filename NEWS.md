NEWS
-----------------------------
  
version 0.2.1.0000

BUG FIXES AND IMPROVEMENTS

* Additional arguments added to `get_targetp` and `get_signalp` functions: `progress` and `attempts`. 
`progress` is a logical value determining whether to show the progress bar, at default set to TRUE. 
`attempts` is an integer value determining the number of repeated attempts if server unresponsive, at default set to 2.
These functions now return finished queues if server becomes unresponsive.  
  
version 0.2.0.0000
-----------------------------
NEW FEATURES

* All functions that operate on multiple input types: `get_big_pi`, `get_espritz`,
`get_hmm`, `get_phobius`, `get_pred_gpi`, `get_signalp`, `get_targetp`, `maab`,
`predict_hyp`, `scan_ag` and `scan_nglc` now do so via the S3 object system which will 
make further extensions of accepted inputs straightforward.

BUG FIXES AND IMPROVEMENTS

* Completely rewrote the code for `scan_ag` which is now simplified and easier to read.
The function is now about 20% slower.

* fixed a bug in `scan_ag` when argument `exclude_ext = "all"` which resulted
in detection of unwanted amino acids in certain sequence arrangements.
The bug occurred only if AG glycomoduls were detected, it did
not introduce any. If you used this function with `exclude_ext = "all"`
it is advisable to rerun these analyses again. 

* removed deprecated functions: `get_signap_file`, `get_targetp_file`,
`get_phobius_file`


version 0.1.0.0004
-----------------------------
NEW FEATURES

* Added `get_espritz` function which queries ESpritz web server for 
  predictions on protein disordered regions.
  
* `plot_prot` function has an additional argument `disorder`, indicating
  should the predicted disordered regions be plotted. Defaults to FALSE.
  
BUG FIXES AND IMPROVEMENTS

* fixed a bug in function `plot_prot` introduced in 0.1.0.0003 which prevented 
  GPIs to be plotted when `gpi = "bigpi"`.


version 0.1.0.0003
-----------------------------
NEW FEATURES

* Added `get_pred_gpi` function which queries PredGPI web server for 
  predictions on GPI presence and omega site location.
  
BUG FIXES AND IMPROVEMENTS
  
* `maab` function argument `get_gpi` has been changed and now accepts strings
  as input: `get_gpi = c("bigpi", "predgpi", "none")` which indicate whether
  to query Big Pi or PredGPI server or not to resolve class ambiguities.
   
* `plot_prot` function argument `gpi` has been changed and now accepts strings
  as input `gpi = c("bigpi", "predgpi", "none")` which indicate whether
  to query Big Pi or PredGPI server or not to plot GPI positions.
  
* fixed a bug in `plot_prot` which caused extracellular regions to start at 0
  instead of 1.

* fixed a bug in `get_big_pi` when shorter than 55 amino acid sequnces 
  were provided.

  
version 0.1.0.0002
-----------------------------
BUG FIXES AND IMPROVEMENTS

* Removed rvest dependency


version 0.1.0.0001
-----------------------------
NEW FEATURES

* Added vignette. 

BUG FIXES AND IMPROVEMENTS

* When using plot_prot with argument dom_sort = "ievalue" the domains with the
  lowest independent e-value will now be correctly plotted on top.
  
* The ratio of x and y axes in plot_prot output is now calculated based on sequence
  length and should provide more consistent diagrams.
  
* get_hmm can now handle protein sequences of arbitrary length by splitting them
  into several shorter overlapping sequences and querying hmmscan. 
 
 
version 0.1.0.0000
-----------------------------

NEW FEATURES

* Added `plot_prot` function which returns a `ggplot2` diagram of protein
  structure based on `hmmscan` domain annotation and several types of
  predictions.

* Added `get_signalp`, `get_targetp` and `get_phobius` which replace
  `get_signalp_file`, `get_targetp_file` and `get_phobius_file`.
  These functions accept R objects as well as FASTA files as input. (#5)

* Added `plot_signalp` function which takes one single letter protein
  sequence and returns a detailed SignalP prediction along with a plot.

* Added `scan_nglc` function which detects motifs for N-glycosylation on
  asparagine residues.
  
* `maab` function now has an additional Bolean argument `get_gpi`
  (default `FALSE`); if `TRUE` `get_big_pi` will be called on all sequences
  that belong to one of the HRGP classes thus resolving class ambiguities
  that depend on GPI knowledge.
  
* Added a `NEWS.md` file to track changes to the package.

BUG FIXES AND IMPROVEMENTS

* Significantly improved the speed of `get_big_pi`. Removed the verbose argument.
  Added a progress bar. (#3)

* Added new arguments to `get_hmm`: `timeout` - time in seconds to wait for
  the server response (default = 10s). `attempts` - number of attempts if
  server unresponsive (default = 2 times). (#3)

* `get_signalp` and `get_targetp` now perform 10 parallel queries instead
  of arbitrary many. This results in higher stability when running on many
  sequences, and a slight reduction in speed. These functions now have a 
  progress bar.
  
* Added progress bar to `get_phobius`.

* Changed internal behavior when `trunc` argument is specified for `get_signalp`.
  Now the truncation is performed prior to sending the files to the server
  resulting in higher efficiency.

* Added `tidy` output for `scan_ag`. It is available when `simplify = FALSE` and
  `tidy = TRUE` in function call. (#1)