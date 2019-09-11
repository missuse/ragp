
<!-- README.md is generated from README.Rmd. Please edit that file -->
ragp - hydroxyproline aware filtering of hydroxyproline rich glycoprotein sequences
===================================================================================

[![Build Status](https://travis-ci.org/missuse/ragp.svg?branch=master)](https://travis-ci.org/missuse/ragp) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

`ragp` is an R package primarily designed for mining and analysis of plant hydroxyproline rich glycoproteins. It incorporates a novel concept with an additional analysis layer where the probability of proline hydroxylation is estimated by a machine learning model. Only proteins predicted to contain hydroxyprolines are further analysed for HRGP characteristic motifs and features. `ragp` can also be used for protein annotation by obtaining predictions for several protein features based on sequence (secretory signals, transmembrane regions, domains, glycosylphosphatidylinositol attachment sites and disordered regions). Additionally ragp provides tools for visualization of the mentioned attributes.

Short example:

``` r
library(ragp)
ids <- c("Q9FLL2", #several uniprot accessions
         "Q9LS14",
         "Q9S7I8",
         "Q9M2Z2",
         "Q9FIN5")

seqs <- unlist(protr::getUniProt(ids)) #download sequences 

p1 <- plot_prot(seqs, #plot sequence features
                ids,
                hyp = FALSE, #do not plot hydroxyprolines
                ag = FALSE) #do not plot ag spans
p1
```

![](https://missuse.github.io/ragp/reference/figures/README-plot_prot-2.svg)

Installation
------------

You can install ragp from github with:

``` r
# install.packages("devtools")
devtools::install_github("missuse/ragp")
```

Or alternatively

``` r
# install.packages("devtools")
devtools::install_github("missuse/ragp",
                         build_vignettes = TRUE)
```

to build vignettes which can be viewed by:

``` r
browseVignettes("ragp")
```

Tutorials
---------

Tutorials on usage of `ragp` functions with examples on how to combine them into meaningful HRGP filtering and analysis pipelines are available at: <https://missuse.github.io/ragp/>

Bug reports
-----------

If you encounter undesired behavior in `ragp` functions or you have ideas how to improve them please open an issue at: <https://github.com/missuse/ragp/issues>

Citation
--------

If you find `ragp` useful in your own research please cite our Glycobiology [paper](https://academic.oup.com/glycob/advance-article-abstract/doi/10.1093/glycob/cwz072/5567434) (Preprint).

> Milan B Dragićević, Danijela M Paunović, Milica D Bogdanović, Slađana I Todorović, Ana D Simonović, ragp: Pipeline for mining of plant hydroxyproline-rich glycoproteins with implementation in R, Glycobiology, cwz072, <https://doi.org/10.1093/glycob/cwz072>

You can get citation info via `citation("ragp")` or by copying the following BibTex entry:

``` bibtex
@article{10.1093/glycob/cwz072,
    author = {Dragićević, Milan B and Paunović, Danijela M and Bogdanović, Milica D and Todorović, Slađana I and Simonović, Ana D},
    title = "{ragp: Pipeline for mining of plant hydroxyproline-rich glycoproteins with implementation in R}",
    journal = {Glycobiology},
    year = {2019},
    month = {09},
    issn = {1460-2423},
    doi = {10.1093/glycob/cwz072},
    url = {https://doi.org/10.1093/glycob/cwz072}
}
```

Acknowledgements
----------------

This software was developed with funding from the Ministry of Education, Science and Technological Development of the Republic of Serbia (Projects TR31019 and OI173024).
