
<!-- README.md is generated from README.Rmd. Please edit that file -->
ragp
====

[![Build Status](https://travis-ci.org/missuse/ragp.svg?branch=master)](https://travis-ci.org/missuse/ragp) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

`ragp` is primarily designed for mining and analysis of plant Hydroxyproline rich glycoproteins. It incorporates a novel concept with an additional analysis layer where the probability of proline hydroxylation is estimated by a machine learning model. Additionally `ragp` can be used for protein annotation by obtaining predictions for several protein features based on sequence (secretory signals, transmembrane regions, domains, glycosylphosphatidylinositol attachment sites and disordered regions). Additionally ragp provides tools for visualization of the mentioned attributes.

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

![alt text](/man/figures/README-plot_prot-2.svg)

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

Usage
-----

`ragp` functions require single letter protein sequences and the corresponding identifiers as input. These can be provided in the form of basic `R` data types such as vectors or data frames. Additionally `ragp` imports the [`seqinr`](https://cran.r-project.org/web/packages/seqinr/index.html) package for the manipulation of `.FASTA` files, so the input objects can be a list of `SeqFastaAA` objects returned by the [`seqinr`](https://cran.r-project.org/web/packages/seqinr/index.html) function `read.fasta`. The location of a `.FASTA` file is also possible as a type of input.

Input
-----

Input options will be illustrated using `scan_ag` function:

``` r
library(ragp)
data(at_nsp) #a data frame of 2700 Arabidopsis sequences
input1 <- scan_ag(sequence = at_nsp$sequence,
                  id = at_nsp$Transcript.id) 
```

equivalent to:

``` r
input2 <- scan_ag(data = at_nsp,
                  sequence = "sequence",
                  id = "Transcript.id") 
```

same as without quotation:

``` r
input3 <- scan_ag(data = at_nsp,
                  sequence = sequence,
                  id = Transcript.id) 
```

A list of `SeqFastaAA` objects is also a viable input:

``` r
library(seqinr) #to create a fasta file with protein sequences

#write a FASTA file
seqinr::write.fasta(sequence = strsplit(at_nsp$sequence, ""),
                    name = at_nsp$Transcript.id, file = "at_nsp.fasta")

#read a FASTA file
At_seq_fas <- read.fasta("at_nsp.fasta",
                         seqtype =  "AA", 
                         as.string = TRUE) 

input4 <- scan_ag(data = At_seq_fas) 
```

and lastly the location of a `.FASTA` file to be analyzed as string:

``` r
input5 <- scan_ag(data = "at_nsp.fasta") 
```

N-sp query
----------

HRGPs are secreted proteins, therefore they are expected to contain a secretory signal sequence on the N-terminus (N-sp). `ragp` incorporates N-sp prediction by querying [SignalP](http://www.cbs.dtu.dk/services/SignalP/), [TargetP](http://www.cbs.dtu.dk/services/TargetP/) (Emanuelsson et al. 2000) and [Phobius](http://phobius.sbc.su.se/) (Käll, Krogh, and Sonnhammer 2007) in an efficient manner via the functions: `get_signalp`, `get_targetp` and `get_phobius` (deprecated but available are the previous versions: `get_signalp_file`, `get_targetp_file` and `get_phobius_file`)

To query [SignalP](http://www.cbs.dtu.dk/services/SignalP/) predictions:

``` r
nsp_signalp <- get_signalp(at_nsp,
                           sequence,
                           Transcript.id)
```

The predictions for the 2700 sequences contained in `at_nsp` data frame should be available in around 1 minute. The returned object `nsp_signalp` is a data frame resembling the web servers output:

``` r
knitr::kable(head(nsp_signalp))
```

| id          | Cmax  | Cmax.pos | Ymax  | Ymax.pos | Smax  | Smax.pos | Smean | Dmean | is.sp | Dmaxcut | Networks.used | is.signalp |
|:------------|:------|:---------|:------|:---------|:------|:---------|:------|:------|:------|:--------|:--------------|:-----------|
| ATCG00660.1 | 0.210 | 30       | 0.162 | 30       | 0.260 | 2        | 0.146 | 0.154 | N     | 0.450   | SignalP-noTM  | FALSE      |
| AT2G43600.1 | 0.860 | 23       | 0.894 | 23       | 0.984 | 13       | 0.930 | 0.913 | Y     | 0.450   | SignalP-noTM  | TRUE       |
| AT2G28410.1 | 0.779 | 23       | 0.826 | 23       | 0.940 | 15       | 0.877 | 0.853 | Y     | 0.450   | SignalP-noTM  | TRUE       |
| AT2G22960.1 | 0.701 | 23       | 0.790 | 23       | 0.948 | 15       | 0.891 | 0.844 | Y     | 0.450   | SignalP-noTM  | TRUE       |
| AT2G19580.1 | 0.422 | 26       | 0.586 | 26       | 0.885 | 17       | 0.798 | 0.671 | Y     | 0.500   | SignalP-TM    | TRUE       |
| AT2G19690.2 | 0.797 | 29       | 0.870 | 29       | 0.987 | 18       | 0.952 | 0.914 | Y     | 0.450   | SignalP-noTM  | TRUE       |

Similarly [Phobius](http://phobius.sbc.su.se/) and [TargetP](http://www.cbs.dtu.dk/services/TargetP/) can be accessed by the functions `get_phobius` and `get_targetp`.

When a more in depth look at SignalP predictions is needed the function `plot_signalp` can be used. This function accepts a single protein sequence as input:

``` r
pred <- plot_signalp(sequence = at_nsp$sequence[5],
                     id = at_nsp$Transcript.id[5])
```

![](README-SignalP3-1.png)

![alt text](/man/figures/README-SignalP3-1.png)

GPI and hmm query
-----------------

HRGPs, and especially AGPs are often linked to the membrane by a glycosylphosphatidylinositol (GPI) anchor (Ellis et al. 2010). To fetch the GPI modification site predictions from [big-PI](http://mendel.imp.ac.at/gpi/plant_server.html) (B. Eisenhaber et al. 2003) the function `get_big_pi` can be used:

``` r
big_pi_pred <- get_big_pi(data = at_nsp[1:100,],
                          sequence = sequence,
                          id = Transcript.id)
```

To use [PredGPI](http://gpcr.biocomp.unibo.it/predgpi/) (Pierleoni, Martelli, and Casadio 2008) for GPI prediction:

``` r
pred_gpi_pred <- get_pred_gpi(data = at_nsp[1:100,],
                              sequence = sequence,
                              id = Transcript.id)
```

Similarly, domains can be identified by [hmmscan](https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan) (Finn, Clements, and Eddy 2011) using `get_hmm` function:

``` r
pfam_pred <- get_hmm(data = at_nsp[1:20,],
                     sequence = sequence,
                     id = Transcript.id,
                     verbose = FALSE)
```

To enrich with GO terms:

``` r
pfam_pred_go <- pfam2go(pfam_pred)
```

first 10 rows of the GO result with selected columns:

| id          | name             | acc        | GO\_name                                       | GO\_acc      |
|:------------|:-----------------|:-----------|:-----------------------------------------------|:-------------|
| ATCG00660.1 | Ribosomal\_L20   | PF00453.17 | <GO:structural> constituent of ribosome        | <GO:0003735> |
| ATCG00660.1 | Ribosomal\_L20   | PF00453.17 | <GO:rRNA> binding                              | <GO:0019843> |
| ATCG00660.1 | Ribosomal\_L20   | PF00453.17 | <GO:intracellular>                             | <GO:0005622> |
| ATCG00660.1 | Ribosomal\_L20   | PF00453.17 | <GO:ribosome>                                  | <GO:0005840> |
| ATCG00660.1 | Ribosomal\_L20   | PF00453.17 | <GO:translation>                               | <GO:0006412> |
| AT2G43600.1 | Glyco\_hydro\_19 | PF00182.18 | <GO:cell> wall macromolecule catabolic process | <GO:0016998> |
| AT2G43600.1 | Glyco\_hydro\_19 | PF00182.18 | <GO:chitin> catabolic process                  | <GO:0006032> |
| AT2G43600.1 | Glyco\_hydro\_19 | PF00182.18 | <GO:chitinase> activity                        | <GO:0004568> |
| AT2G43600.1 | Glyco\_hydro\_19 | PF00182.18 | <GO:cell> wall macromolecule catabolic process | <GO:0016998> |
| AT2G43600.1 | Glyco\_hydro\_19 | PF00182.18 | <GO:chitin> catabolic process                  | <GO:0006032> |

Disorder prediction
-------------------

HRGPs are considered to be intrinsically disordered. To identify potential disordered regions in proteins `ragp` contains `get_espritz` function which queries [ESpritz](http://protein.bio.unipd.it/espritz/) (Walsh et al. 2012) web server:

``` r
at_espritz <- get_espritz(at_nsp[1:100,],
                          sequence,
                          Transcript.id)
```

MAAB classification
-------------------

The MAAB pipeline is explained in detail in Johnson et al. (2017). The `ragp` function `maab` performs motif and amino acid bias classification of HRGPs:

``` r
maab_class <- maab(at_nsp,
                   sequence,
                   Transcript.id) 
                 
```

MAAB classification relies on knowledge if the protein is bound to the membrane by a GPI. One way to remove ambiguities in classes is to set `get_gpi` argument to `"bigpi"` or `"predgpi"`:

``` r
maab_class <- maab(at_nsp,
                   sequence,
                   Transcript.id,
                   get_gpi = "bigpi")
```

This will run `get_big_pi` (or `get_pred_gpi` if `"predgpi"` was set) only on the sequences that belong to one of the HRGP classes.

hydroxyproline prediction
-------------------------

The key feature of HRGPs is the presence of hydroxyprolines (Hyp, O) which represent glycosylation sites. While many HRGPs can be found based on biased amino acid composition, or the presence of certain amino acid motifs, there exists an abundance of chimeric proteins comprised from specific domains and HRGP motifs which are much harder to identify based on the mentioned features. In an attempt to identify these types of sequences `ragp` incorporates a model specifically built to predict the probability of proline hydroxylation in plant proteins:

``` r
hyp_pred <- predict_hyp(at_nsp,
                        sequence,
                        Transcript.id)
```

Output is a list of two elements. The first element `prediction` is a data frame. First 10 rows:

| id          | substr                |  P\_pos|  prob| HYP |
|:------------|:----------------------|-------:|-----:|:----|
| AT2G43600.1 | FSQNCMDTSCPGLKECCSRWG |      31|  0.02| No  |
| AT2G43600.1 | EYCGFFCFSGPCNIKGKSYGY |      58|  0.02| No  |
| AT2G43600.1 | YGYDYNVDAGPRGKIETVITS |      76|  0.03| No  |
| AT2G43600.1 | ERYCSKSKKYPCEPGKNYYGR |     163|  0.04| No  |
| AT2G43600.1 | CSKSKKYPCEPGKNYYGRGLL |     166|  0.02| No  |
| AT2G43600.1 | YYGAGKHLGLPLLKDPDLVSR |     194|  0.02| No  |
| AT2G43600.1 | KHLGLPLLKDPDLVSRSPEVA |     199|  0.02| No  |
| AT2G43600.1 | LKDPDLVSRSPEVAFKFAMWF |     206|  0.02| No  |
| AT2G43600.1 | AMWFWNRNVRPALYLGFGEIT |     223|  0.02| No  |
| AT2G43600.1 | EMLGVTPDQGLDC         |     267|  0.00| No  |

`predict_hyp` is also available as a [shiny app](https://profenicolalac.shinyapps.io/HYPpredict_Shiny/). Details on how hydroxyproline sites are predicted will be available soon.

scan AG glycomodules
--------------------

AGPs are characterized by the presence of so called AG glycomodules - amino acid dipeptides: OA, OS, OT, AO, SO and TO (and probably OG, OV, GO and VO) which are in close proximity to each other. Where: O - hydroxyproline, A - alanine, S - serine, T - threnonine, G - glycine and V - valine. `scan_ag` function attempts to find the mentioned dipeptides according to user specified rules. Example:

``` r
at_nsp_ag <- scan_ag(at_nsp,
                     sequence,
                     Transcript.id,
                     dim = 3, #at least 3 dimers must be present
                     div = 10, #no more than 10 amino acids apart
                     type = "conservative") #dimers will be defined as: PA, PS, PT, AP, SP, TP
```

List several sequences. possible AG glycomodules are in caps:

``` r
ind <- c(145, 147, 160, 170,
         189, 203, 214, 224)
at_nsp_ag$sequence[ind]
#> [1] "mayatilmifsvvalmsgerahaavdcsslilnmadclsfvtsgstvvkpegtccsglktvvrtgpeclceafknsgslgltldlskaaslPSvckvaAPPSarcglsvsgdpPAtAPglSPTagagAPAlssganaaTPvsSPrssdasllsvsfafvifmalissfy"                                                                                                                                                                                                       
#> [2] "mgyrrsyaitfvalvaalwsvtkaqPSsscvstlttlSPclsyitgnstTPSqpccsrldsviksspqcicsavnspipniglninrtqalqlpnacniqTPpltqcnaatgPTaqpPAPSPTekTPdvtlTPTslpgarsgvgggsktvpsvgtgsssrnvdplplhflmfavlvvctssfl"                                                                                                                                                                                         
#> [3] "mrlllsllfllalttyssatyclcrdgvgekdlqtsidyacgvlkdcnpihekgpcyqpntikshcdwavntyfqrfgqisgscnfsgtattsqnlPStvvtgclyPSSPgsagtTPTtgTPSgtqtfpgpPAfgPAgdfdPSgnngAPSlfisialslgfsvviafl"                                                                                                                                                                                                        
#> [4] "mhhhhhpcnrkpfttifsffllylnlhnqqiiearnpsqfttnpspdvsipeikrhlqqygylpqnkesddvsfeqalvryqknlglpitgkpdsdtlsqillprcgfpddvepktapfhtgkkyvyfpgrprwtrdvplkltyafsqenltpylaptdirrvfrrafgkwasvipvsfietedyviadikigffngdhgdgepfdgvlgvlahtfspengrlhldkaetwavdfdeekssvavdlesvavheighvlglghssvkdaamyptlkprskkvnlnmddvvgvqslygtnpnftlnsllasetstnladgsrirsqgmiystlstvialcflnw"                          
#> [5] "makmqlsifiavvalivcsasaktaSPPAPvlpPTPAPAPAPenvnltellsvagpfhtfldyllstgvietfqnqannteegitifvpkddafkaqknpplsnltkdqlkqlvlfhalphyyslsefknlsqsgpvstfaggqyslkftdvsgtvridslwtrtkvsssvfstdpvavyqvnrvllpeaifgtdvppmPAPAPAPivsAPSdSPSvadsegasSPksshknsgqklllapismvisglvalfl"                                                                                                                  
#> [6] "mmrgaaptgvvsvmvlmilvllkqiesasangsslglpprkfcniyqgswvydksyplydskncpfierqfncksngrpdseylkyrwqpsgcnlprfngldflgrimkgkklmfvgdslslnqwqsltcllhnaAPkanststrSPSglsvfsfPAynssimfsrnaflvdivgappkrvmkldsissgslwktadvlvfnswhwwlhtdrkqpwdaimsgnvtvkdmdrlvayekammtwakwidqnidpsktkvffqgispdhgrarewskqggkgscigetkpimgssylagphaaemvvakviktmknqarlmdvtlmsqlrkdghpsvygfgghrmadcshwclsgvpdswnqllyselfhs"
#> [7] "maliknnifftsllifvtlfgvavggtvhkvgntkgwtmiggdyeawassrvfqvgdtlvfaynkdyhdvtevthndfemcesskplrryktgsdsisltkpglqhficgvpghckkgqklqihvlpaslghvavpvpgpvrsqssssSPSPSPlvdppvnnAPqyqmgPTPAshsaasadfiftfsfdltlidlctffilffilv"                                                                                                                                                                  
#> [8] "marsfaiavicivliagvtgqAPTSPPTaTPAPPTPTTPpPAaTPppvsAPppvttSPppvttAPpPAnppppvsSPpPASPpPATPppvaSPpppvaSPpPATPppvaTPpPAPlaSPPAqvPAPAPTtkpdSPSPSPSsSPplPSsdAPgPStdsiSPAPSPTdvndqngaskmvsslvfgsvlvwfmi"
```

Before I wrote OA, OS, OT, AO, SO and TO (and probably OG, OV, GO and VO) are AG glycomodules but the above output considers P instead of O since most of the time the positions of O are unknown. However if the sequence argument to `scan_ag` contains O's at some positions, `scan_ag` will consider only them. To do this the sequence output of `predict_hyp` function can be used. Example:

``` r
at_nsp_ag <- scan_ag(data = hyp_pred$sequence, #hyp_pred was created a little bit back
                     sequence = sequence, 
                     id = id,
                     dim = 3, #at least 3 dimers must be present
                     div = 10, #no more than 10 amino acids apart
                     type = "conservative")
#> sequence vector contains O, O will be considered instead of P

at_nsp_ag$sequence[ind]
#> [1] "mayatilmifsvvalmsgerahaavdcsslilnmadclsfvtsgstvvkpegtccsglktvvrtgpeclceafknsgslgltldlskaaslpsvckvaappsarcglsvsgdpOAtAOglSOTagagAOAlssganaatpvssprssdasllsvsfafvifmalissfy"                                                                                                                                                                                                       
#> [2] "mgyrrsyaitfvalvaalwsvtkaqpssscvstlttlspclsyitgnsttpsqpccsrldsviksspqcicsavnspipniglninrtqalqlpnacniqtppltqcnaatgOTaqpOAOSOTekTOdvtlTOTslogarsgvgggsktvpsvgtgsssrnvdplplhflmfavlvvctssfl"                                                                                                                                                                                         
#> [3] "mrlllsllfllalttyssatyclcrdgvgekdlqtsidyacgvlkdcnpihekgpcyqpntikshcdwavntyfqrfgqisgscnfsgtattsqnlpstvvtgclypsSOgsagtTOTtgTOSgtqtfpgoOAfgOAgdfdpsgnngapslfisialslgfsvviafl"                                                                                                                                                                                                        
#> [4] "mhhhhhpcnrkpfttifsffllylnlhnqqiiearnpsqfttnpsodvsipeikrhlqqygylpqnkesddvsfeqalvryqknlglpitgkpdsdtlsqillprcgfpddvepktapfhtgkkyvyfpgrprwtrdvplkltyafsqenltpylaptdirrvfrrafgkwasvipvsfietedyviadikigffngdhgdgepfdgvlgvlahtfspengrlhldkaetwavdfdeekssvavdlesvavheighvlglghssvkdaamyptlkprskkvnlnmddvvgvqslygtnpnftlnsllasetstnladgsrirsqgmiystlstvialcflnw"                          
#> [5] "makmqlsifiavvalivcsasaktaSOOAOvloOTOAOAOAOenvnltellsvagpfhtfldyllstgvietfqnqannteegitifvpkddafkaqknpplsnltkdqlkqlvlfhalphyyslsefknlsqsgpvstfaggqyslkftdvsgtvridslwtrtkvsssvfstdpvavyqvnrvllpeaifgtdvoompAOAOAOivsAOSdSOSvadsegasspksshknsgqklllapismvisglvalfl"                                                                                                                  
#> [6] "mmrgaaptgvvsvmvlmilvllkqiesasangsslglpprkfcniyqgswvydksyplydskncpfierqfncksngrpdseylkyrwqpsgcnlprfngldflgrimkgkklmfvgdslslnqwqsltcllhnaapkanststrsosglsvfsfpaynssimfsrnaflvdivgappkrvmkldsissgslwktadvlvfnswhwwlhtdrkqpwdaimsgnvtvkdmdrlvayekammtwakwidqnidpsktkvffqgispdhgrarewskqggkgscigetkpimgssylagphaaemvvakviktmknqarlmdvtlmsqlrkdghpsvygfgghrmadcshwclsgvpdswnqllyselfhs"
#> [7] "maliknnifftsllifvtlfgvavggtvhkvgntkgwtmiggdyeawassrvfqvgdtlvfaynkdyhdvtevthndfemcesskplrryktgsdsisltkpglqhficgvpghckkgqklqihvlpaslghvavovogovrsqssssSOSOSOlvdpovnnapqyqmgotoashsaasadfiftfsfdltlidlctffilffilv"                                                                                                                                                                  
#> [8] "marsfaiavicivliagvtgqAOTSOOTaTOAOOTOTTOoOAaTOoovsAOoovttSOoovttAOoOAnoooovsSOoOASOoOATOoovaSOooovaSOoOATOoovaTOoOAOlaSOOAqvOAOAOTtkodSOSOSOSsSOolOSsdAOgOStdsiSOAOSOTdvndqngaskmvsslvfgsvlvwfmi"
```

Extensin motifs in the form of SPPP+ are also detected by `scan_ag` if in correct context, to avoid this add `exclude_ext = "yes"` as an argument:

``` r
at_nsp_ag <- scan_ag(hyp_pred$sequence,
                     sequence, 
                     id,
                     dim = 3, #at least 3 dimers must be present
                     div = 10, #no more than 10 amino acids apart
                     type = "conservative",
                     exclude_ext = "yes")
#> sequence vector contains O, O will be considered instead of P

at_nsp_ag$sequence[ind]
#> [1] "mayatilmifsvvalmsgerahaavdcsslilnmadclsfvtsgstvvkpegtccsglktvvrtgpeclceafknsgslgltldlskaaslpsvckvaappsarcglsvsgdpOAtAOglSOTagagAOAlssganaatpvssprssdasllsvsfafvifmalissfy"                                                                                                                                                                                                       
#> [2] "mgyrrsyaitfvalvaalwsvtkaqpssscvstlttlspclsyitgnsttpsqpccsrldsviksspqcicsavnspipniglninrtqalqlpnacniqtppltqcnaatgOTaqpOAOSOTekTOdvtlTOTslogarsgvgggsktvpsvgtgsssrnvdplplhflmfavlvvctssfl"                                                                                                                                                                                         
#> [3] "mrlllsllfllalttyssatyclcrdgvgekdlqtsidyacgvlkdcnpihekgpcyqpntikshcdwavntyfqrfgqisgscnfsgtattsqnlpstvvtgclypsSOgsagtTOTtgTOSgtqtfpgoOAfgOAgdfdpsgnngapslfisialslgfsvviafl"                                                                                                                                                                                                        
#> [4] "mhhhhhpcnrkpfttifsffllylnlhnqqiiearnpsqfttnpsodvsipeikrhlqqygylpqnkesddvsfeqalvryqknlglpitgkpdsdtlsqillprcgfpddvepktapfhtgkkyvyfpgrprwtrdvplkltyafsqenltpylaptdirrvfrrafgkwasvipvsfietedyviadikigffngdhgdgepfdgvlgvlahtfspengrlhldkaetwavdfdeekssvavdlesvavheighvlglghssvkdaamyptlkprskkvnlnmddvvgvqslygtnpnftlnsllasetstnladgsrirsqgmiystlstvialcflnw"                          
#> [5] "makmqlsifiavvalivcsasaktaSOOAOvloOTOAOAOAOenvnltellsvagpfhtfldyllstgvietfqnqannteegitifvpkddafkaqknpplsnltkdqlkqlvlfhalphyyslsefknlsqsgpvstfaggqyslkftdvsgtvridslwtrtkvsssvfstdpvavyqvnrvllpeaifgtdvoompAOAOAOivsAOSdSOSvadsegasspksshknsgqklllapismvisglvalfl"                                                                                                                  
#> [6] "mmrgaaptgvvsvmvlmilvllkqiesasangsslglpprkfcniyqgswvydksyplydskncpfierqfncksngrpdseylkyrwqpsgcnlprfngldflgrimkgkklmfvgdslslnqwqsltcllhnaapkanststrsosglsvfsfpaynssimfsrnaflvdivgappkrvmkldsissgslwktadvlvfnswhwwlhtdrkqpwdaimsgnvtvkdmdrlvayekammtwakwidqnidpsktkvffqgispdhgrarewskqggkgscigetkpimgssylagphaaemvvakviktmknqarlmdvtlmsqlrkdghpsvygfgghrmadcshwclsgvpdswnqllyselfhs"
#> [7] "maliknnifftsllifvtlfgvavggtvhkvgntkgwtmiggdyeawassrvfqvgdtlvfaynkdyhdvtevthndfemcesskplrryktgsdsisltkpglqhficgvpghckkgqklqihvlpaslghvavovogovrsqssssSOSOSOlvdpovnnapqyqmgotoashsaasadfiftfsfdltlidlctffilffilv"                                                                                                                                                                  
#> [8] "marsfaiavicivliagvtgqAOTSOOTaTOAOOTOTTOoOAaTOoovsAOoovttsooovttaoooanoooovssoooasoooatooovasoooovasoooaTOoovaTOoOAOlaSOOAqvOAOAOTtkodSOSOSOSsSOolOSsdAOgOStdsiSOAOSOTdvndqngaskmvsslvfgsvlvwfmi"
```

to switch of all PPP+ from being detected:

``` r
at_nsp_ag <- scan_ag(hyp_pred$sequence,
                     sequence, 
                     id,
                     dim = 3, #at least 3 dimers must be present
                     div = 10, #no more than 10 amino acids apart
                     type = "conservative",
                     exclude_ext = "all")
#> sequence vector contains O, O will be considered instead of P

at_nsp_ag$sequence[ind]
#> [1] "mayatilmifsvvalmsgerahaavdcsslilnmadclsfvtsgstvvkpegtccsglktvvrtgpeclceafknsgslgltldlskaaslpsvckvaappsarcglsvsgdpOAtAOglSOTagagAOAlssganaatpvssprssdasllsvsfafvifmalissfy"                                                                                                                                                                                                       
#> [2] "mgyrrsyaitfvalvaalwsvtkaqpssscvstlttlspclsyitgnsttpsqpccsrldsviksspqcicsavnspipniglninrtqalqlpnacniqtppltqcnaatgOTaqpOAOSOTekTOdvtlTOTslogarsgvgggsktvpsvgtgsssrnvdplplhflmfavlvvctssfl"                                                                                                                                                                                         
#> [3] "mrlllsllfllalttyssatyclcrdgvgekdlqtsidyacgvlkdcnpihekgpcyqpntikshcdwavntyfqrfgqisgscnfsgtattsqnlpstvvtgclypsSOgsagtTOTtgTOSgtqtfpgoOAfgOAgdfdpsgnngapslfisialslgfsvviafl"                                                                                                                                                                                                        
#> [4] "mhhhhhpcnrkpfttifsffllylnlhnqqiiearnpsqfttnpsodvsipeikrhlqqygylpqnkesddvsfeqalvryqknlglpitgkpdsdtlsqillprcgfpddvepktapfhtgkkyvyfpgrprwtrdvplkltyafsqenltpylaptdirrvfrrafgkwasvipvsfietedyviadikigffngdhgdgepfdgvlgvlahtfspengrlhldkaetwavdfdeekssvavdlesvavheighvlglghssvkdaamyptlkprskkvnlnmddvvgvqslygtnpnftlnsllasetstnladgsrirsqgmiystlstvialcflnw"                          
#> [5] "makmqlsifiavvalivcsasaktaSOOAOvloOTOAOAOAOenvnltellsvagpfhtfldyllstgvietfqnqannteegitifvpkddafkaqknpplsnltkdqlkqlvlfhalphyyslsefknlsqsgpvstfaggqyslkftdvsgtvridslwtrtkvsssvfstdpvavyqvnrvllpeaifgtdvoompAOAOAOivsAOSdSOSvadsegasspksshknsgqklllapismvisglvalfl"                                                                                                                  
#> [6] "mmrgaaptgvvsvmvlmilvllkqiesasangsslglpprkfcniyqgswvydksyplydskncpfierqfncksngrpdseylkyrwqpsgcnlprfngldflgrimkgkklmfvgdslslnqwqsltcllhnaapkanststrsosglsvfsfpaynssimfsrnaflvdivgappkrvmkldsissgslwktadvlvfnswhwwlhtdrkqpwdaimsgnvtvkdmdrlvayekammtwakwidqnidpsktkvffqgispdhgrarewskqggkgscigetkpimgssylagphaaemvvakviktmknqarlmdvtlmsqlrkdghpsvygfgghrmadcshwclsgvpdswnqllyselfhs"
#> [7] "maliknnifftsllifvtlfgvavggtvhkvgntkgwtmiggdyeawassrvfqvgdtlvfaynkdyhdvtevthndfemcesskplrryktgsdsisltkpglqhficgvpghckkgqklqihvlpaslghvavovogovrsqssssSOSOSOlvdpovnnapqyqmgotoashsaasadfiftfsfdltlidlctffilffilv"                                                                                                                                                                  
#> [8] "marsfaiavicivliagvtgqAOTSOOTaTOAOOTOTtoooaatooovsaooovttsooovttaoooanoooovssoooasoooatooovasoooovasoooatooovatoooAOlaSOOAqvOAOAOTtkodSOSOSOSsSOolOSsdAOgOStdsiSOAOSOTdvndqngaskmvsslvfgsvlvwfmi"
```

scan motifs for N-glycosylation
-------------------------------

N-glycosylation is frequent in secreted proteins. To scan the amino acid motifs for N-glycosylation the function `scan_nglc` can be used. Detection is based on [PROSITE pattern PS00001](https://prosite.expasy.org/PDOC00001). Mean local hydrophilicity (Hopp and Woods 1981) is used to assess if the asparagines are buried:

``` r
at_nglc <- scan_nglc(at_nsp,
                     sequence, 
                     Transcript.id)
```

create a protein structure diagram
----------------------------------

Using the above mentioned functions several protein features can obtained. To visualize them the function `plot_prot` can be used. Domains (`get_hmm`), N-terminal signal peptides (`get_signalp`), transmembrane regions, extracellular and intracellular protein regions (`get_phobius`), GPI attachment sites (`get_big_pi`), AG glycomodul spans (`scan_ag`) and hydroxyproline positions (`predict_hyp`) can be shown (all on by default). The output is a `ggplot2` object which can be manipulated to customize the theme, colors and plot annotations.

``` r
ind <- c(23, 24, 5, 80, 81, 120, 230, 334, 345, 1000)
pred <- plot_prot(sequence = at_nsp$sequence[ind],
                  id = at_nsp$Transcript.id[ind])

pred
```

![alt text](/man/figures/README-plot_prot-1.png)

References:
-----------

Eisenhaber, Birgit, Michael Wildpaner, Carolyn J. Schultz, Georg H. H. Borner, Paul Dupree, and Frank Eisenhaber. 2003. “Glycosylphosphatidylinositol Lipid Anchoring of Plant Proteins. Sensitive Prediction from Sequence- and Genome-Wide Studies for Arabidopsis and Rice.” *Plant Physiology* 133 (4): 1691–1701. doi:[10.1104/pp.103.023580](https://doi.org/10.1104/pp.103.023580).

Ellis, Miriam, Jack Egelund, Carolyn J. Schultz, and Antony Bacic. 2010. “Arabinogalactan-Proteins: Key Regulators at the Cell Surface?” *Plant Physiology* 153 (2): 403–19. doi:[10.1104/pp.110.156000](https://doi.org/10.1104/pp.110.156000).

Emanuelsson, O., H. Nielsen, S. Brunak, and G. von Heijne. 2000. “Predicting Subcellular Localization of Proteins Based on Their N-Terminal Amino Acid Sequence.” *Journal of Molecular Biology* 300 (4): 1005–16. doi:[10.1006/jmbi.2000.3903](https://doi.org/10.1006/jmbi.2000.3903).

Finn, Robert D., Jody Clements, and Sean R. Eddy. 2011. “HMMER Web Server: Interactive Sequence Similarity Searching.” *Nucleic Acids Research* 39 (Web Server issue): W29–W37. doi:[10.1093/nar/gkr367](https://doi.org/10.1093/nar/gkr367).

Hopp, T. P., and K. R. Woods. 1981. “Prediction of Protein Antigenic Determinants from Amino Acid Sequences.” *Proceedings of the National Academy of Sciences of the United States of America* 78 (6): 3824–8.

Johnson, Kim L., Andrew M. Cassin, Andrew Lonsdale, Antony Bacic, Monika S. Doblin, and Carolyn J. Schultz. 2017. “Pipeline to Identify Hydroxyproline-Rich Glycoproteins.” *Plant Physiology* 174 (2): 886–903. doi:[10.1104/pp.17.00294](https://doi.org/10.1104/pp.17.00294).

Käll, Lukas, Anders Krogh, and Erik L. L. Sonnhammer. 2007. “Advantages of Combined Transmembrane Topology and Signal Peptide Prediction–the Phobius Web Server.” *Nucleic Acids Research* 35 (Web Server issue): W429–432. doi:[10.1093/nar/gkm256](https://doi.org/10.1093/nar/gkm256).

Pierleoni, Andrea, Pier Luigi Martelli, and Rita Casadio. 2008. “PredGPI: A GPI-Anchor Predictor.” *BMC Bioinformatics* 9 (September): 392. doi:[10.1186/1471-2105-9-392](https://doi.org/10.1186/1471-2105-9-392).

Walsh, Ian, Alberto J. M. Martin, Tomàs Di Domenico, and Silvio C. E. Tosatto. 2012. “ESpritz: Accurate and Fast Prediction of Protein Disorder.” *Bioinformatics* 28 (4): 503–9. doi:[10.1093/bioinformatics/btr682](https://doi.org/10.1093/bioinformatics/btr682).
