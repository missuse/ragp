
<!-- README.md is generated from README.Rmd. Please edit that file -->
ragp
====

The goal of ragp is to facilitate mining of plant hydroxyproline rich glycoproteins (HGRPs), and especially arabinogalactan protein sequences (AGPs) from NGS data. The functions in the package consist of two groups: 1. web server prediction scraping and 2. sequence analyses. Functions in group 2 were made with HGRPs in mind and will be of limited use out of this scope while functions in group 1 are available for a wide range of annotation applications such as domain prediction, cellular localization, presence of trans-membrane regions and GO annotation. All of ragp functions expect protein sequences in single letter code.

Installation
------------

You can install ragp from github with:

``` r
# install.packages("devtools")
devtools::install_github("missuse/ragp")
```

Example
-------

This is a basic example which shows you how to fetch [SignalP](http://www.cbs.dtu.dk/services/SignalP/) predictions:

``` r
library(seqinr) #to create a fasta file with protein sequences
library(ragp)
data(at_nsp) #a data frame of 2700 Arabidopsis protein sequences 

#produce a fasta file from at_nsp data abailible in ragp:
seqinr::write.fasta(sequence = strsplit(at_nsp$sequence, ""),
                    name = at_nsp$Transcript.id, file = "at_nsp.fasta")

file_list <- split_fasta(path_in = "at_nsp.fasta", #path to the FASTA formated file`
                         path_out = "splited_at_nsp", #path to the splited files on which integers and file type will be apended automatically
                         num_seq = 1000) #number of sequences in destination files, usually this will be 10 - 20k

#get the SignalP predictions
signalp_pred_1 <- get_signalp_file(file = file_list[1])
```

First at\_nsp.fasta file was generated from the protein sequences using the library seqinr. This file was then split into smaller fasta files each containing a 1000 sequences using the ragp function split\_fasta. This is for illustration purposes - generally big fasta files should be split to smaller ones containing 10 - 20 thousand sequences. And finally SignalP predictions were obtained using the function get\_signalp\_file. Similarly [Phobius](http://phobius.sbc.su.se/) and [TargetP](http://www.cbs.dtu.dk/services/TargetP/) can be accessed by the functions get\_phobius\_file and get\_targetp\_file.

To fetch the GPI modification site predictions from [big-PI](http://mendel.imp.ac.at/gpi/plant_server.html) the function get\_big\_pi can be used:

``` r
ind <- c(145, 147, 160, 170,
         189, 203, 214, 224) #some indexes
big_pi_pred <- get_big_pi(sequence = at_nsp$sequence[ind],
                          id = at_nsp$Transcript.id[ind],
                          verbose = FALSE)
```

Similarly, domains can be identified by [hmmscan](https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan) using get\_hmm function:

``` r
pfam_pred <- get_hmm(sequence = at_nsp$sequence[1:20], #a vector of protein sequences as strings
                     id = at_nsp$Transcript.id[1:20],
                     verbose = FALSE)
```

To enrich with GO terms:

``` r
pfam_pred_go <- pfam2go(pfam_pred)
```

first 10 rows of the GO result with selected columns:

| id          | name             | acc        | GO\_name                                       | GO\_acc      |
|:------------|:-----------------|:-----------|:-----------------------------------------------|:-------------|
| ATCG00660.1 | Ribosomal\_L20   | PF00453.17 | <GO:translation>                               | <GO:0006412> |
| ATCG00660.1 | Ribosomal\_L20   | PF00453.17 | <GO:intracellular>                             | <GO:0005622> |
| ATCG00660.1 | Ribosomal\_L20   | PF00453.17 | <GO:rRNA> binding                              | <GO:0019843> |
| ATCG00660.1 | Ribosomal\_L20   | PF00453.17 | <GO:structural> constituent of ribosome        | <GO:0003735> |
| ATCG00660.1 | Ribosomal\_L20   | PF00453.17 | <GO:ribosome>                                  | <GO:0005840> |
| AT2G43600.1 | Glyco\_hydro\_19 | PF00182.18 | <GO:chitinase> activity                        | <GO:0004568> |
| AT2G43600.1 | Glyco\_hydro\_19 | PF00182.18 | <GO:cell> wall macromolecule catabolic process | <GO:0016998> |
| AT2G43600.1 | Glyco\_hydro\_19 | PF00182.18 | <GO:chitin> catabolic process                  | <GO:0006032> |
| AT2G43600.1 | Glyco\_hydro\_19 | PF00182.18 | <GO:chitinase> activity                        | <GO:0004568> |
| AT2G43600.1 | Glyco\_hydro\_19 | PF00182.18 | <GO:cell> wall macromolecule catabolic process | <GO:0016998> |

Count specific amino acids or amino acid motifs:

``` r
PAST_bias <- calculate_bias(sequence = at_nsp$sequence, #a vector of protein sequences as strings
                            id = at_nsp$Transcript.id, #a vector of protein identifiers as strings
                            user = c("P", "A", "S", "T"), #amino acids that should be counted - this is the default so it can be omitted in this case
                            simplify = TRUE) #output as data frame instead of a list
```

Filter sequences with total PAST% over 45%:

``` r
PAST45_bias <- PAST_bias[PAST_bias$bias_percent > 45,]
```

Count sequence motifs:

``` r
SPPP_bias <- calculate_bias(sequence = at_nsp$sequence,
                            id = at_nsp$Transcript.id, 
                            user = "SPPP", #motif that should be counted
                            simplify = TRUE)
```

Predict hydroxyproline sites in sequences:

``` r
hyp_pred <- predict_hyp(sequence = at_nsp$sequence, id = at_nsp$Transcript.id)
```

Output is a list of two elements. The first element "prediction"" is a data frame. First 10 rows:

| id          | substr                | P\_pos | prob               | HYP |
|:------------|:----------------------|:-------|:-------------------|:----|
| AT2G43600.1 | FSQNCMDTSCPGLKECCSRWG | 31     | 0.0376294292509556 | No  |
| AT2G43600.1 | EYCGFFCFSGPCNIKGKSYGY | 58     | 0.0523086786270142 | No  |
| AT2G43600.1 | YGYDYNVDAGPRGKIETVITS | 76     | 0.110943540930748  | No  |
| AT2G43600.1 | ERYCSKSKKYPCEPGKNYYGR | 163    | 0.075996607542038  | No  |
| AT2G43600.1 | CSKSKKYPCEPGKNYYGRGLL | 166    | 0.0564411692321301 | No  |
| AT2G43600.1 | YYGAGKHLGLPLLKDPDLVSR | 194    | 0.0403638146817684 | No  |
| AT2G43600.1 | KHLGLPLLKDPDLVSRSPEVA | 199    | 0.0384472720324993 | No  |
| AT2G43600.1 | LKDPDLVSRSPEVAFKFAMWF | 206    | 0.0404469892382622 | No  |
| AT2G43600.1 | AMWFWNRNVRPALYLGFGEIT | 223    | 0.0404815301299095 | No  |
| AT2G28410.1 | TNFALAQDRAPHGLAYETPVA | 27     | 0.0615824684500694 | No  |

predict\_hyp is also availible as a [shiny app](https://profenicolalac.shinyapps.io/HYPpredict_Shiny/). Details on how hydroxyproline sites are predicted will be availible soon.

AGPs are characterized by the presence of so called AGII glycomodules - amino acid dimers: OA, OS, OT, AO, SO and TO (and probably OG, OV, GO and VO) which are in close proximity to each other. Where: O - hydroxyproline, A - alanine, S - serine, T - threnonine, G - glycine and V - valine. scan\_ag function attempts to find the mentioned dimers according to user specified rules. Example:

``` r
at_nsp_ag <- scan_ag(sequence = at_nsp$sequence,
                     id = at_nsp$Transcript.id,
                     dim = 3, #at least 3 dimers must be present
                     div = 10, #no more than 10 amino acids apart
                     type = "conservative") #dimers will be defined as: PA, PS, PT, AP, SP, TP
```

List several sequences. possible AGII glycomodules are in caps:

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

Before I wrote OA, OS, OT, AO, SO and TO (and probably OG, OV, GO and VO) are AGII glycomodules but the above output considers P instead of O since most of the time the positions of O are unknown. However if the sequence argument to scan\_ag contains O's at some positions, scan\_ag will consider only them. To do this the sequence output of predict\_hyp function can be used. Example:

``` r
at_nsp_ag <- scan_ag(sequence = hyp_pred$sequence, #hyp_pred was created a little bit back
                     id = at_nsp$Transcript.id,
                     dim = 3, #at least 3 dimers must be present
                     div = 10, #no more than 10 amino acids apart
                     type = "conservative")
#> sequence vector contains O, O will be considered instead of P

at_nsp_ag$sequence[ind]
#> [1] "mayatilmifsvvalmsgerahaavdcsslilnmadclsfvtsgstvvkpegtccsglktvvrtgpeclceafknsgslgltldlskaaslpsvckvaappsarcglsvsgdoOAtAOglSOTagagAOAlssganaaTOvsSOrssdasllsvsfafvifmalissfy"                                                                                                                                                                                                       
#> [2] "mgyrrsyaitfvalvaalwsvtkaqpssscvstlttlspclsyitgnsttpsqpccsrldsviksspqcicsavnspipniglninrtqalqlpnacniqtppltqcnaatgOTaqoOAOSOTekTOdvtlTOTslogarsgvgggsktvosvgtgsssrnvdplplhflmfavlvvctssfl"                                                                                                                                                                                         
#> [3] "mrlllsllfllalttyssatyclcrdgvgekdlqtsidyacgvlkdcnpihekgpcyqpntikshcdwavntyfqrfgqisgscnfsgtattsqnlpstvvtgclypsSOgsagtTOTtgTOSgtqtfogoOAfgOAgdfdOSgnngapslfisialslgfsvviafl"                                                                                                                                                                                                        
#> [4] "mhhhhhpcnrkpfttifsffllylnlhnqqiiearnpsqfttnpsodvsipeikrhlqqygylpqnkesddvsfeqalvryqknlglpitgkpdsdtlsqillprcgfpddvepktapfhtgkkyvyfpgrprwtrdvplkltyafsqenltpylaptdirrvfrrafgkwasvipvsfietedyviadikigffngdhgdgepfdgvlgvlahtfspengrlhldkaetwavdfdeekssvavdlesvavheighvlglghssvkdaamyptlkprskkvnlnmddvvgvqslygtnpnftlnsllasetstnladgsrirsqgmiystlstvialcflnw"                          
#> [5] "makmqlsifiavvalivcsasaktaSOOAOvloOTOAOAOAOenvnltellsvagpfhtfldyllstgvietfqnqannteegitifvpkddafkaqknpplsnltkdqlkqlvlfhalphyyslsefknlsqsgpvstfaggqyslkftdvsgtvridslwtrtkvsssvfstdpvavyqvnrvllpeaifgtdvpompAOAOAOivsAOSdSOSvadsegasspksshknsgqklllapismvisglvalfl"                                                                                                                  
#> [6] "mmrgaaptgvvsvmvlmilvllkqiesasangsslglpprkfcniyqgswvydksyplydskncpfierqfncksngrpdseylkyrwqpsgcnlprfngldflgrimkgkklmfvgdslslnqwqsltcllhnaapkanststrsosglsvfsfpaynssimfsrnaflvdivgappkrvmkldsissgslwktadvlvfnswhwwlhtdrkqpwdaimsgnvtvkdmdrlvayekammtwakwidqnidpsktkvffqgispdhgrarewskqggkgscigetkpimgssylagphaaemvvakviktmknqarlmdvtlmsqlrkdghpsvygfgghrmadcshwclsgvpdswnqllyselfhs"
#> [7] "maliknnifftsllifvtlfgvavggtvhkvgntkgwtmiggdyeawassrvfqvgdtlvfaynkdyhdvtevthndfemcesskplrryktgsdsisltkpglqhficgvpghckkgqklqihvlpaslghvavovogovrsqssssSOSOSOlvdpovnnapqyqmgotoashsaasadfiftfsfdltlidlctffilffilv"                                                                                                                                                                  
#> [8] "marsfaiavicivliagvtgqAOTSOOTaTOAOOTOTTOoOAaTOoovsAOoovttSOoovttAOoOAnoooovsSOoOASOoOATOoovaSOooovaSOoOATOoovaTOoOAOlaSOOAqvOAOAOTtkodSOSOSOSsSOolOSsdAOgOStdsiSOAOSOTdvndqngaskmvsslvfgsvlvwfmi"
```

Extensin motifs in the form of SPPP+ are also detected by scan\_ag if in correct context, to avoid this add exclude\_ext = "yes" as an argument:

``` r
at_nsp_ag <- scan_ag(sequence = hyp_pred$sequence, 
                     id = at_nsp$Transcript.id,
                     dim = 3, #at least 3 dimers must be present
                     div = 10, #no more than 10 amino acids apart
                     type = "conservative",
                     exclude_ext = "yes")
#> sequence vector contains O, O will be considered instead of P

at_nsp_ag$sequence[ind]
#> [1] "mayatilmifsvvalmsgerahaavdcsslilnmadclsfvtsgstvvkpegtccsglktvvrtgpeclceafknsgslgltldlskaaslpsvckvaappsarcglsvsgdoOAtAOglSOTagagAOAlssganaaTOvsSOrssdasllsvsfafvifmalissfy"                                                                                                                                                                                                       
#> [2] "mgyrrsyaitfvalvaalwsvtkaqpssscvstlttlspclsyitgnsttpsqpccsrldsviksspqcicsavnspipniglninrtqalqlpnacniqtppltqcnaatgOTaqoOAOSOTekTOdvtlTOTslogarsgvgggsktvosvgtgsssrnvdplplhflmfavlvvctssfl"                                                                                                                                                                                         
#> [3] "mrlllsllfllalttyssatyclcrdgvgekdlqtsidyacgvlkdcnpihekgpcyqpntikshcdwavntyfqrfgqisgscnfsgtattsqnlpstvvtgclypsSOgsagtTOTtgTOSgtqtfogoOAfgOAgdfdOSgnngapslfisialslgfsvviafl"                                                                                                                                                                                                        
#> [4] "mhhhhhpcnrkpfttifsffllylnlhnqqiiearnpsqfttnpsodvsipeikrhlqqygylpqnkesddvsfeqalvryqknlglpitgkpdsdtlsqillprcgfpddvepktapfhtgkkyvyfpgrprwtrdvplkltyafsqenltpylaptdirrvfrrafgkwasvipvsfietedyviadikigffngdhgdgepfdgvlgvlahtfspengrlhldkaetwavdfdeekssvavdlesvavheighvlglghssvkdaamyptlkprskkvnlnmddvvgvqslygtnpnftlnsllasetstnladgsrirsqgmiystlstvialcflnw"                          
#> [5] "makmqlsifiavvalivcsasaktaSOOAOvloOTOAOAOAOenvnltellsvagpfhtfldyllstgvietfqnqannteegitifvpkddafkaqknpplsnltkdqlkqlvlfhalphyyslsefknlsqsgpvstfaggqyslkftdvsgtvridslwtrtkvsssvfstdpvavyqvnrvllpeaifgtdvpompAOAOAOivsAOSdSOSvadsegasspksshknsgqklllapismvisglvalfl"                                                                                                                  
#> [6] "mmrgaaptgvvsvmvlmilvllkqiesasangsslglpprkfcniyqgswvydksyplydskncpfierqfncksngrpdseylkyrwqpsgcnlprfngldflgrimkgkklmfvgdslslnqwqsltcllhnaapkanststrsosglsvfsfpaynssimfsrnaflvdivgappkrvmkldsissgslwktadvlvfnswhwwlhtdrkqpwdaimsgnvtvkdmdrlvayekammtwakwidqnidpsktkvffqgispdhgrarewskqggkgscigetkpimgssylagphaaemvvakviktmknqarlmdvtlmsqlrkdghpsvygfgghrmadcshwclsgvpdswnqllyselfhs"
#> [7] "maliknnifftsllifvtlfgvavggtvhkvgntkgwtmiggdyeawassrvfqvgdtlvfaynkdyhdvtevthndfemcesskplrryktgsdsisltkpglqhficgvpghckkgqklqihvlpaslghvavovogovrsqssssSOSOSOlvdpovnnapqyqmgotoashsaasadfiftfsfdltlidlctffilffilv"                                                                                                                                                                  
#> [8] "marsfaiavicivliagvtgqAOTSOOTaTOAOOTOTTOoOAaTOoovsAOoovttsooovttaoooanoooovssoooasoooatooovasoooovasoooaTOoovaTOoOAOlaSOOAqvOAOAOTtkodSOSOSOSsSOolOSsdAOgOStdsiSOAOSOTdvndqngaskmvsslvfgsvlvwfmi"
```

to switch of all PPP+ from beeing detected:

``` r
at_nsp_ag <- scan_ag(sequence = hyp_pred$sequence,
                     id = at_nsp$Transcript.id,
                     dim = 3, #at least 3 dimers must be present
                     div = 10, #no more than 10 amino acids apart
                     type = "conservative",
                     exclude_ext = "all")
#> sequence vector contains O, O will be considered instead of P

at_nsp_ag$sequence[ind]
#> [1] "mayatilmifsvvalmsgerahaavdcsslilnmadclsfvtsgstvvkpegtccsglktvvrtgpeclceafknsgslgltldlskaaslpsvckvaappsarcglsvsgdoOAtAOglSOTagagAOAlssganaaTOvsSOrssdasllsvsfafvifmalissfy"                                                                                                                                                                                                       
#> [2] "mgyrrsyaitfvalvaalwsvtkaqpssscvstlttlspclsyitgnsttpsqpccsrldsviksspqcicsavnspipniglninrtqalqlpnacniqtppltqcnaatgOTaqoOAOSOTekTOdvtlTOTslogarsgvgggsktvosvgtgsssrnvdplplhflmfavlvvctssfl"                                                                                                                                                                                         
#> [3] "mrlllsllfllalttyssatyclcrdgvgekdlqtsidyacgvlkdcnpihekgpcyqpntikshcdwavntyfqrfgqisgscnfsgtattsqnlpstvvtgclypsSOgsagtTOTtgTOSgtqtfogoOAfgOAgdfdOSgnngapslfisialslgfsvviafl"                                                                                                                                                                                                        
#> [4] "mhhhhhpcnrkpfttifsffllylnlhnqqiiearnpsqfttnpsodvsipeikrhlqqygylpqnkesddvsfeqalvryqknlglpitgkpdsdtlsqillprcgfpddvepktapfhtgkkyvyfpgrprwtrdvplkltyafsqenltpylaptdirrvfrrafgkwasvipvsfietedyviadikigffngdhgdgepfdgvlgvlahtfspengrlhldkaetwavdfdeekssvavdlesvavheighvlglghssvkdaamyptlkprskkvnlnmddvvgvqslygtnpnftlnsllasetstnladgsrirsqgmiystlstvialcflnw"                          
#> [5] "makmqlsifiavvalivcsasaktaSOOAOvloOTOAOAOAOenvnltellsvagpfhtfldyllstgvietfqnqannteegitifvpkddafkaqknpplsnltkdqlkqlvlfhalphyyslsefknlsqsgpvstfaggqyslkftdvsgtvridslwtrtkvsssvfstdpvavyqvnrvllpeaifgtdvpompAOAOAOivsAOSdSOSvadsegasspksshknsgqklllapismvisglvalfl"                                                                                                                  
#> [6] "mmrgaaptgvvsvmvlmilvllkqiesasangsslglpprkfcniyqgswvydksyplydskncpfierqfncksngrpdseylkyrwqpsgcnlprfngldflgrimkgkklmfvgdslslnqwqsltcllhnaapkanststrsosglsvfsfpaynssimfsrnaflvdivgappkrvmkldsissgslwktadvlvfnswhwwlhtdrkqpwdaimsgnvtvkdmdrlvayekammtwakwidqnidpsktkvffqgispdhgrarewskqggkgscigetkpimgssylagphaaemvvakviktmknqarlmdvtlmsqlrkdghpsvygfgghrmadcshwclsgvpdswnqllyselfhs"
#> [7] "maliknnifftsllifvtlfgvavggtvhkvgntkgwtmiggdyeawassrvfqvgdtlvfaynkdyhdvtevthndfemcesskplrryktgsdsisltkpglqhficgvpghckkgqklqihvlpaslghvavovogovrsqssssSOSOSOlvdpovnnapqyqmgotoashsaasadfiftfsfdltlidlctffilffilv"                                                                                                                                                                  
#> [8] "marsfaiavicivliagvtgqAOTSOOTaTOAOOTOTtoooaatooovsaooovttsooovttaoooanoooovssoooasoooatooovasoooovasoooatooovatoooAOlaSOOAqvOAOAOTtkodSOSOSOSsSOolOSsdAOgOStdsiSOAOSOTdvndqngaskmvsslvfgsvlvwfmi"
```
