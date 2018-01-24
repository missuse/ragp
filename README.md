
<!-- README.md is generated from README.Rmd. Please edit that file -->
ragp
====

The goal of ragp is to facilitate mining of plant hydroxyproline rich glycoproteins (HGRPs), and especially arabinogalactan protein sequences (AGPs) from NGS data. The functions in the package consist of two groups: 1. web server prediction scraping and 2. sequence analyses. Functions in group 2 were made with HGRPs in mind and will be of limited use out of this scope while functions in group 1 are available for a wide range of annotation applications such as domain prediction, cellular localization, presence of trans-membrane regions and GO annotation.

Installation
------------

You can install ragp from github with:

``` r
# install.packages("devtools")
devtools::install_github("missuse/ragp")
```

Examples
--------

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

First at\_nsp.fasta file was generated from the protein sequences using the library seqinr. Then this file was split to fasta files each containing a 1000 sequences using the ragp function split\_fasta. This is for illustration purposes, generally big fasta files should be split to smaller containing 10 - 20 thousand sequences. And finally SignalP predictions were obtained using the function get\_signalp\_file.

Similary, domains can be identified by [hmmscan](https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan) using get\_hmm function:

``` r
pfam_pred <- get_hmm(sequence = at_nsp$sequence[1:20], #a vector of protein sequences as strings
                     id = at_nsp$Transcript.id[1:20],
                     verbose = F)
```

To enrich with GO terms:

``` r
pfam_pred_go <- pfam2go(pfam_pred)
```

first 10 rows of the GO result with selected columns:

| id          | name             | Pfam\_name       | GO\_name                                       | GO\_acc      |
|:------------|:-----------------|:-----------------|:-----------------------------------------------|:-------------|
| ATCG00660.1 | Ribosomal\_L20   | Ribosomal\_L20   | <GO:translation>                               | <GO:0006412> |
| ATCG00660.1 | Ribosomal\_L20   | Ribosomal\_L20   | <GO:intracellular>                             | <GO:0005622> |
| ATCG00660.1 | Ribosomal\_L20   | Ribosomal\_L20   | <GO:rRNA> binding                              | <GO:0019843> |
| ATCG00660.1 | Ribosomal\_L20   | Ribosomal\_L20   | <GO:structural> constituent of ribosome        | <GO:0003735> |
| ATCG00660.1 | Ribosomal\_L20   | Ribosomal\_L20   | <GO:ribosome>                                  | <GO:0005840> |
| AT2G43600.1 | Glyco\_hydro\_19 | Glyco\_hydro\_19 | <GO:chitinase> activity                        | <GO:0004568> |
| AT2G43600.1 | Glyco\_hydro\_19 | Glyco\_hydro\_19 | <GO:cell> wall macromolecule catabolic process | <GO:0016998> |
| AT2G43600.1 | Glyco\_hydro\_19 | Glyco\_hydro\_19 | <GO:chitin> catabolic process                  | <GO:0006032> |
| AT2G43600.1 | Glyco\_hydro\_19 | Glyco\_hydro\_19 | <GO:chitinase> activity                        | <GO:0004568> |
| AT2G43600.1 | Glyco\_hydro\_19 | Glyco\_hydro\_19 | <GO:cell> wall macromolecule catabolic process | <GO:0016998> |

Count specific amino acids or amino acid motifs:

``` r
PAST_bias <- calculate_bias(sequence = at_nsp$sequence, #a vector of protein sequences as strings
                            id = at_nsp$Transcript.id, #a vector of protein identifiers as strings
                            user = c("P", "A", "S", "T"), #amino acids that should be counted
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
ind <- c(129, 145, 147, 160, 170,
         180, 189, 203, 205, 214, 217, 224)

hyp_pred <- predict_hyp(sequence = at_nsp$sequence[ind], id = at_nsp$Transcript.id[ind])
```

Output is a data frame. First 10 rows of the prediction:

| id          | substr                |  P\_pos|       prob| HYP |
|:------------|:----------------------|-------:|----------:|:----|
| AT2G20700.1 | ISPYCLLSLLPIFLLSGFSLS |      13|  0.0694811| No  |
| AT2G20700.1 | YTIITSRCKGPNYPANVCCSA |      63|  0.0469144| No  |
| AT2G20700.1 | ITSRCKGPNYPANVCCSAFKD |      66|  0.0416294| No  |
| AT2G20700.1 | CCSAFKDFACPFAEVLNDEKN |      80|  0.0376600| No  |
| AT2G20700.1 | FSYINLYGRYPPGIFANMCKE |     107|  0.0599186| No  |
| AT2G20700.1 | SYINLYGRYPPGIFANMCKEG |     108|  0.0571255| No  |
| AT2G20700.1 | QSASATSDSIPRASTTASLAV |     139|  0.0643806| No  |
| AT2G13820.1 | FVTSGSTVVKPEGTCCSGLKT |      50|  0.0379427| No  |
| AT2G13820.1 | SGLKTVVRTGPECLCEAFKNS |      66|  0.0412103| No  |
| AT2G13820.1 | TLDLSKAASLPSVCKVAAPPS |      92|  0.0615301| No  |

Chimeric AGPs contain additional domains apart the arabinogalactan sequence spans, and are usually much harder to find. They are characterized by the presence of so called AGII glycomodules - amino acid dimers: OA, OS, OT, AO, SO and TO (and probably OG, OV, GO and VO) which are in close proximity to each other. Where: O - hydroxyproline, A - alanine, S - serine, T - threnonine, G - glycine and V - valine. scan\_ag function attempts to find the mentioned dimers according to user specified rules. Example:

``` r
at_nsp_ag <- scan_ag(sequence = at_nsp$sequence,
                     id = at_nsp$Transcript.id,
                     dim = 3, #at least 3 dimers must be present
                     div = 10, #no more than 10 amino acids apart
                     type = "conservative") #dimers will be defined as: PA, PS, PT, AP, SP, TP
```

Output is a list. Which sequences contain more than 10 amino acids in dimers:

``` r
index_10 <- which(at_nsp_ag$AG_aa >= 10)
```

List several sequences. Dimers are in caps:

<table>
<colgroup>
<col width="100%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">sequence</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">mlsstsptlssvssfykipissndfsplslslslslfllsspssfslsplpeiftqaffcgsevkkmnhcnlqqnafmsreemmgfdrkdlvvcpkprrvgllannvirplrlhmsqaaadlcdskagaelleiirrkedngtigqllssSPpyfpgSPPSraanplaqdarfrdeklnpiSPnSPflqpysatgfPSPSssssssssrgcvrmkfglnspavrvegfdclnrdrqnssipama</td>
</tr>
<tr class="even">
<td align="left">mkslcfistAPlirSPpfldlSPTslnrftlkiPAslyvshPTrcsisnPSsseellnsnggmsrasisvfggtslnnlkmqvgspislhsinplaklslsdqaflllafivcttsvaftslvitaiptlvamgraatsfakladtarkelpstlaavrlsgmeisdltlelsdlsqditdginksakavqaaeagikqigtlaqqqtlsmieeranlpeislqpvvagaaektshaigsatkrlmniitggnkded</td>
</tr>
<tr class="odd">
<td align="left">marsfaiavicivliagvtgqAPTSPPTaTPAPPTPTTPpPAaTPppvsAPppvttSPppvttAPpPAnppppvsSPpPASPpPATPppvaSPpppvaSPpPATPppvaTPpPAPlaSPPAqvPAPAPTtkpdSPSPSPSsSPplPSsdAPgPStdsiSPAPSPTdvndqvsnlff</td>
</tr>
<tr class="even">
<td align="left">msvslftaftvlslclhtstsefqlstisaAPSflpeAPSsfsasTPAmSPdtSPlfPTPgssemSPSPSessimPTiPSslSPpnpdavTPdpllevSPvgSPlPAsssvclvssqlsslllvllmlllafcsff</td>
</tr>
<tr class="odd">
<td align="left">mgqsqvslflllilvfiygvssttftivnqcsytvwpgllsgagtSPlPTtgfslnPTetrvipiPAawsgriwgrtlctqdattgrftcitgdcgsstvecsgsgaappatlaeftlngangldfydvslvdgynipmtivpqgggdaggvagnctttgcvaelngpcpaqlkvattgaegvacksaceafgTPeyccsgafgTPdtckPSeysqffknacpraysyayddgtstftcggadyvitfcpspnpsvksatkgvqpvavsyskaspnasptlsavfsigvlavaswvmqrvl</td>
</tr>
</tbody>
</table>
