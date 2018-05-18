####################.
# Install packages #-----------------
####################.

devtools::install_github("missuse/ragp")         # contains MAAB
devtools::install_github("shaunpwilkinson/kmer")

library("kmer")
library("ragp")
library("igraph")

file <- "C:\\Users\\T\\OneDrive\\5-Storage\\3. Postdoc Bacic\\Initial data\\FLA\\FLAs\\Initial_masked_FLAs.fa"
file <- "C:\\Users\\T\\OneDrive\\5-Storage\\3. Postdoc Bacic\\Initial data\\FLA\\FLAs\\Putative FLAs (STAP region only).fa"
file.full <- "C:\\Users\\T\\OneDrive\\5-Storage\\3. Postdoc Bacic\\Initial data\\FLA\\FLAs\\Putative FLAs (full seq).fa"

SEQ  <- insect::readFASTA(file,
                         residues = "AA")
SEQ  <- SEQ[stringi::stri_length(SEQ)>=10]

SEQ1 <- Biostrings::readAAStringSet(file)
SEQF <- Biostrings::readAAStringSet(file.full)

#########.
# MAAB+ #--------------------------
#########.

#> FLA annotation -------------
SEQF2 <- lapply(SEQF,as.character)

t  <- Sys.time()
pb <- txtProgressBar(min = 0, max = total, style = 3)
for(i in 1:(length(SEQF2)/10)){
  SEQF3 <- SEQF2[i*10+(1:10)]
  annot <- rbind(annot, ragp::get_hmm(sequence = SEQF3,
                                      id = names(SEQF3),
                                      verbose = FALSE,
                                      sleep = 0.1))
  print(i)
  setTxtProgressBar(pb, i)
  print(Sys.time()-t)
}

annotF <- annot[annot[,2]=="Fasciclin",]
annotF <- annotF[!is.na(annotF[,2]),]

# Number of Fla domains
FLAcount<-NULL
for(i in names(SEQF2)){
  FLAcount[i] <- sum(annotF[,1]==i)
}

#> FLA extraction  ---------------
SEQF4 <- lapply(SEQF2,strsplit,"")
SEQF4 <- lapply(SEQF4,"[[",1)
SEQF4 <- lapply(SEQF4,casefold,upper=FALSE)
# uppercase FLAs
for(i in 1:nrow(annotF)){
  SEQF4[[annotF[i,1]]] [annotF[i,6]:annotF[i,7]] <- casefold(SEQF4[[annotF[i,1]]] [annotF[i,6]:annotF[i,7]],upper=TRUE)
}
as.fasta(SEQF4             ,write = "Full length (annotated).fa")
as.fasta(SEQF4[FLAcount==0],write = "Full length (annotated) 0-FLA.fa")
as.fasta(SEQF4[FLAcount==1],write = "Full length (annotated) 1-FLA.fa")
as.fasta(SEQF4[FLAcount==2],write = "Full length (annotated) 2-FLA.fa")
as.fasta(SEQF4[FLAcount>=3],write = "Full length (annotated) 3+FLA.fa")

counter<-NULL
for(i in 1:length(FLAcount)){
  if(FLAcount[i]>0){
    counter <- append(counter,1:FLAcount[i])
  }
}
FLAmax <- rep(FLAcount,FLAcount)
FLAnames <- paste0(annotF[,1],".FLA.",counter,"/",FLAmax)
SEQF5 <- SEQF4[FLAcount>0]
FLAs  <- list()
for(i in 1:nrow(annotF)){
  start <-annotF[i,6] 
  end   <-annotF[i,7] 
  if(annotF[i,8]<5){
    start <- start - annotF[i,8]          # if starts <5aa short, also collect missing end residues
  }
  if(start<0){                            # prevent negaive start numbers
    start <- 0
  }
  if((128-annotF[i,9])<5){                # if stops <5aa short, also collect missing end residues
    end <- end + (128-annotF[i,9])
  }
  FLAs[[ FLAnames[i] ]] <- SEQF5[[ annotF[i,1] ]][start:end]
}
FLAs <- lapply(FLAs, function(x) x[!is.na(x)])
FLAs <- FLAs[lapply(FLAs, length)>50]
as.fasta(FLAs,write = "ALLFLAS.fa")

FLAsX  <- list()
for(i in 1:nrow(annotF)){
  start <-annotF[i,6] - annotF[i,8]       - 5    # also collect N-ter residues not picked up by model and 5 extra
  end   <-annotF[i,7] + (128-annotF[i,8]) + 5    # also collect C-ter residues not picked up by model and 5 extra
  if(start<0){                            # prevent negaive start numbers
    start <- 0
  }
  FLAsX[[ FLAnames[i] ]] <- SEQF5[[ annotF[i,1] ]][start:end]
  
  start2 <- which(diff(unlist(FLAsX[i])==casefold(unlist(FLAsX[i]),TRUE))==1)[2]
  if(!is.na(start2)){
    FLAsX[[ FLAnames[i] ]] <- FLAsX[[ FLAnames[i] ]] [1:start2]
    }
}
FLAsX <- lapply(FLAsX, function(x) x[!is.na(x)])
FLAsX <- FLAsX[lapply(FLAsX, length)>50]
as.fasta(FLAsX,write = "ALLFLAS plus overhangs.fa")

boxplot(t(head(annotF[order(annotF[,9]-annotF[,8]),8:9],n=2658)))

#> FLA masking  -----------
SEQF6 <- SEQF4
for(i in 1:nrow(annotF)){
  SEQF6[[annotF[i,1]]] [annotF[i,6]:annotF[i,7]] <- "X"
}
as.fasta(SEQF6, write = "Full length (masked).fa") # X-masked FLA domains

#> AGP capitals ------
SEQF7 <- ragp::scan_ag(sequence = SEQF2,
                       id = names(SEQF2),
                       dim = 3,div = 8)$sequence
names <- paste0(">",names(SEQF2))
ord1 <- 2*(1:length(names))-1
ord2 <- 2*(1:length(SEQF7))

cat(c(names,SEQF7)[order(c(ord1,ord2))],sep = "\n",
    file =  "Full length (AGP capitals).fa") # Capitalised AGP regions

##########################.
# Inter-proline distance #--------------------------
##########################.

hist(stringi::stri_length(SEQ),breaks = 1000,xlim = c(0,400))

interP <- lapply(stringi::stri_split(SEQ1,regex = "P"),
                 stringi::stri_length)

hist(interP[[2]],breaks = -1:max(unlist(interP)),plot=FALSE)$count

interP.hist <- lapply (interP,hist,breaks = -1:max(unlist(interP)),plot=FALSE)
interP.hist <- lapply(interP.hist, `[[`, 2)
interP.hist <- t(array(unlist(interP.hist),
                       dim = c(max(unlist(interP))+1,
                               length(SEQ1))))

colnames(interP.hist)<- 0:max(unlist(interP))

interP.hist2 <- cbind(interP.hist[,1:21],
                      rowSums(interP.hist[,22:max(unlist(interP))]))
colnames(interP.hist2)[22]<-">20"

interP.hist2.norm <- interP.hist2/rowSums(interP.hist2)
heatmap(interP.hist2.norm,Colv = NA)

##> k2 -----------
k2 <- kmer::kcount(SEQ,2)
k2 <- k2/rowSums(k2)
k2 <- k2[,order(colSums(k2),decreasing = TRUE)]
k2.sums<-colSums(k2)/sum(k2)
as.matrix(percent(k.sums[1:20],2))

k2p <- k2[,grep("P",colnames(k2))]

plot(100*(colSums(k2)/sum(k2))[1:100],type = "b", pch=16)
txt<-colnames(k2)
txt[21:400]<-""
text(100*colSums(k2)/sum(k2),txt, cex=0.7,pos=4)

##> k3 -------------
k3 <- kmer::kcount(SEQ,3)
k3r<- k3
k3 <- k3/rowSums(k3)
k3 <- k3[,order(colSums(k3),decreasing = TRUE)]
k3.sums<-colSums(k3)/sum(k3)
as.matrix(percent(k3.sums[1:20],2))

k3p <- k3[,grep("P",colnames(k3))]
k3p.mid <- k3[,grep(".P.",colnames(k3))]

plot(100*(colSums(k3)/sum(k3))[1:100],type = "b")
txt<-colnames(k3)
txt[21:400]<-""
text(100*colSums(k3)/sum(k3),txt, cex=0.7,pos=4)


# Heatmanp and dendrogram ---------
topcols<-1:20
k2heat <- heatmap(k3[1:300,topcols],keep.dendro = TRUE)
k2tree <- ape::write.tree(
            ape::as.phylo(
              as.hclust(
                k2heat$Rowv
          )))

k3clust<-(kmer::cluster(SEQ,3))
plot(k3clust)


# Distance network ----------------

library(igraph)

k2dist   <- as.matrix(dist(k2))
k2dist20 <- as.matrix(dist(k2[,1:20]))
k2dist40 <- as.matrix(dist(k2[,1:40]))
k3dist   <- as.matrix(dist(k3))
k3dist20 <- as.matrix(dist(k3[,1:20]))
k3dist40 <- as.matrix(dist(k3[,1:40]))

data <- k3dist20[1:300,1:300]
data[data>=quantile(data,0.1)] <- 0

g <- igraph::graph.adjacency(as.matrix(data), mode="undirected", diag=FALSE, weighted=TRUE)
g.sub <- igraph::decompose.graph(g)
out<- NULL
for (i in 1:length(lengths(g.sub))){
  out <- append(out,length(g.sub[[i]][1]))
}
plot(out, type="h")
g2 <- g.sub[[which(out==max(out))]]


# Plot
plot(g,
     edge.width   = igraph::E(g)$weight,
     vertex.size  = 4,
     vertex.color = "red",
     vertex.label = NA)



# Numericise MSA
name      <- "p450"
folder    <- "C:\\Users\\T\\OneDrive\\0-Sequences\\5-Collaborations\\Drosophila p450s (Rane Batterham)"

res.prop1 <- "C:\\Users\\T\\OneDrive\\0-Sequences\\2-PCA\\0-Raw data and scalers\\Amino acid properties.csv"
res.prop2 <- "C:\\Users\\T\\OneDrive\\0-Sequences\\2-PCA\\0-Raw data and scalers\\Amino acid properties nodisord.csv"

MSA <- "Supp datannotations_allCYPAA.fasta"

# Seqspace --------------------------
#> Analysis -------------

res.prop1 <- "C:\\Users\\T\\OneDrive\\0-Sequences\\2-PCA\\0-Raw data and scalers\\Amino acid properties.csv"
MSA       <- "C:\\Users\\T\\OneDrive\\5-Storage\\3. Postdoc Bacic\\Initial data\\FLA\\Seq alignments\\ALLFLAS plus overhangs aligned.TrimAl99.fa"

#SAPCA <- PCA_MSA (MSA = "Defensins_database_aln2.fa",
SAPCA <- PCA_MSA (MSA        = MSA,
                  res.prop   = res.prop1,
                  cys        = 0,
                  # bootstrap  = 10,
                  clusterPCs = 1:40, 
                  clusters   = 1:20,
                  model      = "VVV")

Flanames2       <- gsub(paste0("X",FLAnames),replacement = ".", pattern = "[/]")
counter2        <- counter
names(counter2) <- Flanames2
counter2        <- counter2[SAPCA$numerical.alignment$seq.names]

FLAmax2         <- FLAmax
names(FLAmax2)  <- Flanames2
FLAmax2         <- FLAmax2[SAPCA$numerical.alignment$seq.names]
  
#> Plot ------------------
colours<-palette(c("white",       #1 A
                   "blue",     #2 B
                   "red",        #3 N
                   "darkgrey"))      #4 T


plot_modelfit(SAPCA,legend = FALSE)
plot_3Dclusters(SAPCA,
                plotPCs    = c(1,2,3),
                col = FLAmax2[SAPCA$numerical.alignment$seq.names])

rgl::par3d()->pp
rgl::par3d(pp)


# FLAnnotate! ----------------

sequences <- SEQF[1:100]

gpis <- NULL
pbt    <- txtProgressBar(min = 0, max = length(sequences), style = 3)
pbw    <- winProgressBar(min = 0, max = length(sequences), title = "HMM progress")

for(i in 1:length(sequences)){
  seqsubset <- sequences[i]
  gpis <- rbind(gpis, ragp::get_big_pi(sequence = seqsubset,
                                       id = names(seqsubset),
                                       verbose = FALSE,
                                       sleep = 0))
  setTxtProgressBar(pbt, i)
  setWinProgressBar(pbw, i, title= paste("GPI scan progress:",
                                         round(i/length(sequences)*100, 0),
                                         "%      (",
                                         names(sequences[i]),
                                         ")"
  ))
}
close(pbw)
rownames(gpis)<-gpis$id
signal.full <- ragp::get_phobius_file(file.full)
rownames(signal.full)<-signal.full$Name

a2 <- head(annot,44)
s2 <- SEQF[unique(a2$id)]
g2 <- gpis[unique(a2$id),]
sig2 <- signal.full [unique(a2$id),]

plot_domains(sequences   = s2,
             annotations = a2,
             gpis        = g2,
             signal      = sig2)      
