# Generating microsatellite calls and annotations for STRATIFY for a different reference genome

## Step 1: Running Tandem Repeats Finder (TRF)
1. Install TRF 4.09.1 from <a href="https://github.com/Benson-Genomics-Lab/TRF/releases/tag/v4.09.1">here</a>, if you haven't done so already.
2. Load fasta files for reference genome from <a href="https://hgdownload.soe.ucsc.edu/downloads.html#human">USCS Genome Browser</a> or elsewhere and split by chromosome (if not already split), using the following command:
```bash
awk '/^>/{sub(" .*",""); filename=$0; sub("^>","",filename); s=filename ".fasta"} {print > s}' [FILENAME.fasta]
```
3. Remove unused fasta files (e.g. HLA*.fasta, chrUn*.fasta, *_random.fasta, *_alt.fasta, chrEBV.fasta, etc.) using the ```rm``` command. Note: these files will be different depending on the reference genome used, so make sure to list all of the files split by chromosome by running
```bash
ls chr*.fasta
```
to see all the files output in step 2, and make sure to discard any files that you do not wish to include in your analyses.
4. Use TRF to find microsatellites on each chromosome. Note: it is best to run these commands in parallel on a computing cloud to speed up the process.
```bash
for i in `ls chr*.fasta`
do
  trf409.linux64 $i 2 7 7 80 10 36 6 -l 6 -h
done
```
5. Consolidate separate TRF outputs for each chromosome into a single file for the entire genome.
```bash
chr=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
for k in "${chr[@]}"
do
  awk 'NR>15 {print "chr" "'"$k"'" " " $0}' chr$k.fasta.2.7.7.80.10.36.6.dat >> tempoutput
done
```
6. Ensure that chromosomes are sorted in numerical order as they go into the consolidated file.
```bash
sort -V tempoutput > [FILENAME].dat
rm tempoutput
```
## Step 2: Import data into R & save as RDS for ease of future loading in and out of R to add annotations
```R
data = read.table("[FILENAME].dat")
header = c("chrom","chrom.start","chrom.stop","period.size","copy.num","consensus.size","per.match","per.indel","score","A","C","G","T","entropy","consensus","sequence")
colnames(data) <- header
saveRDS(data,"[FILENAME].RDS")
```
## Step 3: Remove microsatellites that overlap other microsatellites by >90%
1. Load in packages required for subsequent chunks:
```R
#Note: BiocManager is required for installation of all other packages listed here
library(BiocManager)
library(GenomicRanges)
library(plyranges)
library(Biostrings)
library(rtracklayer)
library(preprocessCore)
library(BSgenome.Hsapiens.UCSC.hg19) #or whatever other reference genome you are using
```
2. Write and run custom function to remove microsatellites that overlap other microsatellites by >90%:
```R
RemoveTRFOverlaps <- function(data.grange){
#find overlaping regions with granges object and itself, and calculate percent overlap.
hits <- findOverlaps(data.grange, ignore.strand=TRUE)
overlaps <- pintersect(data.grange[queryHits(hits)], data.grange[subjectHits(hits)])
percentOverlaps <- width(overlaps) / width(data.grange[subjectHits(hits)])

#Keep only intersections with > 0.9 fraction overlap and remove intersections between an element and itself
hits <- hits[percentOverlaps > 0.9 & (queryHits(hits)!=subjectHits(hits))]

#Reduce list so there is only 1 row for each intersection, instead of 2. For example, if there are query and subject hits of 1,5 and 5,1, it should be reduced to just 1,5.
hits.df <- data.frame(query=queryHits(hits),subject=subjectHits(hits))
hits.df <- t(apply(hits.df,1,sort))
hits.df <- unique(hits.df)

##Filter step 1: For each hit pair, find the element with the LOWER score.
item1_scores <- score(data.grange)[hits.df[,1]]
item2_scores <- score(data.grange)[hits.df[,2]]
item1_scores_lower <- item1_scores < item2_scores
item2_scores_lower <- item1_scores > item2_scores

#Retrieve all the item indices for later removal from data.grange
data.grange_itemstoremove <- unique(c(hits.df[item1_scores_lower,1],hits.df[item2_scores_lower,2]))

#Remove comparisons processed in this step from hits.df so they are not used in the next filter step.  
hits.df <- hits.df[!(item1_scores_lower | item2_scores_lower),]

##Filter step 2: hits.df now contains hits where both elements have the same score. For each hit pair, find the element with the LARGER period size.
item1_periodsizes <- data.grange$period.size[hits.df[,1]]
item2_periodsizes <- data.grange$period.size[hits.df[,2]]
item1_periodsizes_larger <- item1_periodsizes > item2_periodsizes
item2_periodsizes_larger <- item1_periodsizes < item2_periodsizes

#Retrieve all the item indices for later removal from data.grange
data.grange_itemstoremove <- c(data.grange_itemstoremove,unique(c(hits.df[item1_periodsizes_larger,1],hits.df[item2_periodsizes_larger,2])))

#Remove comparisons processed in this step from hits.df so they are not used in the next filter step.  
hits.df <- hits.df[!(item1_periodsizes_larger | item2_periodsizes_larger),]

##Filter step 3: hits.df now contains hits where both elements have the same score and the same period size. For each hit pair, find the element that is SMALLER in size.
item1_sizes <- width(data.grange)[hits.df[,1]]
item2_sizes <- width(data.grange)[hits.df[,2]]
item1_sizes_smaller <- item1_sizes < item2_sizes
item2_sizes_smaller <- item1_sizes > item2_sizes

#Retrieve all the item indices for later removal from data.grange
data.grange_itemstoremove <- c(data.grange_itemstoremove,unique(c(hits.df[item1_sizes_smaller,1],hits.df[item2_sizes_smaller,2])))

#Remove comparisons processed in this step from hits.df so they are not used in the next filter step.  
hits.df <- hits.df[!(item1_sizes_smaller | item2_sizes_smaller),]

##Filter step 4: hits.df now contains hits where both elements have the same score, period size, and element size. For each hit pair, find the element that has the LOWER entropy. We want to remove these because these are more perfect microsatellites and we should err on the side of less perfect microsatellites for a given locus.
item1_entropies <- data.grange$entropy[hits.df[,1]]
item2_entropies <- data.grange$entropy[hits.df[,2]]
item1_entropies_smaller <- item1_entropies < item2_entropies
item2_entropies_smaller <- item1_entropies > item2_entropies

#Retrieve all the item indices for later removal from data.grange
data.grange_itemstoremove <- c(data.grange_itemstoremove,unique(c(hits.df[item1_entropies_smaller,1],hits.df[item2_entropies_smaller,2])))

#Remove comparisons processed in this step from hits.df so they are not used in the next filter step.  
hits.df <- hits.df[!(item1_entropies_smaller | item2_entropies_smaller),]

##Remove items from the dataset that were found by the above filters
data.grange <- data.grange[-data.grange_itemstoremove]

return(data.grange)
}

#format data as a granges (the necessary data format for the function above)
data <- makeGRangesFromDataFrame(data, keep.extra.columns = T)
#run custom function above
data <- RemoveTRFOverlaps(data)
```
3. Check for >90% overlap remaining after running the function in step 2:
```R
#initialize overlap.remaining dataframe
overlap.remaining <- data.frame(matrix(ncol = 0, nrow = 0))
#find any microsatellites that overlap with another microsatellite and calculate percent overlap
hits <- findOverlaps(data, ignore.strand=TRUE)
overlaps <- pintersect(data[queryHits(hits)], data[subjectHits(hits)])
percentOverlaps <- width(overlaps) / width(data[subjectHits(hits)])
#filter for microsatellites that overlap with another microsatellite by more than 90% and is not the query microsatellite itself
hits <- hits[percentOverlaps > 0.9 & (queryHits(hits)!=subjectHits(hits))]
overlap.remaining <- pintersect(data[queryHits(hits)],
                                data[subjectHits(hits)])
if (!(isEmpty(overlap.remaining))){
  print("Overlap Remaining")
} else {
  print("No Overlap Remaining")
}
```
## Step 4: Define & annotate motif families
1. Generate a master dictionary of all possible motifs for period sizes 1-6, with reduction to a motif family for all reverse complements and circular permutations of a given motif. Then, invert the dictionary for rapid indexing:
```R
#initialize 6 different "motif dictionaries" for motif lengths 1-6
motif.dictionaries <- c("dictionary.1","dictionary.2","dictionary.3","dictionary.4","dictionary.5","dictionary.6")
for (ii in 1:length(motif.dictionaries)){
  #for each motif dictionary, get the motif length of the dictionary and find all possible combinations of nucleotides, given the specified motif length. Put those options in a variable called "motif.library"
  nucleotides <- c("A", "C", "T", "G")
  motif.length <- ii
  result <- as.matrix(expand.grid(lapply(numeric(motif.length), 
                                         function(x) nucleotides)),ncol=motif.length)
  motif.library <- apply(result,1,function(x) paste0(x,collapse=""))
  
  #for each motif in the motif library, define all circular permutations and reverse complements of the motif (including reverse complements of circular permutations). put all circular permutations and reverse complements of the motif in a "motif.set" variable. If the selected motif isn't already somewhere in the "motif.dictionary" data table, add the motif set to the motif.dictionary and rename the row of the motif dictionary to the search motif (this check prevents duplicate entries, eg. A row for an c("A","T") motif set. and a row for a c("T","A") motif set).
  motif.dictionary <- NULL
  for (jj in 1:length(motif.library)){
    motif <- motif.library[jj]
    circular_permutations = NULL
    rev_comps = NULL
    permutation <- DNAString(motif)
    for (kk in 1:nchar(motif)){
      move = substr(permutation,nchar(permutation),nchar(permutation))
      dont.move = substr(permutation,1,nchar(permutation)-1)
      permutation <- DNAString(paste0(move,dont.move))
      rev.comp <- reverseComplement(permutation)
      permutation <- toString(permutation)
      rev.comp <- toString(rev.comp)
      circular_permutations <- c(circular_permutations,permutation)
      rev_comps <- c(rev_comps,rev.comp)
      }
    motif.set <- c(circular_permutations,rev_comps)
    if (!(motif %in% motif.dictionary)){
      motif.dictionary <- rbind(motif.dictionary,motif.set)
      row.num <- nrow(motif.dictionary)
      rownames(motif.dictionary)[row.num] <- motif
    }
  }
  #once you have gone through all of the motifs in the motif library and added them and their circular permutations and reverse complements to the motif dictionary, assign the completed motif dictionary its proper name based on the motif length, and proceed to make the next motif dictionary until all 6 dictionaries are completed.
  assign(motif.dictionaries[ii], motif.dictionary)
}

#Each row should have all of the circular permutations and reverse complements of the "motif family" of a given size, but other dictionaries will have  repetitions of some of the smaller motifs that are also considered equivalent motifs for our motif families (e.g. "A" = "AA" = "AAA" = "T" = "TT" = "TTT" etc.). To address this, the following code loops through each of the six dictionaries. For each dictionary, it goes row by row and assesses whether the appropriately-sized repetition of the base motif is present in any of the other dictionaries. If one or more repetitions of the base motif are found in any of the other dictionaries, then that row in the queried dictionary is appended to the row in the search dictionary (the one the "base motif" came from) and it is removed from its original dictionary. After the loop has finished going through all of the motifs in all of the lines of the first dictionary, it then assesses repetitions of the motifs in the remaining lines of the second dictionary (note: there will be fewer lines in each subsequent dictionary because lines with repetitions of motifs that occur in dictionaries of smaller motifs will be appended to those dictionaries and removed from the dictionaries of larger motifs) and so on, until it has found all of the repetitions of the base motif in all of the other dictionaries. Then, that line is added to the "master.dictionary" list and is named after the "base motif". If the line is passed through without any revisions, it is still added to the "master.dictionary" list and named after the "base motif".
master.dictionary <- list()
for (ii in 1:length(motif.dictionaries)){
  dictionary.a <- get(motif.dictionaries[ii])
  other.dictionaries <- motif.dictionaries[which(motif.dictionaries != motif.dictionaries[ii])]
  for (jj in 1:length(dictionary.a[,1])){
    revised.line = dictionary.a[jj,]
    for (kk in 1:length(dictionary.a[1,])){
      base.motif = as.character(dictionary.a[jj,kk])
      rep = strrep(base.motif,6)
      for (mm in 1:length(other.dictionaries)){
        dictionary.b <- get(other.dictionaries[mm])
        key = substr(rep,1,nchar(as.character(dictionary.b[1,1])))
        if (is.na(key)){
          key = key
        } else if (nchar(key) <= nchar(base.motif)){
          key = key
        } else if (nchar(key)%%nchar(base.motif) == 0){
          key = key
        } else {
          key = NA
        }
        if (is.na(key)){
          next
        } else if (any(dictionary.b == key)){
          revised.line <- append(revised.line,dictionary.b[rownames(dictionary.b)[
            which(dictionary.b == key, arr.ind=TRUE)[1,1]],])
          dictionary.b <- dictionary.b[-c(which(dictionary.b == key, arr.ind=TRUE)[1,1]),]
        } else {
          revised.line = revised.line
        }
        assign(other.dictionaries[mm],dictionary.b)
      }
      if (is.na(base.motif)){
          next
        } else {
          master.dictionary[[base.motif]] <- revised.line
          dictionary.a[jj,] <- NA
          assign(motif.dictionaries[ii],dictionary.a)
        }
    }
  }
}

#After the single "Master Dictionary" has been generated, it is then inverted for ease of indexing, as it is faster and easier to use this way
inverted.master.dictionary <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 0))
for(i in 1:length(names(master.dictionary))){
  motif.vector <- unique(master.dictionary[[i]])
  motif.family <- rep(names(master.dictionary)[i],1,
                      length(motif.vector))
  temp <- as.data.frame(motif.family)
  rownames(temp) <- motif.vector
  inverted.master.dictionary <- rbind(inverted.master.dictionary,temp)
  }
```
2. Annotate motif using the consensus motif generated by TRF as reference 
```R
data <- data.frame(data)
data$motif <- as.character(data$consensus)
```
3. Reduce 7bp motifs to 6bp motifs to assign a motif family, using the inverted master dictionary
```R
dictionary.idx <- as.character(rownames(inverted.master.dictionary))
data$motif[which(data$consensus.size == 7)] <- 
  sapply(data$motif[which(data$consensus.size == 7)], function(x) 
    dictionary.idx[which.min(adist(x,dictionary.idx))])
```
4. Remove any entries with motifs larger than 7bp
```R
data <- data[which(!(data$consensus.size >= 8)),]
```
5. Annotate motif family
```R
data$motif.family <- NA
data$motif.family <- sapply(data$motif,function(x)
  as.character(inverted.master.dictionary[which(rownames(inverted.master.dictionary) == x),1]))
data[[names(data)[ii]]] <- data
```
## Step 5: Calculate & annotate uninterrupted length and uninterrupted copy number of each microsatellite
```R
#get a set of motifs and format them to search for maximum uninterrupted length
#Note: we need to use the consensus motif, rather than the reduced motif or motif family, since that is what TRF used to call matches
motifs <- mapply(function(x) paste("(",x,")","+",sep = ""),x=data$consensus)
#get maximum uninterrupted repetition of each microsatellite's motif
length.uninterrupted <- mapply(function(x,y) max(attr(gregexpr(x,y)[[1]],"match.length")),
                               x = motifs,
                               y = data$sequence)

data$length.uninterrupted <- length.uninterrupted
data$length.uninterrupted[which(data$length.uninterrupted == -1)] <- 0

#divide by the length of the motif to get the uninterrupted copy number
uninterrupted.copy.num <- data$length.uninterrupted/data$consensus.size

data$uninterrupted.copy.num <- uninterrupted.copy.num
```

## Step 6: Annotate 150bp 5' and 3' flanking regions for each microsatellite and get sequences for the 60bp and 90bp 5' and 3' flanking regions for future calculations
1. Define the size of the flanking sequences we would like to print in our dataframe
```R
n = 150
```
2. Define the sizes of the flanks that we would like to use for annotation calculations
```R
flank.calc.1 = 90
flank.calc.2 = 60
```
3. Annotate 5' flanking sequence:
```R
#get the first microsatellite entries from chromosomes that end before full 5' flank can be captured
chrom.ends.5.prime <- data.frame(seqnames = NULL, end = NULL)
for(ii in unique(data$seqnames)){
  first.usat.start <- min(data[which(data$seqnames == ii),]$start)
  coords <- c(ii, first.usat.start)
  if(first.usat.start-n < 1){
    chrom.ends.5.prime <- rbind(chrom.ends.5.prime,coords)
  }
}
colnames(chrom.ends.5.prime) <- c("seqnames","start")
chrom.ends.5.prime$start <- as.numeric(chrom.ends.5.prime$start)

#get all other microsatellites
usat.ends <- data[,c("seqnames","end")]
usat.ends <- anti_join(usat.ends,chrom.ends.5.prime, by = c("seqnames","end"))

#calculate coordinates of 5' flank for microsatellites that are not the last microsatellites on a chromosome that ends before the full 5' flank of length n can be captured.
flank.5.prime <- data.frame(data$seqnames,data$start-n,data$start-1)

#add coordinates of 5' flank to the end of the chromosome for microsatellites that are the last microsatellites on a chromosome that ends before the full flank of length n can be captured.
flank.5.prime <- rbind(flank.5.prime,c(chrom.ends.5.prime$seqnames, 1, chrom.ends.5.prime$start-1))

#get 5' flank sequence annotation
colnames(flank.5.prime) <- c("seqnames","start","end")
flank.5.prime <- makeGRangesFromDataFrame(flank.5.prime)
flank.5.prime.seq <- getSeq(Hsapiens, names = flank.5.prime)
data$flank.5.prime.seq <- as.character(flank.5.prime.seq)
```
4. Get but don't annotate 90bp 5' flanking sequence for future annotations:
```R
#get the first microsatellite entries from chromosomes that end before full 5' flank can be captured
chrom.ends.5.prime <- data.frame(seqnames = NULL, end = NULL)
for(ii in unique(data$seqnames)){
  first.usat.start <- min(data[which(data$seqnames == ii),]$start)
  coords <- c(ii, first.usat.start)
  if(first.usat.start-flank.calc.1 < 1){
    chrom.ends.5.prime <- rbind(chrom.ends.5.prime,coords)
  }
}
colnames(chrom.ends.5.prime) <- c("seqnames","start")
chrom.ends.5.prime$start <- as.numeric(chrom.ends.5.prime$start)

#get all other microsatellites
usat.ends <- data[,c("seqnames","end")]
usat.ends <- anti_join(usat.ends,chrom.ends.5.prime, by = c("seqnames","end"))

#calculate coordinates of 5' flank for microsatellites that are not the last microsatellites on a chromosome that ends before the full 5' flank of length n can be captured.
flank.5.prime.90bp <- data.frame(data$seqnames,data$start-flank.calc.1,data$start-1)

#add coordinates of 5' flank to the end of the chromosome for microsatellites that are the last microsatellites on a chromosome that ends before the full flank of length n can be captured.
flank.5.prime.90bp <- rbind(flank.5.prime.90bp,c(chrom.ends.5.prime$seqnames, 1, chrom.ends.5.prime$start-1))

#get 5' flank sequence
colnames(flank.5.prime.90bp) <- c("seqnames","start","end")
flank.5.prime.90bp <- makeGRangesFromDataFrame(flank.5.prime.90bp)
flank.5.prime.90bp.seq <- getSeq(Hsapiens, names = flank.5.prime.90bp)
```
5. Get but don't annotate 60bp 5' flanking sequence for future annotations:
```R
#get the first microsatellite entries from chromosomes that end before full 5' flank can be captured
chrom.ends.5.prime <- data.frame(seqnames = NULL, end = NULL)
for(ii in unique(data$seqnames)){
  first.usat.start <- min(data[which(data$seqnames == ii),]$start)
  coords <- c(ii, first.usat.start)
  if(first.usat.start-flank.calc.2 < 1){
    chrom.ends.5.prime <- rbind(chrom.ends.5.prime,coords)
  }
}
colnames(chrom.ends.5.prime) <- c("seqnames","start")
chrom.ends.5.prime$start <- as.numeric(chrom.ends.5.prime$start)

#get all other microsatellites
usat.ends <- data[,c("seqnames","end")]
usat.ends <- anti_join(usat.ends,chrom.ends.5.prime, by = c("seqnames","end"))

#calculate coordinates of 5' flank for microsatellites that are not the last microsatellites on a chromosome that ends before the full 5' flank of length n can be captured.
flank.5.prime.60bp <- data.frame(data$seqnames,data$start-flank.calc.2,data$start-1)

#add coordinates of 5' flank to the end of the chromosome for microsatellites that are the last microsatellites on a chromosome that ends before the full flank of length n can be captured.
flank.5.prime.60bp <- rbind(flank.5.prime.60bp,c(chrom.ends.5.prime$seqnames, 1, chrom.ends.5.prime$start-1))

#get 5' flank sequence
colnames(flank.5.prime.60bp) <- c("seqnames","start","end")
flank.5.prime.60bp <- makeGRangesFromDataFrame(flank.5.prime.60bp)
flank.5.prime.60bp.seq <- getSeq(Hsapiens, names = flank.5.prime.60bp)
```
6. Annotate 3' flanking sequence:
```R
#get the last microsatellite entries from chromosomes that end before full 3' flank can be captured
chrom.ends.3.prime <- data.frame(seqnames = NULL, end = NULL)
for(ii in unique(data$seqnames)){
  last.usat.end <- max(data[which(data$seqnames == ii),]$end)
  chrom.length <- length(Hsapiens[[ii]])
  coords <- c(ii, last.usat.end, chrom.length)
  if(chrom.length < last.usat.end+n){
    chrom.ends.3.prime <- rbind(chrom.ends.3.prime,coords)
  }
}
colnames(chrom.ends.3.prime) <- c("seqnames","end","chrom.length")
chrom.ends.3.prime$end <- as.numeric(chrom.ends.3.prime$end)
chrom.ends.3.prime$chrom.length <- as.numeric(chrom.ends.3.prime$chrom.length)

#get all other microsatellites
usat.ends <- data[,c("seqnames","end")]
usat.ends <- anti_join(usat.ends,chrom.ends.3.prime, by = c("seqnames","end"))

#calculate coordinates of 3' flank for microsatellites that are not the last microsatellites on a chromosome that ends before the full flank of length n can be captured.
flank.3.prime <- data.frame(usat.ends$seqnames,usat.ends$end+1,usat.ends$end+n)
#add coordinates of 3' flank to the end of the chromosome for microsatellites that are the last microsatellites on a chromosome that ends before the full flank of length n can be captured.
flank.3.prime <- rbind(flank.3.prime,c(chrom.ends.3.prime$seqnames, chrom.ends.3.prime$end+1, chrom.ends.3.prime$chrom.length))

#get 3' flank sequence annotation
colnames(flank.3.prime) <- c("seqnames","start","end")
flank.3.prime <- makeGRangesFromDataFrame(flank.3.prime)
flank.3.prime.seq <- getSeq(Hsapiens, names = flank.3.prime)
data$flank.3.prime.seq <- as.character(flank.3.prime.seq)
```
7. Get but don't annotate 90bp 3' flanking sequence for future annotations:
```R
#get the last microsatellite entries from chromosomes that end before full 3' flank can be captured
chrom.ends.3.prime <- data.frame(seqnames = NULL, end = NULL)
for(ii in unique(data$seqnames)){
  last.usat.end <- max(data[which(data$seqnames == ii),]$end)
  chrom.length <- length(Hsapiens[[ii]])
  coords <- c(ii, last.usat.end, chrom.length)
  if(chrom.length < last.usat.end+flank.calc.1){
    chrom.ends.3.prime <- rbind(chrom.ends.3.prime,coords)
  }
}
colnames(chrom.ends.3.prime) <- c("seqnames","end","chrom.length")
chrom.ends.3.prime$end <- as.numeric(chrom.ends.3.prime$end)
chrom.ends.3.prime$chrom.length <- as.numeric(chrom.ends.3.prime$chrom.length)

#get all other microsatellites
usat.ends <- data[,c("seqnames","end")]
usat.ends <- anti_join(usat.ends,chrom.ends.3.prime, by = c("seqnames","end"))

#calculate coordinates of 3' flank for microsatellites that are not the last microsatellites on a chromosome that ends before the full flank of length n can be captured.
flank.3.prime.90bp<- data.frame(usat.ends$seqnames,usat.ends$end+1,usat.ends$end+flank.calc.1)
#add coordinates of 3' flank to the end of the chromosome for microsatellites that are the last microsatellites on a chromosome that ends before the full flank of length n can be captured.
flank.3.prime.90bp<- rbind(flank.3.prime,c(chrom.ends.3.prime$seqnames, chrom.ends.3.prime$end+1, chrom.ends.3.prime$chrom.length))

#get 3' flank sequence annotation
colnames(flank.3.prime) <- c("seqnames","start","end")
flank.3.prime.90bp<- makeGRangesFromDataFrame(flank.3.prime)
flank.3.prime.seq <- getSeq(Hsapiens, names = flank.3.prime)
data$flank.3.prime.seq <- as.character(flank.3.prime.seq)
```
8. Get but don't annotate 60bp 3' flanking sequence for future annotations:
```R
#get the last microsatellite entries from chromosomes that end before full 3' flank can be captured
chrom.ends.3.prime <- data.frame(seqnames = NULL, end = NULL)
for(ii in unique(data$seqnames)){
  last.usat.end <- max(data[which(data$seqnames == ii),]$end)
  chrom.length <- length(Hsapiens[[ii]])
  coords <- c(ii, last.usat.end, chrom.length)
  if(chrom.length < last.usat.end+flank.calc.2){
    chrom.ends.3.prime <- rbind(chrom.ends.3.prime,coords)
  }
}
colnames(chrom.ends.3.prime) <- c("seqnames","end","chrom.length")
chrom.ends.3.prime$end <- as.numeric(chrom.ends.3.prime$end)
chrom.ends.3.prime$chrom.length <- as.numeric(chrom.ends.3.prime$chrom.length)

#get all other microsatellites
usat.ends <- data[,c("seqnames","end")]
usat.ends <- anti_join(usat.ends,chrom.ends.3.prime, by = c("seqnames","end"))

#calculate coordinates of 3' flank for microsatellites that are not the last microsatellites on a chromosome that ends before the full flank of length n can be captured.
flank.3.prime.60bp<- data.frame(usat.ends$seqnames,usat.ends$end+1,usat.ends$end+flank.calc.2)
#add coordinates of 3' flank to the end of the chromosome for microsatellites that are the last microsatellites on a chromosome that ends before the full flank of length n can be captured.
flank.3.prime.60bp<- rbind(flank.3.prime,c(chrom.ends.3.prime$seqnames, chrom.ends.3.prime$end+1, chrom.ends.3.prime$chrom.length))

#get 3' flank sequence annotation
colnames(flank.3.prime) <- c("seqnames","start","end")
flank.3.prime.60bp<- makeGRangesFromDataFrame(flank.3.prime)
flank.3.prime.seq <- getSeq(Hsapiens, names = flank.3.prime)
data$flank.3.prime.seq <- as.character(flank.3.prime.seq)
```
## Step 7: Calculate & annotate the percent of the microsatellite body that is overlapped by other microsatellites. Also calculate & annotate the percent of the 5' and 3' flanks of the microsatellite that are overlapped by other microsatellites.
1. Calculate the percent of the microsatellite body that is overlapped by other microsatellites:
```R
#make a granges object from the microsatellites dataframe
data.grange <- makeGRangesFromDataFrame(data)
#use findOverlaps to get the coordinates of regions within the granges object of our microsatellites dataframe that overlap with other regions in the granges object
hits <- findOverlaps(data.grange, ignore.strand=TRUE)
#exclude hits that have the same query and subject coordinates (i.e. where the microsatellite overlaps with itself)
hits <- hits[queryHits(hits)!=subjectHits(hits)]
#use pintersect to generate a granges object with only the intersecting regions between the microsatellites that overlap
overlaps <- pintersect(data.grange[queryHits(hits)], data.grange[subjectHits(hits)])
#divide the width of the overlap region by the width of the entire subject microsatellite to get the percent of the subject microsatellite that is overlapped by another microsatellite
percentOverlaps <- width(overlaps) / width(data.grange[subjectHits(hits)])
#convert the hits object to a dataframe
hits.df <- data.frame(hits)
```
2. Annotate the percent of the microsatellite body that is overlapped by other microsatellites:
```R
#initialize a column to store microsatellite overlap annotations. Fill it with zeros.
data$msat.overlap <- 0
#replace the rows that match the subjectHits coordinates (the second column in the hits.df dataframe) with the calculated percent overlap for those microsatellites
data$msat.overlap[hits.df[,2]] <- percentOverlaps
```
3. Calculate the percent of the 90bp 5' flank of the microsatellite that is overlapped by other microsatellites:
```R
hits <- findOverlaps(flank.5.prime.90bp, data.grange, ignore.strand=TRUE)
overlaps <- pintersect(flank.5.prime.90bp[queryHits(hits)], data.grange[subjectHits(hits)])
percentOverlaps <- width(overlaps) / width(flank.5.prime.90bp[queryHits(hits)])
```
4. Annotate the percent of the 90bp 5' flank of the microsatellite that is overlapped by other microsatellites:
```R
hits.df <- data.frame(hits)
data$overlap.5.prime.90bp <- 0
data$overlap.5.prime.90bp[hits.df[,1]] <- percentOverlaps
```
5. Calculate the percent of the 60bp 5' flank of the microsatellite that is overlapped by other microsatellites:
```R
hits <- findOverlaps(flank.5.prime.60bp, data.grange, ignore.strand=TRUE)
overlaps <- pintersect(flank.5.prime.60bp[queryHits(hits)], data.grange[subjectHits(hits)])
percentOverlaps <- width(overlaps) / width(flank.5.prime.60bp[queryHits(hits)])
```
6. Annotate the percent of the 60bp 5' flank of the microsatellite that is overlapped by other microsatellites:
```R
hits.df <- data.frame(hits)
data$overlap.5.prime.60bp <- 0
data$overlap.5.prime.60bp[hits.df[,1]] <- percentOverlaps
```
7. Calculate the percent of the 90bp 3' flank of the microsatellite that is overlapped by other microsatellites:
```R
hits <- findOverlaps(flank.3.prime.90bp, data.grange, ignore.strand=TRUE)
overlaps <- pintersect(flank.3.prime.90bp[queryHits(hits)], data.grange[subjectHits(hits)])
percentOverlaps <- width(overlaps) / width(flank.3.prime.90bp[queryHits(hits)])
```
8. Annotate the percent of the 90bp 3' flank of the microsatellite that is overlapped by other microsatellites:
```R
hits.df <- data.frame(hits)
data$overlap.3.prime.90bp <- 0
data$overlap.3.prime.90bp[hits.df[,1]] <- percentOverlaps
```
9. Calculate the percent of the 60bp 3' flank of the microsatellite that is overlapped by other microsatellites:
```R
hits <- findOverlaps(flank.3.prime.60bp, data.grange, ignore.strand=TRUE)
overlaps <- pintersect(flank.3.prime.60bp[queryHits(hits)], data.grange[subjectHits(hits)])
percentOverlaps <- width(overlaps) / width(flank.3.prime.60bp[queryHits(hits)])
```
10. Annotate the percent of the 60bp 3' flank of the microsatellite that is overlapped by other microsatellites:
```R
hits.df <- data.frame(hits)
data$overlap.3.prime.60bp <- 0
data$overlap.3.prime.60bp[hits.df[,1]] <- percentOverlaps
```
## Step 8: Calculate & annotate the distance to the next nearest microsatellite upstream and downstream
```R
#make granges object from microsatellites dataframe
data.granges <- makeGRangesFromDataFrame(data, keep.extra.columns = T)
#call upstream hits by using the follow function in the GenomicRanges package. "follow" returns the index of the range in the subject (data.granges) that is directly followed by the range being evaluated. Remove any indices that do not have a preceding range listed.
upstream.hits <- data.granges[follow(data.granges)[-c(which(is.na(follow(data.granges))))]]
#get the regions in the dataset that have regions that preceed them
data.upstream <- data.granges[-c(which(is.na(follow(data.granges))))]
#use the distance function from the GenomicRanges package to calculate the distance between each microsatellite (data.upstream) and the microsatellite directly preceeding it (upstream.hits)
dist.upstream <- GenomicRanges::distance(data.upstream,upstream.hits)
#initialize column for nearest upstream microsatellite with NA in each row
data$nearest.usat.5 <- NA
#fill in the distance to the next nearest upstream microsatellite for entries that have microsatellites that preceed them.
data$nearest.usat.5[-c(which(is.na(follow(data.granges))))] <- as.numeric(dist.upstream)

#repeat the same process as above for microsatellites that follow each microsatellite using the "precede" function in the GenomicRanges package
downstream.hits <- data.granges[precede(data.granges)[-c(which(is.na(precede(data.granges))))]]
data.downstream <- data.granges[-c(which(is.na(precede(data.granges))))]
dist.downstream <- GenomicRanges::distance(data.downstream,downstream.hits)
data$nearest.usat.3 <- NA
data$nearest.usat.3[-c(which(is.na(precede(data.granges))))] <- as.numeric(dist.downstream)
```
## Step 9: Calculate & annotate the % GC content for the microsatellite motif, the microsatellite body, the 5' flank of the microsatellite, and the 3' flank of the microsatellite
```R
#use the letetrFrequency function from the Biostrings package to get the number of times the letters "G" and "C" occur in each microsatellite's motif family. Divide that by the length of the motif family to get the % GC content of the motif for each microsatellite
data$GC.motif <- letterFrequency(DNAStringSet(data$motif.family),
                                                "GC")/nchar(data$motif.family)

#repeat this process for the actual microsatellite sequence
data$GC.msat <- letterFrequency(DNAStringSet(data$sequence), 
                                               "GC")/data$width
flank.calc.1 = 90
flank.calc.2 = 60

#repeat this process for the 90bp 5' flanking sequence of the microsatellite
data$GC.90bp.flank.5 <-
  letterFrequency(DNAStringSet(flank.5.prime.90bp.seq), 
                                               "GC")/flank.calc.1

#repeat this process for the 90bp 3' flanking sequence of the microsatellite
data$GC.90bp.flank.3 <-
  letterFrequency(DNAStringSet(flank.3.prime.90bp.seq), 
                                               "GC")/flank.calc.1

#repeat this process for the 60bp 5' flanking sequence of the microsatellite
data$GC.60bp.flank.5 <-
  letterFrequency(DNAStringSet(flank.5.prime.60bp.seq), 
                                               "GC")/flank.calc.2

#repeat this process for the 60bp 3' flanking sequence of the microsatellite
data$GC.60bp.flank.3 <-
  letterFrequency(DNAStringSet(flank.3.prime.60bp.seq), 
                                               "GC")/flank.calc.2
```
## Step 10: Sort and export microsatellite body and 60bp and 90bp 5' and 3' flank regions as BED files for further calculations in linux
1. Sort microsatellite dataframe and replace unsorted dataframe
```R
data.sorted <- data[with(data,
                         order(data$seqnames,
                               as.numeric(data$start),
                               as.numeric(data$end))),]
rm(data)
```
2. Export microsatellite data as bed file for use with bedtools/UCSC tools functions for further annotations
```R
#export microsatellite data as a bed file
data.granges <- makeGRangesFromDataFrame(data.sorted)
names(data.granges) <- 1:nrow(data.sorted)
export(data.granges,con = "[FILENAME].bed",format = "BED")
```
3. Export 90bp 5' flank data as bed file for use with bedtools/UCSC tools functions for further annotations
```R
flank.calc.1 = 90

#get the first microsatellite entries from chromosomes that end before full 5' flank can be captured
chrom.ends.5.prime <- data.frame(seqnames = NULL, end = NULL)
for(ii in unique(data.sorted$seqnames)){
  first.usat.start <- min(data.sorted[which(data.sorted$seqnames == ii),]$start)
  coords <- c(ii, first.usat.start)
  if(first.usat.start-flank.calc.1 < 1){
    chrom.ends.5.prime <- rbind(chrom.ends.5.prime,coords)
  }
}
colnames(chrom.ends.5.prime) <- c("seqnames","start")
chrom.ends.5.prime$start <- as.numeric(chrom.ends.5.prime$start)

#get all other microsatellites
usat.ends <- data.sorted[,c("seqnames","end")]
usat.ends <- anti_join(usat.ends,chrom.ends.5.prime, by = c("seqnames","end"))

#calculate coordinates of 5' flank for microsatellites that are not the last microsatellites on a chromosome that ends before the full 5' flank of length n can be captured.
flank.5.prime.90bp <- data.frame(data.sorted$seqnames,data.sorted$start-flank.calc.1,data.sorted$start-1)

#add coordinates of 5' flank to the end of the chromosome for microsatellites that are the last microsatellites on a chromosome that ends before the full flank of length n can be captured.
flank.5.prime.90bp <- rbind(flank.5.prime.90bp,c(chrom.ends.5.prime$seqnames, 1, chrom.ends.5.prime$start-1))

#export 90bp 5' flank
flank.5.prime.90bp.granges <- 
  makeGRangesFromDataFrame(flank.5.prime.90bp)
names(flank.5.prime.90bp.granges) <- 1:nrow(data.sorted)
export(flank.5.prime.90bp.granges,
       con = “[FILENAME].bed,format = "BED")
```
4. Export 60bp 5' flank data as bed file for use with bedtools/UCSC tools functions for further annotations
```R
flank.calc.2 = 60

#get the first microsatellite entries from chromosomes that end before full 5' flank can be captured
chrom.ends.5.prime <- data.frame(seqnames = NULL, end = NULL)
for(ii in unique(data.sorted$seqnames)){
  first.usat.start <- min(data.sorted[which(data.sorted$seqnames == ii),]$start)
  coords <- c(ii, first.usat.start)
  if(first.usat.start-flank.calc.2 < 1){
    chrom.ends.5.prime <- rbind(chrom.ends.5.prime,coords)
  }
}
colnames(chrom.ends.5.prime) <- c("seqnames","start")
chrom.ends.5.prime$start <- as.numeric(chrom.ends.5.prime$start)

#get all other microsatellites
usat.ends <- data.sorted[,c("seqnames","end")]
usat.ends <- anti_join(usat.ends,chrom.ends.5.prime, by = c("seqnames","end"))

#calculate coordinates of 5' flank for microsatellites that are not the last microsatellites on a chromosome that ends before the full 5' flank of length n can be captured.
flank.5.prime.60bp <- data.frame(data.sorted$seqnames,data.sorted$start-flank.calc.2,data.sorted$start-1)

#add coordinates of 5' flank to the end of the chromosome for microsatellites that are the last microsatellites on a chromosome that ends before the full flank of length n can be captured.
flank.5.prime.60bp <- rbind(flank.5.prime.60bp,c(chrom.ends.5.prime$seqnames, 1, chrom.ends.5.prime$start-1))

#export 60bp 5' flank
flank.5.prime.60bp.granges <- 
  makeGRangesFromDataFrame(flank.5.prime.60bp)
names(flank.5.prime.60bp.granges) <- 1:nrow(data.sorted)
export(flank.5.prime.60bp.granges,
       con = “[FILENAME].bed,format = "BED")
```
5. Export 90bp 3' flank data as bed file for use with bedtools/UCSC tools functions for further annotations
```R
flank.calc.1 = 90

#get the first microsatellite entries from chromosomes that end before full 3' flank can be captured
chrom.ends.3.prime <- data.frame(seqnames = NULL, end = NULL)
for(ii in unique(data.sorted$seqnames)){
  first.usat.start <- min(data.sorted[which(data.sorted$seqnames == ii),]$start)
  coords <- c(ii, first.usat.start)
  if(first.usat.start-flank.calc.1 < 1){
    chrom.ends.3.prime <- rbind(chrom.ends.3.prime,coords)
  }
}
colnames(chrom.ends.3.prime) <- c("seqnames","start")
chrom.ends.3.prime$start <- as.numeric(chrom.ends.3.prime$start)

#get all other microsatellites
usat.ends <- data.sorted[,c("seqnames","end")]
usat.ends <- anti_join(usat.ends,chrom.ends.3.prime, by = c("seqnames","end"))

#calculate coordinates of 3' flank for microsatellites that are not the last microsatellites on a chromosome that ends before the full 3' flank of length n can be captured.
flank.3.prime.90bp <- data.frame(data.sorted$seqnames,data.sorted$start-flank.calc.1,data.sorted$start-1)

#add coordinates of 3' flank to the end of the chromosome for microsatellites that are the last microsatellites on a chromosome that ends before the full flank of length n can be captured.
flank.3.prime.90bp <- rbind(flank.3.prime.90bp,c(chrom.ends.3.prime$seqnames, 1, chrom.ends.3.prime$start-1))

#export 90bp 3' flank
flank.3.prime.90bp.granges <- 
  makeGRangesFromDataFrame(flank.3.prime.90bp)
names(flank.3.prime.90bp.granges) <- 1:nrow(data.sorted)
export(flank.3.prime.90bp.granges,
       con = “[FILENAME].bed,format = "BED")
```
6. Export 60bp 3' flank data as bed file for use with bedtools/UCSC tools functions for further annotations
```R
flank.calc.2 = 60

#get the first microsatellite entries from chromosomes that end before full 3' flank can be captured
chrom.ends.3.prime <- data.frame(seqnames = NULL, end = NULL)
for(ii in unique(data.sorted$seqnames)){
  first.usat.start <- min(data.sorted[which(data.sorted$seqnames == ii),]$start)
  coords <- c(ii, first.usat.start)
  if(first.usat.start-flank.calc.2 < 1){
    chrom.ends.3.prime <- rbind(chrom.ends.3.prime,coords)
  }
}
colnames(chrom.ends.3.prime) <- c("seqnames","start")
chrom.ends.3.prime$start <- as.numeric(chrom.ends.3.prime$start)

#get all other microsatellites
usat.ends <- data.sorted[,c("seqnames","end")]
usat.ends <- anti_join(usat.ends,chrom.ends.3.prime, by = c("seqnames","end"))

#calculate coordinates of 3' flank for microsatellites that are not the last microsatellites on a chromosome that ends before the full 3' flank of length n can be captured.
flank.3.prime.60bp <- data.frame(data.sorted$seqnames,data.sorted$start-flank.calc.2,data.sorted$start-1)

#add coordinates of 3' flank to the end of the chromosome for microsatellites that are the last microsatellites on a chromosome that ends before the full flank of length n can be captured.
flank.3.prime.60bp <- rbind(flank.3.prime.60bp,c(chrom.ends.3.prime$seqnames, 1, chrom.ends.3.prime$start-1))

#export 60bp 3' flank
flank.3.prime.60bp.granges <- 
  makeGRangesFromDataFrame(flank.3.prime.60bp)
names(flank.3.prime.60bp.granges) <- 1:nrow(data.sorted)
export(flank.3.prime.60bp.granges,
       con = “[FILENAME].bed,format = "BED")
```
## Step 11: Generate & annotate mappability of each microsatellite using the UMAP k24, k36, k50, and k100, available here: https://bismap.hoffmanlab.org/
1. [optional] There are pre-made bedgraph files mapped to the hg38, hg19, and mm10 reference genomes. If you need to use another reference genome, you will need to use the UCSC ```liftOver``` tool to liftover the annotations to the reference genome of your choice.
2. Unzip and sort bedgraph files for use
```bash
gunzip -c k24.umap.bedgraph.gz | sort -k1,1 -k2,2n > k24.umap.sorted.bedgraph"
gunzip -c k36.umap.bedgraph.gz | sort -k1,1 -k2,2n > k36.umap.sorted.bedgraph"
gunzip -c k50.umap.bedgraph.gz | sort -k1,1 -k2,2n > k50.umap.sorted.bedgraph"
gunzip -c k100.umap.bedgraph.gz | sort -k1,1 -k2,2n > k100.umap.sorted.bedgraph"
```
3. Use ```bedGraphToBigWig``` tool from UCSC tools to convert each umap mappability bedgraph file to a bigwig format
```bash
bedGraphToBigWig k24.umap.sorted.bedgraph [genome] k24.umap.sorted.bw
bedGraphToBigWig k36.umap.sorted.bedgraph [genome] k36.umap.sorted.bw
bedGraphToBigWig k50.umap.sorted.bedgraph [genome] k50.umap.sorted.bw
bedGraphToBigWig k100.umap.sorted.bedgraph [genome] k100.umap.sorted.bw
```
4. Use the ```bigWigAverageOverBed``` function from the UCSC tools to calculate the average mappability across each microsatellite in the microsatellite bed file for each mappability track
```bash
bigWigAverageOverBed k24.umap.sorted.bw [FILENAME].bed [FILENAME].k24.umap.sorted.tab
bigWigAverageOverBed k36.umap.sorted.bw [FILENAME].bed [FILENAME].k36.umap.sorted.tab
bigWigAverageOverBed k50.umap.sorted.bw [FILENAME].bed [FILENAME].k50.umap.sorted.tab
bigWigAverageOverBed k100.umap.sorted.bw [FILENAME].bed [FILENAME].k100.umap.sorted.tab
```
5. Repeat step 4 to average over the 90bp 5' flank bed file
6. Repeat step 4 to average over the 60bp 5' flank bed file
7. Repeat step 4 to average over the 90bp 3' flank bed file
8. Repeat step 4 to average over the 60bp 3' flank bed file
9. Load microsatellite mappability annotations into R and add to the microsatellite annotations dataframe using the following commands:
```R
umap.k24.msat = system("awk '{print $5}' /Path/to/[FILENAME].k24.umap.sorted.tab", intern=T)
umap.k24.msat <- as.numeric(umap.k24.msat)
data.sorted$umap.k24.msat <- umap.k24.msat

umap.k36.msat = system("awk '{print $5}' /Path/to/[FILENAME].k36.umap.sorted.tab", intern=T)
umap.k36.msat <- as.numeric(umap.k36.msat)
data.sorted$umap.k36.msat <- umap.k36.msat

umap.k50.msat = system("awk '{print $5}' /Path/to/[FILENAME].k50.umap.sorted.tab", intern=T)
umap.k50.msat <- as.numeric(umap.k50.msat)
data.sorted$umap.k50.msat <- umap.k50.msat

umap.k100.msat = system("awk '{print $5}' /Path/to/[FILENAME].k100.umap.sorted.tab", intern=T)
umap.k100.msat <- as.numeric(umap.k100.msat)
data.sorted$umap.k100.msat <- umap.k100.msat
```
10. Repeat step 9 to add the mappability annotations for the 90bp 5' flank
11. Repeat step 9 to add the mappability annotations for the 60bp 5' flank
12. Repeat step 9 to add the mappability annotations for the 90bp 3' flank
13. Repeat step 9 to add the mappability annotations for the 60bp 3' flank
## Step 12: Generate replication timing annotations for GM12878 lymphoblastoid cell line & annotate microsatellites with this data
1. Use ```nano``` command to write a script called "replace.zeros.awk"
2. copy the following code into the script and save:
```bash
#before doing anything, ensure that the field separator and the output field separator is set to "\t" (tab)
BEGIN { OFS = FS = "\t"}

#for all rows loop through all fields from 4 until you have reached the total number of fields in the file and if you encounter a 0, replace it with 0.005. (we start with field 4 because fields 1-3 will describe the region coordinates in a bed file (the file type with which we are working) and we do not want to change those)
NR > 0 {
    for (i = 4; i <= NF; ++i) {
        if ($i == "0") {
            $i = "0.005";
        }
    }
}
#print the edited datafile when finished
{ print }
```
3. Run ```chmod +x replace.zeros.awk``` to make the script executable
4. Download the 4D Nucleome data for the GM12878 lymphoblastoid cell line from <a href="https://data.4dnucleome.org/browse/?experimentset_type=replicate&type=ExperimentSetReplicate">here</a>
5. [optional] Note: these data are mapped to the hg38 reference genome only. To use a different reference genome, you will need to liftover the annotations to your preferred reference genome using the UCSC ```liftOver``` tool before preceeding to the next step.
6. Sort datafiles containing replication timing data and rename them intuitively
```bash
sort -k1,1 -k2,2n 4DNFI6TILWWX.bedGraph > 4DN-early-rep1_sorted.bedGraph
sort -k1,1 -k2,2n 4DNFIIMJQ8NT.bedGraph > 4DN-early-rep2_sorted.bedGraph
sort -k1,1 -k2,2n 4DNFIS7J9B9X.bedGraph > 4DN-late-rep1_sorted.bedGraph
sort -k1,1 -k2,2n 4DNFIT26294Y.bedGraph > 4DN-late-rep2_sorted.bedGraph
```
7. Use the unionbedg function from bedtools to combine all of the replication datafiles into a single file with a separate field (column) for each data entry. Fill in the regions of the genome without data in one or multiple fields with zeros for the missing fields. Pipe that output into the replace.zeros.awk script that we wrote to replace any zeros with 0.005 so that we can perform log calculations with these data. Then, add two columns onto the datatable that takes the average of the two replicate columnns for the early and late samples, then add another column that takes the log2(Early/Late) samples. Finally combine just the first three columns (the region coordinates) and the last column (the log2(Early/Late) samples) and save it as a bedGraph file.
```bash
bedtools unionbedg -i 4DN-late-rep1_sorted.bedGraph 4DN-late-rep2_sorted.bedGraph 4DN-early-rep1_sorted.bedGraph 4DN-early-rep2_sorted.bedGraph -empty -g [CHROM LENGTHS FILE] | awk -f replace.zeros.awk | awk '{print $0 "\t" ($4+$5)/2 "\t" ($6+$7)/2}' | awk '{print $0 "\t" log($9/$8)/log(2)}' | awk '{print $1 "\t" $2 "\t" $3 "\t" $10}' > 4DN-GM12878-RT-annotation_sorted.bedGraph
```
8. Convert processed replication timing datafile from a bedGraph to a bigWig using ```bedGraphToBigWig``` from the UCSC tools
```bash
bedGraphToBigWig 4DN-GM12878-RT-annotation_sorted.bedGraph [genome] 4DN-GM12878-RT-annotation_sorted.bw
```
9. Use ```bigWigAverageOverBed``` to find the average replication timing for each microsatellite in our microsatellite dataset
```bash
bigWigAverageOverBed 4DN-GM12878-RT-annotation_sorted.bw [FILENAME].bed [FILENAME].4DN.RT.sorted.tab
```
10. In R, add 4D nucleome replication timing annotation to the microsatellite annotations dataframe:
```R
reptiming.4DN = 
  system("awk '{print $5}'Path/To/[FILENAME].4DN.RT.sorted.tab", intern=T)
reptiming.4DN <- as.numeric(reptiming.4DN)
data.sorted$reptiming.4DN <- reptiming.4DN
```
## Step 13: Generate replication timing annotations for neural progenitor cells using repli-chip data from the <a href="https://www.encodeproject.org/">ENCODE Project</a>
The following methods are based on <a href="https://www.nature.com/articles/nprot.2017.148#Sec12">https://www.nature.com/articles/nprot.2017.148#Sec12</a>, procedure section, data analysis subsection, starting from step 81 A Analyzing single-end data part iii (for multiple samples)
1. Download the repli-chip data from the ENCODE Project <a href="https://www.encodeproject.org/experiments/ENCSR497CCF/">here</a>
2. [optional] If using a reference genome other than hg19, do the following steps to liftover the annotations to your reference genome of choice. If using hg19 as your reference genome, skip to step 6.
3. Convert bigwig to bedgraph format using UCSC ```bigWigToBedGraph``` tool.
4. Sort bedgraph file using ```sort -k1,1 -k2,2n [FILENAME].bedgraph```.
5. Run sorted bedgraph through the UCSC ```liftOver``` tool.
6. Use the ```bedClip``` tool from UCSC tools to remove any data that are called in regions of the genome that do not exist
7. Use the bedtools ```unionbedg``` function to combine both annotations into a single .txt file and fill in missing data in either or both field with "NA"
```bash
bedtools unionbedg -i [FILENAME1].sorted.clipped.bedGraph [FILENAME2].sorted.clipped.bedGraph -empty -g [CHROM LENGTH] -filler "NA" > merge_RT.txt
```
8. Next, in R, normalize the data between replicates and run a Loess smoothing (Locally Weighted Scatterplot Smoothing) across each chromosome, so that almost all loci in the genome have replication timing estimates
```R
#load necessary library for data normalizing and smoothing
library(preprocessCore)
#read in datatable generated in bash
merge <- read.table("/Path/To/merge_RT.txt",header=FALSE)
#get only columns that contain replication timing data (i.e. not region coordinates), convert them to a matrix, and store them in a merge_values variable
merge_values <- as.matrix(merge[,4:ncol(merge)])
#get the values from each column in our merge matrix (containing just the replication timing data) and reformat it so that it is just a single column of values
ad <- stack(merge[,4:ncol(merge)])$values
#normalize the replication timing data for each replicate separately based on the distribution of values of the combined replicates replication timing (ad variable)
norm_data <- normalize.quantiles.use.target(merge_values,ad)
#combine the normalized data back with the coordinates of the data regions
merge_norm <- data.frame(merge[,1:3],norm_data)
#rename the columns of the normalized dataframe
colnames(merge_norm) <- c("chr","start","stop","ENCFF469TYS.[REF GENOME].sorted.bedGraph","ENCFF923TOC.[REF GENOME].sorted.bedGraph")

#output quantile-normalized data by going through each column that records replication timing for each replicate and printing the first three columns and the column containing the replication timing data for that replicate for rows that have no missing data. Put that dataset into a bedGraph file with the extension ".qnorm.bedGraph" to differentiate it from the non-normalized data
for(i in 4:ncol(merge_norm))
  {write.table(merge_norm[complete.cases(merge_norm[,i]), c(1,2,3,i)], 
               gsub(".bedGraph", ".qnorm.bedGraph", colnames(merge_norm)[i]), 
               sep= "\t",
               row.names=FALSE, 
               quote=FALSE, 
               col.names = FALSE)}

#get chromosomes except the unmapped chromosomes (containing “_” in their name) which can be problematic for Loess smoothing, and the chromosomes Y and M (mitochondrial)
chrs=grep(levels(as.factor(merge_norm$chr)),pattern= "[_YM]",invert=TRUE,value=TRUE)

#initialize list for Loess smoothed data
AllLoess=list()

#for each replicate (column after the first three columns, which contain locus coordinates) generate Loess-smoothed data by running the following loop
for(i in 1:(ncol(merge_norm)-3)){
  #intialize a dataframe as the first element of the AllLoess list
  AllLoess[[i]]=data.frame();
  #as you go through each replicate's dataset, print which dataset you are currently processing
  cat("Current dataset:",colnames(merge_norm)[i+3],"\n");
  #for each chromosome recorded in the dataset (excluding unmapped chromosomes, and chromosomes Y and M), subset the dataset based on that chromosome, fit a Loess regression (RTla) to the replication timing data using a span of 300kb/length of the given chromosome, and append the fitted values, based on that regression, to the region coordinates. Then, get the regions in the subset for which there is no observed replication timing (i.e. the replication timing is listed as "NA"), and use the Loess regression (RTla) to predict what the replication timing values would be based on the start location of the regions with missing replication timing. Finally, rename the columns of the missing values dataframe so that it can then be rebound to the regions that did not have missing value, and then add it to the AllLoess list.
  for(Chr in chrs){
    RTb=subset(merge_norm, merge_norm$chr==Chr);
    lspan= 300000/(max(RTb$start)-min(RTb$start));
    cat("Current chrom:", Chr, "\n");
    RTla=loess(RTb[,i+3] ~ RTb$start, span=lspan);
    
    RTl=data.frame(c(rep(Chr,times=RTla$n)), RTla$x,
                   merge_norm[which(merge_norm$chr==Chr & 
                                      merge_norm$start %in% RTla$x),3],RTla$fitted);
    
    missing.values <- subset(RTb, is.na(RTb[,i+3]));
    fit <- transform(missing.values$start, y.pred = predict(RTla, missing.values$start))
    missing.values[,i+3] <- fit$y.pred
    missing.values <- missing.values[which(!(is.na(missing.values[,i+3]))),c(1,2,3,i+3)]
    colnames(missing.values) <- colnames(RTl)
    RTl <- rbind(RTl,missing.values)
    
    colnames(RTl)=c("chr","start","end",colnames(RTb)[i+3]);
    if(length(AllLoess[[i]])!=0){
      AllLoess[[i]]=rbind(AllLoess[[i]],RTl)};
    if(length(AllLoess[[i]])==0){
      AllLoess[[i]]=RTl};
  }
}

#once you have Loess smoothed the replication timing data for all replicates, export the data with a new extension (".Loess.bedGraph" instead of ".bedGraph") so that we can intersect it with our microsatellite dataset and generate the microsatellite repliaction timing annotations in bash.
for(i in 1:length(AllLoess)){
  write.table(AllLoess[[i]][complete.cases(AllLoess[[i]]),], 
              gsub(".bedGraph", ".Loess.bedGraph", 
                   colnames(AllLoess[[i]]))[4], 
              sep= "\t", 
              row.names=FALSE, 
              quote=FALSE, 
              col.names = FALSE)}
```
9. Sort each Loess-smoothed replication timing bedgraph file by chromosome and then by start position using the ```sort -k1,1 -k2,2n``` command to ensure it is compatible with our microsatellite data
10. Convert each bedGraph file to a bigwig file so that we can the use the UCSC tool ```bigWigAverageOverBed``` to get the average replication timing annotation for each microsatellite in our dataset
```bash
bedGraphToBigWig ENCFF469TYS.[REF GENOME].Loess.sorted.bedGraph [genome] ENCFF469TYS.hg19.Loess.sorted.bw
bedGraphToBigWig ENCFF923TOC.[REF GENOME].Loess.sorted.bedGraph [genome] ENCFF923TOC.hg19.Loess.sorted.bw
```
11. Use the UCSC tool ```bigWigAverageOverBed``` to get the average replication timing annotation for each microsatellite in our dataset
```bash
bigWigAverageOverBed ENCFF469TYS.[REF GENOME].Loess.sorted.bw [FILENAME].sorted.bed [FILENAME].ENCFF469TYS.Loess.sorted.tab
bigWigAverageOverBed ENCFF923TOC.[REF GENOME].Loess.sorted.bw [FILENAME].sorted.bed [FILENAME].ENCFF923TOC.Loess.sorted.tab
```
12. In R, add NPC replication timing annotation to the microsatellite annotations dataframe:
```R
reptiming.NPC.ENCFF469TYS = 
  system("awk '{print $5}' Path/To/[FILENAME].ENCFF469TYS.Loess.sorted.tab", 
         intern=T)
reptiming.NPC.ENCFF469TYS <- as.numeric(reptiming.NPC.ENCFF469TYS)
data.sorted$reptiming.NPC.ENCFF469TYS <- reptiming.NPC.ENCFF469TYS

reptiming.NPC.ENCFF923TOC = 
  system("awk '{print $5}' Path/To/[FILENAME].ENCFF923TOC.Loess.sorted.tab", 
         intern=T)
reptiming.NPC.ENCFF923TOC <- as.numeric(reptiming.NPC.ENCFF923TOC)
data.sorted$reptiming.NPC.ENCFF923TOC <- reptiming.NPC.ENCFF923TOC
```
## Step 13: Annotate Bioskryb coverage of each microsatellite and its 60bp and 90bp 5' and 3' flanking regions
1. Download data from this GitHub repository.  Note: these data are mapped to the hg38 reference genome only. To use a different reference genome, you will need to liftover the annotations to your preferred reference genome using the <a href="https://academic.oup.com/bioinformatics/article/30/7/1006/234947?login=false">CrossMap</a> tool before preceeding to the next step.
2. Use the UCSC tool ```bigWigAverageOverBed``` to get the average Bioskryb coverage annotation for each microsatellite in our dataset
```bash
bigWigAverageOverBed BioSkryb.allsc.genomecov.chr1-22X.normalized.sorted.bw [FILENAME].sorted.bed [FILENAME].bioskryb.sorted.tab
```
3. Repeat step 2 to average over the 90bp 5' flank bed file
4. Repeat step 2 to average over the 60bp 5' flank bed file
5. Repeat step 2 to average over the 90bp 3' flank bed file
6. Repeat step 2 to average over the 60bp 3' flank bed file
7. In R, add Bioskryb coverage annotation to the microsatellite annotations dataframe:
```R
msat.bioskryb = 
  system("awk '{print $5}' ./[FILENAME].bioskryb.sorted.tab", 
         intern=T)
msat.bioskryb <- as.numeric(msat.bioskryb)
data.sorted$msat.bioskryb <- msat.bioskryb

flank.5.90bp.bioskryb = 
  system("awk '{print $5}' ./[FILENAME].flank.5.prime.90bp.bioskryb.sorted.tab", 
         intern=T)
flank.5.90bp.bioskryb <- as.numeric(flank.5.90bp.bioskryb)
data.sorted$flank.5.90bp.bioskryb <- flank.5.90bp.bioskryb

flank.5.60bp.bioskryb = 
  system("awk '{print $5}' ./[FILENAME].flank.5.prime.60bp.bioskryb.sorted.tab", 
         intern=T)
flank.5.60bp.bioskryb <- as.numeric(flank.5.60bp.bioskryb)
data.sorted$flank.5.60bp.bioskryb <- flank.5.60bp.bioskryb

flank.3.90bp.bioskryb = 
  system("awk '{print $5}' ./[FILENAME].flank.3.prime.90bp.bioskryb.sorted.tab", 
         intern=T)
flank.3.90bp.bioskryb <- as.numeric(flank.3.90bp.bioskryb)
data.sorted$flank.3.90bp.bioskryb <- flank.3.90bp.bioskryb

flank.3.60bp.bioskryb = 
  system("awk '{print $5}' ./[FILENAME].flank.3.prime.60bp.bioskryb.sorted.tab", 
         intern=T)
flank.3.60bp.bioskryb <- as.numeric(flank.3.60bp.bioskryb)
data.sorted$flank.3.60bp.bioskryb <- flank.3.60bp.bioskryb
```
## Step 14: Train & run the Gymrek lab's MUTEA model, using autosomal intergenic training data from the Gymrek lab and the annotations we generated that best explained the variance in the test dataset ("uninterrupted length" and "period size") to generate mutation rate estimates for each microsatellite
This step recapitulates the MUTEA linear model from this paper: <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5679271/">https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5679271/</a>
1. Download training data from <a href="https://github.com/gymreklab/mutea-autosomal/tree/master/data">https://github.com/gymreklab/mutea-autosomal/tree/master/data</a>
2. [Optional] These data are mapped to the hg19 reference genome only. To use a different reference genome, you will need to liftover the annotations to your preferred reference genome using the UCSC ```liftOver``` tool before preceeding to the next step.
3. Define filter 1, which filters microsatellites "with standard errors equal to 0, >0.1, or undefined (usually indicating the lower optimization boundary of 10−8 was reached)" (Gymrek et al., 2017)
```R
SetFilterModel1 <- function(df){
  df$colname <- sapply(df$ml_mu_stderr_adj, function(x) isTRUE(x <= 0 | x >= 0.1))
  return(df$colname)
}
```
4. Define filter parameters
```R
minlen = 20
features = c("uninterrupted_length","period")
```
5. Remove duplicates from training data
```bash
sort -k1,1 -k2,2n -k3,3n autosomal_perlocus_train_intergenic.bed | uniq > temp.bed
cp temp.bed autosomal_perlocus_train_intergenic.bed
rm temp.bed
```
6. load training data into R
```R
datafile = read.table("./autosomal_perlocus_train_intergenic.bed", header = T)
```
7. Filter training file based on adjusted standard error
```R
datafile$ml_mu_stderr_adj <- datafile$ml_mu_stderr/abs(datafile$ml_mu)
datafile <- datafile[which(datafile$ml_mu_stderr_adj != 0),]
```
8. Define new features in training data
```R
datafile$period <- nchar(as.character(datafile$motif))
datafile$pbyl <- datafile$uninterrupted_length*1.0/datafile$period
```
9. Restrict training data to di-, tri-, and tetramers. Also restrict based on "feature filter" and length of microsatellite
```R
datafile <- datafile[which(datafile$period %in% c(2,3,4,5,6)),]
datafile <- datafile[which(datafile$featurefilter == "False"),]
datafile <- datafile[which(datafile$length >= minlen),]
```
10. Set filter field
```R
datafile$filter1 <- SetFilterModel1(datafile)
```
11. Split data into train and test
```R
base_X <- datafile[,features]
base_Y <- datafile[,"ml_mu"]
base_filter <- datafile[,"filter1"]
base_weights <- sapply(datafile[,"ml_mu_stderr_adj"], function(x) 1/(x^2))

set.seed(200713)
train.idx <- sample(1:length(base_X[,1]),round(length(base_X[,1])*0.75))
train_X <- base_X[train.idx,]
test_X <- base_X[-(train.idx),]
train_Y <- base_Y[train.idx]
test_Y <- base_Y[-(train.idx)]
train_filter <- base_filter[train.idx]
test_filter <- base_filter[-(train.idx)]
train_weights <- base_weights[train.idx]
test_weights <- base_weights[-(train.idx)]
```
12. Filter training data
```R
train_X <- train_X[which(train_filter == FALSE),]
train_Y <- train_Y[which(train_filter == FALSE)]
train_weights <- train_weights[which(train_filter == FALSE)]
```
13. Build linear regression model
```R
df.lm <- cbind(train_X, train_Y)
model1 <- lm(train_Y ~ uninterrupted_length + period, data = df.lm)
```
14. Apply model to our data
```R
data.sorted$motif.size <- nchar(data.sorted$motif.family)
model.data <- data.frame(data.sorted$length.uninterrupted,
                         data.sorted$motif.size)
colnames(model.data) <- c("uninterrupted_length","period")
data.sorted$ml_mu <- round(predict(model1,model.data),2)
```
## Step 15: Export and compress data for use with shiny app
1. Export data from R as .csv file
```R
write.csv(data.sorted,"./[FILENAME].final.csv")
```
2. Compress file using ```zip``` function in the command line
```bash
zip ./[FILENAME].final.csv > ./[FILENAME].final.csv.zip
```

Congratulations! You now have an annotated set of microsatellites that is compatible with STRATIFY. See [here](README.md) for instructions on how to use this with STRATIFY.
