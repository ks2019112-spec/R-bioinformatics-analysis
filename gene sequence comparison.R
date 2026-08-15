#gene sequence comparison

install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("msa")
BiocManager::install("DECIPHER")
BiocManager::install("ape")
BiocManager::install("seqinr")
library(Biostrings)
library(msa)
library(DECIPHER)
library(ape)
library(seqinr)
setwd("C:/Users/pc/Desktop/R_gene_project")
getwd()
list.files()
seq1 <- DNAString("ATGCGTACGTAGCTAGCTAG")
seq2 <- DNAString("ATGCGTACGTAGATAGCTAG")
seq3 <- DNAString("ATGCGTTCGTAGCTAGCTAA")
seq1
seq2
seq3
list.files()
dna<-readDNAStringSet("gene_sequences.fasta")
print(dna)
names(dna)
dna[[1]]
data<-read.csv("gene_sequences.csv")
print(data)
dna<-DNAStringset(data$Sequence)
dna
alignment <- Biostrings::pairwiseAlignment(dna[[1]], dna[[2]])
BiocManager::install("pwalign")
library(pwalign)

alignment <- pairwiseAlignment(dna[[1]], dna[[2]])
alignment
score(alignment)
aln<-msa(dna)
print(aln)
msapPettyPrint(aln,output="pdf",file="Alignment.pdf")
msaPrettyPrint(aln, output="asis")
#DistanceMatrix
# Distance Matrix
dist <- DistanceMatrix(dna)
dist

install.packages("ape")
library(ape)
tree<-nj(dist)
tree
plot(tree)
alphabetFrequency(dna)
letterFrequency(dna, c("A","T","G","C"))
letterFrequency(dna, "GC", as.prob=TRUE)
consensusString(dna)
writeXStringSet(dna, "output.fasta")
write.csv(as.matrix(dist), "distance_matrix.csv")
list.files()
