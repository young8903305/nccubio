######################################
# the reference code of program2 
######################################

######################################
# initial
######################################
source("https://bioconductor.org/biocLite.R")
biocLite()
library("Biostrings",verbose=F,quietly=T)


# read parameters
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Rscript hw2_105753024_v2.R --input test.fasta --score PAM250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta", call.=FALSE)
}

gap_open <- 0
gap_extend <- 0

# parse parameters
i <- 1 
while(i < length(args))
{
  if(args[i] == "--input"){
    input_file <- args[i+1]
    i <- i + 1
  }else if(args[i] == "--score"){
    score_file <- args[i+1]
    i <- i + 1
  }else if(args[i] == "--aln"){
    aln_mode <- args[i+1]
    i <- i + 1
  }else if(args[i] == "--gap_open"){
    gap_open <- as.integer(args[i+1]) 
    i <- i + 1
  }else if(args[i] == "--gap_extend"){
    gap_extend <- as.integer(args[i+1])
    i <- i + 1    
  }else if(args[i] == "--output"){
    output_file <- args[i+1]
    i <- i + 1
  }else{
    stop(paste("Unknown flag", args[i]), call.=FALSE)
  }
  i <- i + 1
}

print("PARAMETERS")
print(paste("input file         :", input_file))
print(paste("output file        :", output_file))
print(paste("score file         :", score_file))
print(paste("aln_mode           :", aln_mode))
print(paste("gap open penalty   :", gap_open))
print(paste("gap extend penalty :", gap_extend))

gap_penalty <- as.integer(gap_extend)

######################################
# main
######################################
# read fasta file
ff <- readAAStringSet(input_file)
seq_name = names(ff)
sequence = paste(ff)

# aln length
aln_length1 <- nchar(sequence[1])
aln_length2 <- nchar(sequence[2])

# read score file
s_m <- read.table(score_file)
score_matrix <- as.matrix(s_m)

F <- matrix(NA, aln_length1 + 1, aln_length2 + 1) # row*col, alignment matrix

if(aln_mode == "global"){ ######### global alignment #########
  
  # initialize first row and first col with gap penalty decreasing
  F[,1] <- seq(0, aln_length1 * gap_penalty, gap_penalty)
  F[1,] <- seq(0, aln_length2 * gap_penalty, gap_penalty)
  
  for(i in 2:(nrow(F)) ){
    for(j in 2:(ncol(F)) ){
      match <- F[i-1,j-1] + score_matrix[substring(sequence[1], i-1, i-1), substring(sequence[2], j-1, j-1)]
      gapUp <- F[i-1,j] + gap_penalty
      gapLeft<- F[i,j-1] + gap_penalty
      F[i,j] <- max( c(match, gapUp, gapLeft) )
    }
  }

  # null string
  rlt1 <- ""
  rlt2 <- ""
  
  i <- nrow(F)
  j <- ncol(F)
  # trace back
  while(i > 1 && j > 1){
    if(F[i,j] - score_matrix[substring(sequence[1], i-1, i-1), substring(sequence[2], j-1, j-1)] == F[i-1,j-1]){
      rlt1 <- paste0(substring(sequence[1], i-1, i-1), rlt1)
      rlt2 <- paste0(substring(sequence[2], j-1, j-1), rlt2)
      i <- i - 1
      j <- j - 1
    }else if(F[i,j] - gap_penalty == F[i,j-1]){
      rlt1 <- paste0("-", rlt1)
      rlt2 <- paste0(substring(sequence[2], j-1, j-1), rlt2)
      j <- j - 1
    }else if(F[i,j] - gap_penalty == F[i-1,j]){
      rlt1 <- paste0(substring(sequence[1], i-1, i-1), rlt1)
      rlt2 <- paste0("-", rlt2)
      i <- i - 1
    }
  }
  # tail of sequence1 and sequence2
  if(j > 1){
    while(j > 1){
      rlt1 = paste0("-", rlt1)
      rlt2 = paste0(substring(sequence[2], j-1, j-1), rlt2)
      j <- j - 1
    }
  }else if(i > 1){
    while(i > 1){
      rlt1 = paste0(substring(sequence[1], i-1, i-1), rlt1)
      rlt2 = paste0("-", rlt2)
      i <- i - 1
    }
  }
}else if(aln_mode == "local"){ ######### local alignment #########
  #initial matrix with int 0
  for(i in 0 : nrow(F)){
    for(j in 0 : ncol(F)){
      F[i, j] <- 0
    }
  }
  
  max <- 0
  for(i in 2:(nrow(F)) ){
    for(j in 2:(ncol(F)) ){
      match <- F[i-1,j-1]+score_matrix[substring(sequence[1], i-1, i-1), substring(sequence[2], j-1, j-1)]
      gapUp <- F[i-1,j] + gap_penalty
      gapLeft<- F[i,j-1] + gap_penalty
      F[i,j] <- max(c(match, gapUp, gapLeft, 0))
      if(F[i, j] > max){
        index_i <- i
        index_j <- j
        max <- F[i, j]
      }
    }
  }
  
  # null string
  rlt1 <- ""
  rlt2 <- ""
  
  #trace back
  i <- index_i
  j <- index_j
  while(max != 0){
    if(F[i,j] - score_matrix[substring(sequence[1], i-1, i-1), substring(sequence[2], j-1, j-1)] == F[i-1,j-1]){
      rlt1 <- paste0(substring(sequence[1], i-1, i-1), rlt1)
      rlt2 <- paste0(substring(sequence[2], j-1, j-1), rlt2)
      i <- i - 1
      j <- j - 1
    }else if(F[i,j] - gap_penalty == F[i,j-1]){
      rlt1 <- paste0("-", rlt1)
      rlt2 <- paste0(substring(sequence[2], j-1, j-1), rlt2)
      j <- j - 1
      
    }else if(F[i,j] - gap_penalty == F[i-1,j]){
      rlt1 <- paste0(substring(sequence[1], i-1, i-1), rlt1)
      rlt2 <- paste0("-", rlt2)
      i <- i - 1
    }
    max <- F[i, j]
  }
}

print(rlt1)
print(rlt2)

ff[1] <- rlt1
ff[2] <- rlt2
print(ff)

# output
writeXStringSet(ff, output_file)