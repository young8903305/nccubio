# read PAM1 from data
pam1<-read.table("pam1.txt",header = TRUE)    #read file without seperate

# check PAM1 data
dim(pam1)    # show/set the dimension
str(pam1)    # show the structure


# construct PAM250 from PAM1
multiple  <- function(matrix, n){    # function of matrix multiply
  if(n == 1)
    return (matrix)
  else
    return (matrix %*% multiple(matrix, n-1))
}

p <- as.matrix( pam1 ) / 10000    # change pam1 data.frame into matrix
p250 <- round( multiple(p, 250) * 100, digits = 0 )    # result * 100, then round to digit 0

# output PAM250 as a file
write.table(p250, file = "pam250.txt", sep = "\t", quote = FALSE, append = FALSE, na = "NA")    #seperate with \t, no quote, can overwrite the same file

