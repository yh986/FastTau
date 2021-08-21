library(Rcpp)
sourceCpp("IntegerTau.cpp")

#Functions to prepare real vectors into integer vectors for IntegerTau.cpp's TauSVIV function

#SV() sorts a vector X and gives information about the sorting;
#such that the first value in the output gives the index of the lowest value in X;
#the second value in the output gives the index of the next highest value in X, and so on.
#Tied values are represented by negative adjacent values in the output;
#and NaN (missing data) are reported as 0s at the end of the output.
SV <- function(x) {
   Rorder <- order(x, decreasing = F, na.last = NA)
   sort <- integer(length(x))
   adjacentdiff <- diff(x[Rorder])
   ties <- adjacentdiff == 0
   dties <- diff(ties)
   dties_complete <- c(ties[1], dties != 0, ties[length(ties)])
   Rorder[dties_complete] <- -Rorder[dties_complete]
   sort[1:length(Rorder)] <- Rorder
   sort[is.na(sort)] <- 0 #This is to deal with the all-Na-input case
   sort
}

#IV() makes an integer representation of a vector Y:
#In Y, NaNs become 0, the lowest finite value becomes '1', the next higher value becomes '2', etc;
#Finally, the very first value in Y codes for the position of the largest value.
#Requires a 'sorted' vector representation, sv, to know the rankings of each index of the original vector!
IV <- function(sv) {
   integerized <- integer(length(sv) + 1) #First value indicates index of largest
   j <- 1; i <- 1; streak <- FALSE
   while (j <= length(sv) && sv[j] != 0) {
     if (sv[j] < 0) {
       integerized[1 - sv[j]] <- i
       streak <- !streak
       if (!streak) {
         i <- i + 1
       }
     }
     else {
       integerized[1 + sv[j]] <- i
       if (!streak) {
         i <- i + 1
       }
     }
     j <- j + 1
   }
   #j is now one over, unless sv is all 0s
   if (j == 1) {
     integerized[1] <- 1
   }
   else
   {
     integerized[1] <- abs(sv[j - 1]) + 1 
   }
   integerized
 }

#This computes the Tau correlation between vectors or matrices A and B!
#If given matrices, then it will do all possible correlations between
#the columns of the inputs, and output a matrix.

Tau <- function(A, B = NULL) {
  if (!is.atomic(A) || !is.atomic(B)) {
    stop("Need to provide atomic vectors or matrices")
  }
  
  if (is.null(dim(A))) { #vector, not matrix
    A <- matrix(A)
  }
  
  if (!is.null(B) && is.null(dim(B))) #Non-null vector
  {
    B <- matrix(B)
  }
  
  if (!is.null(B) && nrow(A) != nrow(B)) {
    stop("Unequal lengths")
  }
  
  svA <- apply(A, 2, SV)
  
  if (!is.null(B)) {
    svB <- apply(B, 2, SV)
  }
  else
  {
    svB <- svA
  }
  
  ivB <- apply(svB, 2, IV)
  
  corout <- matrix(nrow = ncol(svA), ncol = ncol(svB))
  
  for (i in 1:ncol(svA)) {
    corcol <- svA[, i]
    for(j in 1:ncol(ivB)){
      corout[i, j] <- TauSVIV(corcol, ivB[, j])
    }
  }
  
  if (nrow(corout) == 1) {
    corout <- as.vector(corout)
  }
  corout
}
