#' Align two sequences globally using Hirschberg's algorithm
#' 
#' @param X DNAString object representing NT or AA sequence to be aligned.
#' @param Y DNAString object representing NT or AA sequence to be aligned.
#' @param alignment A list of DNAString objects with alignment of input sequences.
#' @param match An integer value of a score for matching bases.
#' @param mismatch An integer value of a score for mismatching bases.
#' @param gap An integer value of a penalty for gap insertion.
#' @returns A list of DNAString objects with alignment of input sequences.
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
library(Biostrings)



NeedlemanWunsch <- function(A, B, match = 1, mismatch = -1, gap = -2) {
  # Create matrices to store scores and traceback
  n <- nchar(A) + 1
  m <- nchar(B) + 1
  score_matrix <- matrix(0, nrow = n, ncol = m)
  traceback_matrix <- matrix(0, nrow = n, ncol = m)
  
  # Initialize the first row and column with gap penalties
  for (i in 2:n) {
    score_matrix[i, 1] <- score_matrix[i - 1, 1] + gap
    traceback_matrix[i, 1] <- "???"
  }
  for (j in 2:m) {
    score_matrix[1, j] <- score_matrix[1, j - 1] + gap
    traceback_matrix[1, j] <- "???"
  }
  
  # Fill in the score matrix and traceback matrix
  for (i in 2:n) {
    for (j in 2:m) {
      charA <- substr(A, i - 1, i - 1)
      charB <- substr(B, j - 1, j - 1)
      
      match_score <- ifelse(charA == charB, match, mismatch)
      
      diagonal <- score_matrix[i - 1, j - 1] + match_score
      up <- score_matrix[i - 1, j] + gap
      left <- score_matrix[i, j - 1] + gap
      
      score_matrix[i, j] <- max(diagonal, up, left)
      
      if (score_matrix[i, j] == diagonal) {
        traceback_matrix[i, j] <- "???"
      } else if (score_matrix[i, j] == up) {
        traceback_matrix[i, j] <- "???"
      } else {
        traceback_matrix[i, j] <- "???"
      }
    }
  }
  
  # Traceback to form the alignment
  i <- n
  j <- m
  aligned_A <- ""
  aligned_B <- ""
  
  while (i > 1 || j > 1) {
    if (i > 1 && j > 1 && traceback_matrix[i, j] == "???") {
      aligned_A <- paste0(substr(A, i - 1, i - 1), aligned_A)
      aligned_B <- paste0(substr(B, j - 1, j - 1), aligned_B)
      i <- i - 1
      j <- j - 1
    } else if (i > 1 && traceback_matrix[i, j] == "???") {
      aligned_A <- paste0(substr(A, i - 1, i - 1), aligned_A)
      aligned_B <- paste0("-", aligned_B)
      i <- i - 1
    } else if (j > 1 && traceback_matrix[i, j] == "???") {
      aligned_A <- paste0("-", aligned_A)
      aligned_B <- paste0(substr(B, j - 1, j - 1), aligned_B)
      j <- j - 1
    }
  }
  
  alignment <- list(aligned_A, aligned_B)
  
  # Return alignment and score matrix
  return(alignment)
}

A <- "GATTACA"
B <- "GCATGCU"

result <- NeedlemanWunsch(A, B)



needleman = function(seq1, seq2, match, mismatch, match, gap){
  
  # Stop conditions
  stopifnot(gap <= 0) # check if penalty negative
  stopifnot(mismatch <= 0)  # check if penalty negative
  stopifnot(match >= 0)  # check if score positive
  
  # Initialize col and rownames for matrices
  len1 = nchar(seq1); len2 = nchar(seq2) # Save number of chars in each sequence
  seq1 = unlist(strsplit(seq1, split = "")) # convert seq to character vector
  seq2 = unlist(strsplit(seq2, split = "")) # convert seq to character vector
  
  # Initialize matrix M (for scores)
  M = matrix(0, nrow = len1 + 1, ncol = len2 + 1) # Initialize matrix
  rownames(M) = c("-", seq1) # assign seq chars to matrix names
  colnames(M) = c("-", seq2) # assign seq chars to matrix names
  M[1, ] = cumsum(c(0, rep(gap, len2))) # Fill 1st row with gap penalites
  M[, 1] = cumsum(c(0, rep(gap, len1))) # Fill 1st col with gap penalites
  
  # Initialize matrix D (for directions)
  D = matrix(0, nrow = len1 + 1, ncol = len2 + 1) # Initialize matrix
  rownames(D) = c("-", seq1) # assign seq chars to matrix names
  colnames(D) = c("-", seq2) # assign seq chars to matrix names
  D[1, ] = rep("hor") # Fill 1st row with "hor" for horizontal moves
  D[, 1] = rep("ver") # Fill 1st col with "ver" for vertical moves
  type = c("dia", "hor", "ver") # Lookup vector
  
  # Compute scores and save moves
  for (i in 2:(len1 + 1)){# for every (initially zero) row
    for (j in 2:(len2 + 1)){# for every (initially zero) col
      hor = M[i, j - 1] + gap # horizontal move = gap for seq1
      ver = M[i - 1, j] + gap # vertical move = gap for seq2 
      dia = ifelse(rownames(M)[i] == colnames(M)[j], # diagonal = ifelse(chars equal, match, mismatch) 
                   M[i - 1, j - 1] + match, 
                   M[i - 1, j - 1] + mismatch)
      M[i, j] = max(dia, hor, ver) # Save current (best) score in M
      D[i, j] = type[which.max(c(dia, hor, ver))] # Save direction of move in D
    }
  } 
  
  # Backtracing
  align1 = c(); align2 = c() # Note: length of final alignments is unknown at this point
  
  while(i > 1 && j > 1){
    
    if(D[i, j] == "dia") {
      align1 = c(rownames(M)[i], align1)
      align2 = c(colnames(M)[j], align2)
      j = j - 1; i = i - 1  # update indices
    } else if (D[i, j] == "ver") {
      align1 = c(rownames(M)[i], align1)
      align2 = c("-", align2) # vertical movement = gap for seq2
      i = i - 1 # update indices
    } else if (D[i, j] == "hor") {
      align1 = c("-", align1) # horizontal movement = gap for seq1
      align2 = c(colnames(M)[j], align2) 
      j = j - 1 # update indices
    } 
    
  }
  
  # Prepare output
  return(list(aligned_seqs = matrix(c(align1, align2), byrow = TRUE, nrow = 2),
              score = M[nrow(M), ncol(M)], score_matrix = M, movement_matrix = D))
  
}

HirschbergTemplate <- function(X, Y, alignment, match, mismatch, gap){
    
    first_alignment_row <- alignment[[1]] # initialize the first row of alignment
    second_alignment_row <- alignment[[2]] # initialize the second row of alignment
  
    if (length(X) == 0) { # length of X is equal to zero
        for (y in Y) { # for each character in Y
            first_alignment_row <- append(first_alignment_row, "-")
            second_alignment_row <- append(second_alignment_row, y)
        }
        alignment <- c(DNAStringSet(first_alignment_row), DNAStringSet(second_alignment_row))
        print(alignment)
    } else if (length(Y) == 0) { # length of Y is equal to zero
        for (x in X) { # for each character in X
            first_alignment_row <- append(first_alignment_row, x)
            second_alignment_row <- append(second_alignment_row, "-")
        }
        alignment <- c(DNAStringSet(first_alignment_row), DNAStringSet(second_alignment_row))
        print(alignment)
    } else if (length(X) == 1 & length(Y) == 1) { # length of X and Y is equal to 1
        first_alignment_row <- # add character from X
        second_alignment_row <- # add character from Y
        alignment <- c(DNAStringSet(first_alignment_row), DNAStringSet(second_alignment_row))
        print(alignment)
    } else {
        xlen <- length(X)
        xmid <- floor(xlen / 2)
        ylen <- length(Y)
        
        matrix1 <- needleman(first_alignment_row[1:xmid], second_alignment_row, match, mismatch, gap)
        matrix2 <- needleman(first_alignment_row[xmid+1:xlen], second_alignment_row, match, mismatch, gap)
        
        # NW score for the first half of X and the whole Y
        first_score <- matrix1[nrow(matrix1), ]
        # NW score for the second half of X and the whole Y (both are reversed)
        second_score <- matrix2[nrow(matrix2), ]

        ymid <- which.max(second_score) - 1 # index of division for Y

        # The first half
        if (ymin == 0) { # index of division for Y is equal to 0
            # call Hirschberg function for the first half of X and for an empty DNAString object
            alignment <- HirschbergTemplate(first_alignment_row[1:xmid], DNAString(""), alignment, match, mismatch, gap)
        } else {
            # call Hirschberg function for the first half of X and for the first part of Y
            alignment <- HirschbergTemplate(first_alignment_row[1:xmid], second_alignment_row[1:ymin], alignment, match, mismatch, gap)
        }
        
        # The second half
        if ((xmid + 1) > xlen) { # X cannot be further divided
            # call Hirschberg function for an empty DNAString object and the second half of Y
            alignment <- HirschbergTemplate(DNAString(""), second_alignment_row[ymin+1:ylen], alignment, match, mismatch, gap)
        } else if ((ymid + 1) > ylen) { # Y cannot be further divided
            # call Hirschberg function for the second half of X and for an empty DNAString object
            alignment <- HirschbergTemplate(first_alignment_row[xmid+1:xlen], DNAString(""), alignment, match, mismatch, gap)
        } else {
            # call Hirschberg function for the second half of X and the second part of Y
            alignment <- HirschbergTemplate(first_alignment_row[xmid+1:xlen], second_alignment_row[ymin+1:ylen],alignment, match, mismatch, gap)
        }
    }
    return(alignment)
}


X = DNAString("TACGAGGCA")
Y = DNAString("ACGGA")
match = 3
mismatch = 1
gap = -3

alignment <- matrix(0, nrow = length(X), ncol = length(Y))


outputH <- HirschbergTemplate(X, Y, alignment, match, mismatch, gap)

