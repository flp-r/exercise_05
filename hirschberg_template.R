#' Align two sequences globally using Hirschberg's algorithm
#' 
#' @param X DNAString object representing NT or AA sequence to be aligned.
#' @param Y DNAString object representing NT or AA sequence to be aligned.
#' @param alignment A list of DNAString objects with alignment of input sequences.
#' @param match An integer value of a score for matching bases.
#' @param mismatch An integer value of a score for mismatching bases.
#' @param gap An integer value of a penalty for gap insertion.
#' @returns A list of DNAString objects with alignment of input sequences.
HirschbergTemplate <- function(X, Y, alignment, match, mismatch, gap){
    
    first_alignment_row <- alignment[[1]] # initialize the first row of alignment
    second_alignment_row <- alignment[[2]] # initialize the second row of alignment
  
    if () { # length of X is equal to zero
        for () { # for each character in Y
            first_alignment_row <- # add gap
            second_alignment_row <- # add character from Y
        }
        alignment <- c(DNAStringSet(first_alignment_row), DNAStringSet(second_alignment_row))
        print(alignment)
    } else if () { # length of Y is equal to zero
        for () { # for each character in X
            first_alignment_row <- # add character from X
            second_alignment_row <- # add gap
        }
        alignment <- c(DNAStringSet(first_alignment_row), DNAStringSet(second_alignment_row))
        print(alignment)
    } else if () { # length of X and Y is equal to 1
        first_alignment_row <- # add character from X
        second_alignment_row <- # add character from Y
        alignment <- c(DNAStringSet(first_alignment_row), DNAStringSet(second_alignment_row))
        print(alignment)
    } else {
        xlen <- # length of X
        xmid <- # half of the length of X
        ylen <- # length of Y

        # NW score for the first half of X and the whole Y
        first_score <-
        # NW score for the second half of X and the whole Y (both are reversed)
        second_score <-

        ymid <- # index of division for Y

        # The first half
        if () { # index of division for Y is equal to 0
            # call Hirschberg function for the first half of X and for an empty DNAString object
            alignment <-
        } else {
            # call Hirschberg function for the first half of X and for the first part of Y
            alignment <-
        }
        
        # The second half
        if ((xmid + 1) > xlen) { # X cannot be further divided
            # call Hirschberg function for an empty DNAString object and the second half of Y
            alignment <-
        } else if ((ymid + 1) > ylen) { # Y cannot be further divided
            # call Hirschberg function for the second half of X and for an empty DNAString object
            alignment <-
        } else {
            # call Hirschberg function for the second half of X and the second part of Y
            alignment <-
        }
    }
    return(alignment)
}
