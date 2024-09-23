#' Align two sequences globally using Hirschberg's algorithm
#' 
#' @param X DNAString object representing NT or AA sequence to be aligned.
#' @param Y DNAString object representing NT or AA sequence to be aligned.
#' @param align A list of DNAString objects with alignment of input sequences.
#' @param match An integer value of a score for matching bases.
#' @param mismatch An integer value of a score for mismatching bases.
#' @param gap An integer value of a penalty for gap insertion.
#' @returns A list of DNAString objects with alignment of input sequences.
HirschbergTemplate <- function(X, Y, align, match, mismatch, gap){
    
    first_align_row <- align[[1]] # initialize the first row of alignment
    second_align_row <- align[[2]] # initialize the second row of alignment
  
  
    if # length of X is equal to zero
    {
        for # for each character in Y
        {
            first_align_row <- # add gap
            second_align_row <- # add character from Y
        }
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else if # length of Y is equal to zero
    {
        for # for each character in X
        {
            first_align_row <- # add character from X
            second_align_row <- # add gap
        }
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else if # length of X and Y is equal to 1
    {
        first_align_row <- # add character from X
        second_align_row <- # add character from Y
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else
    {
        xlen <- # length of X
        xmid <- # half of the length of X
        ylen <- # length of Y
        
        left_score <- # NW score for the first half of X and the whole Y
        right_score <- # NW score for the second half of X and the whole Y (both are reversed)
        ymid <- # index of division for Y

        # The first half
        if # index of division for Y is equal to 0
        {
            align <- # call Hirschberg function for the first half of X and for an empty DNAString object
        }
        else
        {
            align <- # call Hirschberg function for the first half of X and for the first part of Y
        }
        
        # The second half
        if ((xmid + 1) > xlen) # X cannot be further divided
        {
            align <- # call Hirschberg function for an empty DNAString object and the second half of Y
        }
        else if ((ymid + 1) > ylen) # Y cannot be further divided
        {
            align <- # call Hirschberg function for the second half of X and for an empty DNAString object
        }
        else 
        {
            align <- # call hirschberg function for the second half of X and the second part of Y
        }
    }
    return(align)
}