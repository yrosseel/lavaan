# deprecated (or renamed) lavaan functions
# created 29 March 2015, pre 0.5-18 by YR

vech <- function(S, diagonal = TRUE) {
    .Deprecated("lav_matrix_vech", package = "lavaan")
    lav_matrix_vech(S = S, diagonal = diagonal)
}

vechr <- function(S, diagonal = TRUE) {
    .Deprecated("lav_matrix_vechr", package = "lavaan")
    lav_matrix_vechr(S = S, diagonal = diagonal)
}

vechu <- function(S, diagonal = TRUE) {
    .Deprecated("lav_matrix_vechu", package = "lavaan")
    lav_matrix_vechu(S = S, diagonal = diagonal)
}

vechru <- function(S, diagonal = TRUE) {
    .Deprecated("lav_matrix_vechru", package = "lavaan")
    lav_matrix_vechru(S = S, diagonal = diagonal)
}


vech.reverse <- function(x, diagonal = TRUE) {
    .Deprecated("lav_matrix_vech_reverse",  package = "lavaan")
    lav_matrix_vech_reverse(x = x, diagonal = diagonal)
}

vechru.reverse <- function(x, diagonal = TRUE) {
    .Deprecated("lav_matrix_vechru_reverse",  package = "lavaan")
    lav_matrix_vechru_reverse(x = x, diagonal = diagonal)
}

vechr.reverse <- function(x, diagonal = TRUE) {
    .Deprecated("lav_matrix_vechr_reverse",  package = "lavaan")
    lav_matrix_vechr_reverse(x = x, diagonal = diagonal)
}

vechu.reverse <- function(x, diagonal = TRUE) {
    .Deprecated("lav_matrix_vechu_reverse",  package = "lavaan")
    lav_matrix_vechu_reverse(x = x, diagonal = diagonal)
}

lower2full <- function(x, diagonal = TRUE) {
      .Deprecated("lav_matrix_lower2full",  package = "lavaan")
    lav_matrix_lower2full(x = x, diagonal = diagonal)
}

upper2full <- function(x, diagonal = TRUE) {
    .Deprecated("lav_matrix_upper2full",  package = "lavaan")
    lav_matrix_upper2full(x = x, diagonal = diagonal)
}

duplicationMatrix <- function(n = 1L) {
    .Deprecated("lav_matrix_duplication",  package = "lavaan")
    lav_matrix_duplication(n = n)
}

commutationMatrix <- function(m = 1L, n = 1L) {
    .Deprecated("lav_matrix_commutation",  package = "lavaan")
    lav_matrix_commutation(m = m, n = n)
}

sqrtSymmetricMatrix <- function(S) {
    .Deprecated("lav_matrix_symmetric_sqrt", package = "lavaan")
    lav_matrix_symmetric_sqrt(S = S)
}
