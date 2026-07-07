#
# added by procedure shortening function names
#
                                                                 # nolint start

lav_matrix_duplication_ginv_pre_post <- function(A = matrix(0, 0, 0)) {
  lav_deprecated("lav_mat_dup_ginv_pre_post", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_dup_ginv_pre_post)
  eval(sc, parent.frame())
}
lav_matrix_orthogonal_complement <- function(A = matrix(0, 0, 0)) {
  lav_deprecated("lav_mat_ortho_complement", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_ortho_complement)
  eval(sc, parent.frame())
}
lav_func_gradient_complex <- function(
                                      func,
                                      x,
                                      h = .Machine$double.eps,
                                      ...,
                                      fallback.simple = TRUE) {
  lav_deprecated("lav_func_grad_complex", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_func_grad_complex)
  eval(sc, parent.frame())
}
lav_matrix_duplication_ginv_post <- function(A = matrix(0, 0, 0)) {
  lav_deprecated("lav_mat_dup_ginv_post", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_dup_ginv_post)
  eval(sc, parent.frame())
}
lav_func_gradient_simple <- function(
                                     func,
                                     x,
                                     h = sqrt(.Machine$double.eps),
                                     ...) {
  lav_deprecated("lav_func_grad_simple", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_func_grad_simple)
  eval(sc, parent.frame())
}
lav_matrix_vech_row_idx <- function(n = 1L, diagonal = TRUE) {
  lav_deprecated("lav_mat_vech_row_idx", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vech_row_idx)
  eval(sc, parent.frame())
}
lav_matrix_vech_col_idx <- function(n = 1L, diagonal = TRUE) {
  lav_deprecated("lav_mat_vech_col_idx", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vech_col_idx)
  eval(sc, parent.frame())
}
lav_matrix_antidiag_idx <- function(n = 1L) {
  lav_deprecated("lav_mat_antidiag_idx", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_antidiag_idx)
  eval(sc, parent.frame())
}
lav_matrix_duplication_pre_post <- function(A = matrix(0, 0, 0)) {
  lav_deprecated("lav_mat_dup_pre_post", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_dup_pre_post)
  eval(sc, parent.frame())
}
lav_matrix_duplication_ginv_pre <- function(A = matrix(0, 0, 0)) {
  lav_deprecated("lav_mat_dup_ginv_pre", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_dup_ginv_pre)
  eval(sc, parent.frame())
}
lav_matrix_commutation_pre_post <- function(A = matrix(0, 0, 0)) {
  lav_deprecated("lav_mat_com_pre_post", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_com_pre_post)
  eval(sc, parent.frame())
}
lav_partable_unrestricted <- function(
                                      lavobject = NULL,
                                      lavdata = NULL,
                                      lavpta = NULL,
                                      lavoptions = NULL,
                                      lavsamplestats = NULL,
                                      lavh1 = NULL,
                                      sample.cov = NULL,
                                      sample.mean = NULL,
                                      sample.slopes = NULL,
                                      sample.th = NULL,
                                      sample.th.idx = NULL,
                                      sample.cov.x = NULL,
                                      sample.mean.x = NULL) {
  lav_deprecated("lav_pt_unrestricted", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_pt_unrestricted)
  eval(sc, parent.frame())
}
lav_partable_independence <- function(
                                      lavobject = NULL,
                                      lavdata = NULL,
                                      lavpta = NULL,
                                      lavoptions = NULL,
                                      lavsamplestats = NULL,
                                      lavh1 = NULL,
                                      sample.cov = NULL,
                                      sample.mean = NULL,
                                      sample.slopes = NULL,
                                      sample.th = NULL,
                                      sample.th.idx = NULL,
                                      sample.cov.x = NULL,
                                      sample.mean.x = NULL) {
  lav_deprecated("lav_pt_independence", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_pt_independence)
  eval(sc, parent.frame())
}
lav_matrix_vechru_idx <- function(n = 1L, diagonal = TRUE) {
  lav_deprecated("lav_mat_vechru_idx", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vechru_idx)
  eval(sc, parent.frame())
}
lav_matrix_upper2full <- function(x, diagonal = TRUE) {
  lav_deprecated("lav_mat_upper2full", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_upper2full)
  eval(sc, parent.frame())
}
lav_matrix_lower2full <- function(x, diagonal = TRUE) {
  lav_deprecated("lav_mat_lower2full", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_lower2full)
  eval(sc, parent.frame())
}
lav_matrix_commutation_mn_pre <- function(
                                          A,
                                          m = 1L,
                                          n = 1L) {
  lav_deprecated("lav_mat_com_mn_pre", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_com_mn_pre)
  eval(sc, parent.frame())
}
lav_samplestats_from_data <- function(
                                      lavdata = NULL,
                                      lavoptions = NULL,
                                      WLS.V = NULL,
                                      NACOV = NULL) {
  lav_deprecated("lav_samp_from_data", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_samp_from_data)
  eval(sc, parent.frame())
}
lav_matrix_vechru_reverse <- function(x, diagonal = TRUE) {
  lav_deprecated("lav_mat_vechru_rev", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vechru_rev)
  eval(sc, parent.frame())
}
lav_matrix_vechr_idx <- function(n = 1L, diagonal = TRUE) {
  lav_deprecated("lav_mat_vechr_idx", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vechr_idx)
  eval(sc, parent.frame())
}
lav_matrix_vechu_idx <- function(n = 1L, diagonal = TRUE) {
  lav_deprecated("lav_mat_vechu_idx", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vechu_idx)
  eval(sc, parent.frame())
}
lav_matrix_diagh_idx <- function(n = 1L) {
  lav_deprecated("lav_mat_diagh_idx", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_diagh_idx)
  eval(sc, parent.frame())
}
lav_partable_attributes <- function(partable, pta = NULL) {
  lav_deprecated("lav_pt_attributes", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_pt_attributes)
  eval(sc, parent.frame())
}
lav_matrix_vechr_reverse <- function(x, diagonal = TRUE) {
  lav_deprecated("lav_mat_vechr_rev", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vechr_rev)
  eval(sc, parent.frame())
}
lav_matrix_vechu_reverse <- function(x, diagonal = TRUE) {
  lav_deprecated("lav_mat_vechu_rev", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vechu_rev)
  eval(sc, parent.frame())
}
lav_matrix_vech_idx <- function(n = 1L, diagonal = TRUE) {
  lav_deprecated("lav_mat_vech_idx", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vech_idx)
  eval(sc, parent.frame())
}
lav_matrix_diag_idx <- function(n = 1L) {
  lav_deprecated("lav_mat_diag_idx", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_diag_idx)
  eval(sc, parent.frame())
}
lav_matrix_duplication_post <- function(A = matrix(0, 0, 0)) {
  lav_deprecated("lav_mat_dup_post", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_dup_post)
  eval(sc, parent.frame())
}
lav_matrix_duplication_ginv <- function(n = 1L) {
  lav_deprecated("lav_mat_dup_ginv", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_dup_ginv)
  eval(sc, parent.frame())
}
lav_matrix_commutation_post <- function(A = matrix(0, 0, 0)) {
  lav_deprecated("lav_mat_com_post", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_com_post)
  eval(sc, parent.frame())
}
lav_matrix_symmetric_sqrt <- function(S = matrix(0, 0, 0)) {
  lav_deprecated("lav_mat_sym_sqrt", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_sym_sqrt)
  eval(sc, parent.frame())
}
lav_matrix_vech_reverse <- function(x, diagonal = TRUE) {
  lav_deprecated("lav_mat_vech_rev", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vech_rev)
  eval(sc, parent.frame())
}
lav_matrix_duplication_pre <- function(A = matrix(0, 0, 0)) {
  lav_deprecated("lav_mat_dup_pre", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_dup_pre)
  eval(sc, parent.frame())
}
lav_matrix_commutation_pre <- function(A = matrix(0, 0, 0)) {
  lav_deprecated("lav_mat_com_pre", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_com_pre)
  eval(sc, parent.frame())
}
lav_partable_complete <- function(partable = NULL, start = TRUE) {
  lav_deprecated("lav_pt_complete", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_pt_complete)
  eval(sc, parent.frame())
}
lav_matrix_vechru <- function(S, diagonal = TRUE) {
  lav_deprecated("lav_mat_vechru", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vechru)
  eval(sc, parent.frame())
}
lav_partable_constraints_def <- function(
                                         partable,
                                         con = NULL,
                                         debug = FALSE,
                                         txtOnly = FALSE,
                                         warn = TRUE) {
  lav_deprecated("lav_pt_con_def", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_pt_con_def)
  eval(sc, parent.frame())
}
lav_partable_constraints_ceq <- function(
                                         partable,
                                         con = NULL,
                                         debug = FALSE,
                                         txtOnly = FALSE) {
  lav_deprecated("lav_pt_con_ceq", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_pt_con_ceq)
  eval(sc, parent.frame())
}
lav_partable_constraints_ciq <- function(
                                         partable,
                                         con = NULL,
                                         debug = FALSE,
                                         txtOnly = FALSE) {
  lav_deprecated("lav_pt_con_ciq", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_pt_con_ciq)
  eval(sc, parent.frame())
}
lav_partable_from_lm <- function(
                                 object,
                                 est = FALSE,
                                 label = FALSE,
                                 as.data.frame. = FALSE) {
  lav_deprecated("lav_pt_from_lm", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_pt_from_lm)
  eval(sc, parent.frame())
}
lav_constraints_parse <- function(
                                  partable = NULL,
                                  constraints = NULL,
                                  theta = NULL,
                                  debug = FALSE) {
  lav_deprecated("lav_con_parse", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_con_parse)
  eval(sc, parent.frame())
}
lav_matrix_vechr <- function(S, diagonal = TRUE) {
  lav_deprecated("lav_mat_vechr", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vechr)
  eval(sc, parent.frame())
}
lav_matrix_vechu <- function(S, diagonal = TRUE) {
  lav_deprecated("lav_mat_vechu", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vechu)
  eval(sc, parent.frame())
}
lav_matrix_bdiag <- function(...) {
  lav_deprecated("lav_mat_bdiag", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_bdiag)
  eval(sc, parent.frame())
}
lav_matrix_trace <- function(..., check = TRUE) {
  lav_deprecated("lav_mat_trace", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_trace)
  eval(sc, parent.frame())
}
lav_partable_labels <- function(
                                partable,
                                blocks = c("group", "level"),
                                group.equal = "",
                                group.partial = "",
                                type = "user") {
  lav_deprecated("lav_pt_labels", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_pt_labels)
  eval(sc, parent.frame())
}
# used in bain
getParameterLabels <- lav_alias("lav_partable_labels")

lav_matrix_vecr <- function(A) {
  lav_deprecated("lav_mat_vecr", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vecr)
  eval(sc, parent.frame())
}
lav_matrix_vech <- function(S, diagonal = TRUE) {
  lav_deprecated("lav_mat_vech", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vech)
  eval(sc, parent.frame())
}
lav_partable_merge <- function(
                               pt1 = NULL,
                               pt2 = NULL,
                               remove.duplicated = FALSE,
                               fromLast = FALSE,
                               warn = TRUE) {
  lav_deprecated("lav_pt_merge", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_pt_merge)
  eval(sc, parent.frame())
}
lav_matrix_vec <- function(A) {
  lav_deprecated("lav_mat_vec", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_vec)
  eval(sc, parent.frame())
}
lav_matrix_duplication <- function(n = 1L) {
  lav_deprecated("lav_mat_dup", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_dup)
  eval(sc, parent.frame())
}
lav_matrix_commutation <- function(m = 1L, n = 1L) {
  lav_deprecated("lav_mat_com", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_com)
  eval(sc, parent.frame())
}
lav_matrix_cov <- function(Y, Mu = NULL) {
  lav_deprecated("lav_mat_cov", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_mat_cov)
  eval(sc, parent.frame())
}
lav_partable_ndat <- function(partable) {
  lav_deprecated("lav_pt_ndat", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_pt_ndat)
  eval(sc, parent.frame())
}
lav_partable_npar <- function(partable) {
  lav_deprecated("lav_pt_npar", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_pt_npar)
  eval(sc, parent.frame())
}
lav_partable_add <- function(partable = NULL, add = list()) {
  lav_deprecated("lav_pt_add", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_pt_add)
  eval(sc, parent.frame())
}
lav_partable_df <- function(partable) {
  lav_deprecated("lav_pt_df", times = 0L) # --> for now no warning
  sc <- sys.call()
  names(sc) <- lav_snake_case(names(sc))
  sc[[1L]] <- quote(lavaan::lav_pt_df)
  eval(sc, parent.frame())
}                                                             # nolint end
