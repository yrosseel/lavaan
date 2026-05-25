\name{lav_long_names}
\alias{lav_matrix_duplication_ginv_pre_post}
\alias{lav_matrix_orthogonal_complement}
\alias{lav_func_gradient_complex}
\alias{lav_matrix_duplication_ginv_post}
\alias{lav_func_gradient_simple}
\alias{lav_matrix_vech_row_idx}
\alias{lav_matrix_vech_col_idx}
\alias{lav_matrix_antidiag_idx}
\alias{lav_matrix_duplication_pre_post}
\alias{lav_matrix_duplication_ginv_pre}
\alias{lav_matrix_commutation_pre_post}
\alias{lav_partable_unrestricted}
\alias{lav_partable_independence}
\alias{lav_matrix_vechru_idx}
\alias{lav_matrix_upper2full}
\alias{lav_matrix_lower2full}
\alias{lav_matrix_commutation_mn_pre}
\alias{lav_samplestats_from_data}
\alias{lav_matrix_vechru_reverse}
\alias{lav_matrix_vechr_idx}
\alias{lav_matrix_vechu_idx}
\alias{lav_matrix_diagh_idx}
\alias{lav_partable_attributes}
\alias{lav_matrix_vechr_reverse}
\alias{lav_matrix_vechu_reverse}
\alias{lav_matrix_vech_idx}
\alias{lav_matrix_diag_idx}
\alias{lav_matrix_duplication_post}
\alias{lav_matrix_duplication_ginv}
\alias{lav_matrix_commutation_post}
\alias{lav_matrix_symmetric_sqrt}
\alias{lav_matrix_vech_reverse}
\alias{lav_matrix_duplication_pre}
\alias{lav_matrix_commutation_pre}
\alias{lav_partable_complete}
\alias{lav_matrix_vechru}
\alias{lav_partable_constraints_def}
\alias{lav_partable_constraints_ceq}
\alias{lav_partable_constraints_ciq}
\alias{lav_partable_from_lm}
\alias{lav_constraints_parse}
\alias{lav_matrix_vechr}
\alias{lav_matrix_vechu}
\alias{lav_matrix_bdiag}
\alias{lav_matrix_trace}
\alias{lav_partable_labels}
\alias{lav_matrix_print}
\alias{lav_matrix_vecr}
\alias{lav_matrix_vech}
\alias{lav_partable_merge}
\alias{lav_model_partable}
\alias{lav_matrix_vec}
\alias{lav_matrix_duplication}
\alias{lav_matrix_commutation}
\alias{lav_matrix_cov}
\alias{lav_partable_ndat}
\alias{lav_partable_npar}
\alias{lav_partable_add}
\alias{lav_partable_df}
\alias{lav_scores}

\title{deprecated functions with long names}
\description{Functions which are deprecated because of adaption of shorter names in May, 2026}
\usage{
lav_matrix_duplication_ginv_pre_post(A = matrix(0, 0, 0))
lav_matrix_orthogonal_complement(A = matrix(0, 0, 0))
lav_func_gradient_complex(
                          func,
                          x,
                          h = .Machine$double.eps,
                          ...,
                          fallback.simple = TRUE)
lav_matrix_duplication_ginv_post(A = matrix(0, 0, 0))
lav_func_gradient_simple(
                         func,
                         x,
                         h = sqrt(.Machine$double.eps),
                         ...)
lav_matrix_vech_row_idx(n = 1L, diagonal = TRUE)
lav_matrix_vech_col_idx(n = 1L, diagonal = TRUE)
lav_matrix_antidiag_idx(n = 1L)
lav_matrix_duplication_pre_post(A = matrix(0, 0, 0))
lav_matrix_duplication_ginv_pre(A = matrix(0, 0, 0))
lav_matrix_commutation_pre_post(A = matrix(0, 0, 0))
lav_partable_unrestricted(
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
                          sample.mean.x = NULL)
lav_partable_independence(
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
                          sample.mean.x = NULL)
lav_matrix_vechru_idx(n = 1L, diagonal = TRUE)
lav_matrix_upper2full(x, diagonal = TRUE)
lav_matrix_lower2full(x, diagonal = TRUE)
lav_matrix_commutation_mn_pre(
                              A,
                              m = 1L,
                              n = 1L)
lav_samplestats_from_data(
                          lavdata = NULL,
                          lavoptions = NULL,
                          WLS.V = NULL,
                          NACOV = NULL)
lav_matrix_vechru_reverse(x, diagonal = TRUE)
lav_matrix_vechr_idx(n = 1L, diagonal = TRUE)
lav_matrix_vechu_idx(n = 1L, diagonal = TRUE)
lav_matrix_diagh_idx(n = 1L)
lav_partable_attributes(partable, pta = NULL)
lav_matrix_vechr_reverse(x, diagonal = TRUE)
lav_matrix_vechu_reverse(x, diagonal = TRUE)
lav_matrix_vech_idx(n = 1L, diagonal = TRUE)
lav_matrix_diag_idx(n = 1L)
lav_matrix_duplication_post(A = matrix(0, 0, 0))
lav_matrix_duplication_ginv(n = 1L)
lav_matrix_commutation_post(A = matrix(0, 0, 0))
lav_matrix_symmetric_sqrt(S = matrix(0, 0, 0))
lav_matrix_vech_reverse(x, diagonal = TRUE)
lav_matrix_duplication_pre(A = matrix(0, 0, 0))
lav_matrix_commutation_pre(A = matrix(0, 0, 0))
lav_partable_complete(partable = NULL, start = TRUE)
lav_matrix_vechru(S, diagonal = TRUE)
lav_partable_constraints_def(
                             partable,
                             con = NULL,
                             debug = FALSE,
                             txtOnly = FALSE,
                             warn = TRUE)
lav_partable_constraints_ceq(
                             partable,
                             con = NULL,
                             debug = FALSE,
                             txtOnly = FALSE)
lav_partable_constraints_ciq(
                             partable,
                             con = NULL,
                             debug = FALSE,
                             txtOnly = FALSE)
lav_partable_from_lm(
                     object,
                     est = FALSE,
                     label = FALSE,
                     as.data.frame. = FALSE)
lav_constraints_parse(
                      partable = NULL,
                      constraints = NULL,
                      theta = NULL,
                      debug = FALSE)
lav_matrix_vechr(S, diagonal = TRUE)
lav_matrix_vechu(S, diagonal = TRUE)
lav_matrix_bdiag(...)
lav_matrix_trace(..., check = TRUE)
lav_partable_labels(
                    partable,
                    blocks = c("group", "level"),
                    group.equal = "",
                    group.partial = "",
                    type = "user")
lav_matrix_vecr(A)
lav_matrix_vech(S, diagonal = TRUE)
lav_partable_merge(
                   pt1 = NULL,
                   pt2 = NULL,
                   remove.duplicated = FALSE,
                   fromLast = FALSE,
                   warn = TRUE)
lav_model_partable(
                   model = NULL,
                   meanstructure = FALSE,
                   int.ov.free = FALSE,
                   int.lv.free = FALSE,
                   marker.int.zero = FALSE,
                   orthogonal = FALSE,
                   orthogonal.y = FALSE,
                   orthogonal.x = FALSE,
                   orthogonal.efa = FALSE,
                   std.lv = FALSE,
                   correlation = FALSE,
                   composites = TRUE,
                   effect.coding = "",
                   conditional.x = FALSE,
                   fixed.x = FALSE,
                   parameterization = "delta",
                   constraints = NULL,
                   ceq.simple = FALSE,
                   auto = FALSE,
                   model.type = "sem",
                   auto.fix.first = FALSE,
                   auto.fix.single = FALSE,
                   auto.var = FALSE,
                   auto.cov.lv.x = FALSE,
                   auto.cov.y = FALSE,
                   auto.th = FALSE,
                   auto.delta = FALSE,
                   auto.efa = FALSE,
                   varTable = NULL,
                   ngroups = 1L,
                   nthresholds = NULL,
                   group.equal = NULL,
                   group.partial = NULL,
                   group.w.free = FALSE,
                   debug = FALSE,
                   warn = TRUE,
                   as.data.frame. = TRUE)
lav_matrix_vec(A)
lav_matrix_duplication(n = 1L)
lav_matrix_commutation(m = 1L, n = 1L)
lav_matrix_cov(Y, Mu = NULL)
lav_partable_ndat(partable)
lav_partable_npar(partable)
lav_partable_add(partable = NULL, add = list())
lav_partable_df(partable)
lav_scores(
           object,
           scaling = FALSE,
           ignore.constraints = FALSE,
           remove.duplicated = TRUE,
           remove.empty.cases = TRUE)
}
\arguments{
\item{...}{See argument \code{...} in the corresponding function with shortened name.}
\item{A}{See argument \code{a} in the corresponding function with shortened name.}
\item{add}{See argument \code{add} in the corresponding function with shortened name.}
\item{as.data.frame.}{See argument \code{as_data_frame} in the corresponding function with shortened name.}
\item{auto}{See argument \code{auto} in the corresponding function with shortened name.}
\item{auto.cov.lv.x}{See argument \code{auto_cov_lv_x} in the corresponding function with shortened name.}
\item{auto.cov.y}{See argument \code{auto_cov_y} in the corresponding function with shortened name.}
\item{auto.delta}{See argument \code{auto_delta} in the corresponding function with shortened name.}
\item{auto.efa}{See argument \code{auto_efa} in the corresponding function with shortened name.}
\item{auto.fix.first}{See argument \code{auto_fix_first} in the corresponding function with shortened name.}
\item{auto.fix.single}{See argument \code{auto_fix_single} in the corresponding function with shortened name.}
\item{auto.th}{See argument \code{auto_th} in the corresponding function with shortened name.}
\item{auto.var}{See argument \code{auto_var} in the corresponding function with shortened name.}
\item{blocks}{See argument \code{blocks} in the corresponding function with shortened name.}
\item{ceq.simple}{See argument \code{ceq_simple} in the corresponding function with shortened name.}
\item{check}{See argument \code{check} in the corresponding function with shortened name.}
\item{composites}{See argument \code{composites} in the corresponding function with shortened name.}
\item{con}{See argument \code{con} in the corresponding function with shortened name.}
\item{conditional.x}{See argument \code{conditional_x} in the corresponding function with shortened name.}
\item{constraints}{See argument \code{constraints} in the corresponding function with shortened name.}
\item{correlation}{See argument \code{correlation} in the corresponding function with shortened name.}
\item{debug}{See argument \code{debug} in the corresponding function with shortened name.}
\item{diagonal}{See argument \code{diagonal} in the corresponding function with shortened name.}
\item{effect.coding}{See argument \code{effect_coding} in the corresponding function with shortened name.}
\item{est}{See argument \code{est} in the corresponding function with shortened name.}
\item{fallback.simple}{See argument \code{fallback_simple} in the corresponding function with shortened name.}
\item{fixed.x}{See argument \code{fixed_x} in the corresponding function with shortened name.}
\item{fromLast}{See argument \code{from_last} in the corresponding function with shortened name.}
\item{func}{See argument \code{func} in the corresponding function with shortened name.}
\item{group.equal}{See argument \code{group_equal} in the corresponding function with shortened name.}
\item{group.partial}{See argument \code{group_partial} in the corresponding function with shortened name.}
\item{group.w.free}{See argument \code{group_w_free} in the corresponding function with shortened name.}
\item{h}{See argument \code{h} in the corresponding function with shortened name.}
\item{ignore.constraints}{See argument \code{ignore_constraints} in the corresponding function with shortened name.}
\item{int.lv.free}{See argument \code{int_lv_free} in the corresponding function with shortened name.}
\item{int.ov.free}{See argument \code{int_ov_free} in the corresponding function with shortened name.}
\item{label}{See argument \code{label} in the corresponding function with shortened name.}
\item{lavdata}{See argument \code{lavdata} in the corresponding function with shortened name.}
\item{lavh1}{See argument \code{lavh1} in the corresponding function with shortened name.}
\item{lavobject}{See argument \code{lavobject} in the corresponding function with shortened name.}
\item{lavoptions}{See argument \code{lavoptions} in the corresponding function with shortened name.}
\item{lavpta}{See argument \code{lavpta} in the corresponding function with shortened name.}
\item{lavsamplestats}{See argument \code{lavsamplestats} in the corresponding function with shortened name.}
\item{m}{See argument \code{m} in the corresponding function with shortened name.}
\item{marker.int.zero}{See argument \code{marker_int_zero} in the corresponding function with shortened name.}
\item{meanstructure}{See argument \code{meanstructure} in the corresponding function with shortened name.}
\item{model}{See argument \code{model} in the corresponding function with shortened name.}
\item{model.type}{See argument \code{model_type} in the corresponding function with shortened name.}
\item{Mu}{See argument \code{mu} in the corresponding function with shortened name.}
\item{n}{See argument \code{n} in the corresponding function with shortened name.}
\item{NACOV}{See argument \code{nacov} in the corresponding function with shortened name.}
\item{ngroups}{See argument \code{ngroups} in the corresponding function with shortened name.}
\item{nthresholds}{See argument \code{nthresholds} in the corresponding function with shortened name.}
\item{object}{See argument \code{object} in the corresponding function with shortened name.}
\item{orthogonal}{See argument \code{orthogonal} in the corresponding function with shortened name.}
\item{orthogonal.efa}{See argument \code{orthogonal_efa} in the corresponding function with shortened name.}
\item{orthogonal.x}{See argument \code{orthogonal_x} in the corresponding function with shortened name.}
\item{orthogonal.y}{See argument \code{orthogonal_y} in the corresponding function with shortened name.}
\item{parameterization}{See argument \code{parameterization} in the corresponding function with shortened name.}
\item{partable}{See argument \code{partable} in the corresponding function with shortened name.}
\item{pt1}{See argument \code{pt1} in the corresponding function with shortened name.}
\item{pt2}{See argument \code{pt2} in the corresponding function with shortened name.}
\item{pta}{See argument \code{pta} in the corresponding function with shortened name.}
\item{remove.duplicated}{See argument \code{remove_duplicated} in the corresponding function with shortened name.}
\item{remove.empty.cases}{See argument \code{remove_empty_cases} in the corresponding function with shortened name.}
\item{S}{See argument \code{s} in the corresponding function with shortened name.}
\item{sample.cov}{See argument \code{sample_cov} in the corresponding function with shortened name.}
\item{sample.cov.x}{See argument \code{sample_cov_x} in the corresponding function with shortened name.}
\item{sample.mean}{See argument \code{sample_mean} in the corresponding function with shortened name.}
\item{sample.mean.x}{See argument \code{sample_mean_x} in the corresponding function with shortened name.}
\item{sample.slopes}{See argument \code{sample_slopes} in the corresponding function with shortened name.}
\item{sample.th}{See argument \code{sample_th} in the corresponding function with shortened name.}
\item{sample.th.idx}{See argument \code{sample_th_idx} in the corresponding function with shortened name.}
\item{scaling}{See argument \code{scaling} in the corresponding function with shortened name.}
\item{start}{See argument \code{start} in the corresponding function with shortened name.}
\item{std.lv}{See argument \code{std_lv} in the corresponding function with shortened name.}
\item{theta}{See argument \code{theta} in the corresponding function with shortened name.}
\item{txtOnly}{See argument \code{txt_only} in the corresponding function with shortened name.}
\item{type}{See argument \code{type} in the corresponding function with shortened name.}
\item{varTable}{See argument \code{var_table} in the corresponding function with shortened name.}
\item{warn}{See argument \code{warn} in the corresponding function with shortened name.}
\item{WLS.V}{See argument \code{wls_v} in the corresponding function with shortened name.}
\item{x}{See argument \code{x} in the corresponding function with shortened name.}
\item{Y}{See argument \code{y} in the corresponding function with shortened name.}
}
\details{
Function names are shortened by following replacements:

\describe{
  \item{_mvnorm}{_mvn}
  \item{_cluster}{_cl}
  \item{_missing}{_mi}
  \item{_matrix}{_mat}
  \item{_duplication}{_dup}
  \item{_commutation}{_com}
  \item{_information}{_info}
  \item{_samplestats}{_samp}
  \item{_satorra_bentler}{_sb}
  \item{_SatorraBentler}{_sb}
  \item{_yuan_bentler}{_yb}
  \item{_inverted}{_inv}
  \item{_estimate}{_est}
  \item{_symmetric}{_sym}
  \item{_lavaan_step}{_step}
  \item{_univariate}{_uni}
  \item{_fit_measures}{_fit}
  \item{_constraints}{_con}
  \item{_scores}{_sc}
  \item{_reverse}{_rev}
  \item{_vech_sigma}{_sigma}
  \item{_vechsigma}{_sigma}
  \item{_object_inspect}{_inspect}
  \item{_orthogonal}{_ortho}
  \item{_gradient}{_grad}
  \item{_partable}{_pt}
}
}
