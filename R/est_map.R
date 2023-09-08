#' @export
est_map_R <- function(m,
                      n.mrk,
                      n.ind,
                      emit,
                      rf_vec,
                      tol = 10e-3) {
  rf_vec[rf_vec < 10e-10] <- 10e-10
  res.temp <- .Call("est_hmm_map",
                    m,
                    n.mrk,
                    n.ind,
                    emit,
                    rf_vec,
                    tol,
                    verbose = TRUE,
                    PACKAGE = "highprecHMM")
  res.temp
}



#' @export
est_map_R_no_log <- function(m,
                      n.mrk,
                      n.ind,
                      emit,
                      rf_vec,
                      tol = 10e-3) {
  rf_vec[rf_vec < 10e-10] <- 10e-10
  res.temp <- .Call("est_hmm_map_no_log",
                    m,
                    n.mrk,
                    n.ind,
                    emit,
                    rf_vec,
                    tol,
                    verbose = TRUE,
                    PACKAGE = "highprecHMM")
  res.temp
}
