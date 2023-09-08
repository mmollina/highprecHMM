#' @export
calc_genoprob_R <- function(m,
                            n.mrk,
                            n.ind,
                            mrknames,
                            indnames,
                            emit,
                            rf_vec) {
  rf_vec[rf_vec < 10e-10] <- 10e-10
  res.temp <- .Call("calc_genprob",
                    m,
                    n.mrk,
                    n.ind,
                    emit,
                    rf_vec,
                    as.numeric(rep(0, choose(m, m/2)^2 * n.mrk * n.ind)),
                    verbose = TRUE,
                    PACKAGE = "highprecHMM")
    dim(res.temp[[1]])<-c(choose(m,m/2)^2,n.mrk,n.ind)
    dimnames(res.temp[[1]])<-list(kronecker(apply(combn(letters[1:m],m/2),2, paste, collapse=""),
                                            apply(combn(letters[(m+1):(2*m)],m/2),2, paste, collapse=""), paste, sep=":"), mrknames, indnames)
    map <- cumsum(c(0, mappoly::imf_h(rf_vec)))
    names(map) <- mrknames
    structure(list(probs = res.temp[[1]], map = map),
              class = "mappoly.genoprob")
}

#' @export
calc_genoprob_R_no_log <- function(m,
                            n.mrk,
                            n.ind,
                            mrknames,
                            indnames,
                            emit,
                            rf_vec) {
  rf_vec[rf_vec < 10e-10] <- 10e-10
  res.temp <- .Call("calc_genprob_no_log",
                    m,
                    n.mrk,
                    n.ind,
                    emit,
                    rf_vec,
                    as.numeric(rep(0, choose(m, m/2)^2 * n.mrk * n.ind)),
                    verbose = TRUE,
                    PACKAGE = "highprecHMM")
  dim(res.temp[[1]])<-c(choose(m,m/2)^2,n.mrk,n.ind)
  dimnames(res.temp[[1]])<-list(kronecker(apply(combn(letters[1:m],m/2),2, paste, collapse=""),
                                          apply(combn(letters[(m+1):(2*m)],m/2),2, paste, collapse=""), paste, sep=":"), mrknames, indnames)
  map <- cumsum(c(0, mappoly::imf_h(rf_vec)))
  names(map) <- mrknames
  structure(list(probs = res.temp[[1]], map = map),
            class = "mappoly.genoprob")
}
