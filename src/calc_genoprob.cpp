#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "hmm_utils.h"
#include <thread>
#include <chrono>
#define THRESH 200.0

using namespace std;
using namespace Rcpp;


/* FUNCTION: calc_genprob
 * -----------------------------------------------------
 */
RcppExport SEXP calc_genprob(SEXP ploidyR,
                                 SEXP n_mrkR,
                                 SEXP n_indR,
                                 SEXP log_emitR,
                                 SEXP rfR,
                                 SEXP probsR,
                                 SEXP verboseR)
{
  //convert input to C++ types
  int m = Rcpp::as<int>(ploidyR);
  int n_mrk = Rcpp::as<int>(n_mrkR);
  int n_ind = Rcpp::as<int>(n_indR);
  int n_gen = pow(nChoosek(m, m/2),2);
  Rcpp::NumericVector log_emit(log_emitR);
  std::vector<double> alpha(log_emit.size());
  std::vector<double> beta(log_emit.size());
  std::vector<double> T;
  std::vector<double> rf = Rcpp::as<std::vector<double> >(rfR);
  int verbose = Rcpp::as<int>(verboseR);
  std::vector<double> probs = Rcpp::as<std::vector<double> >(probsR);
  int k, k1,count = 0, count2 = 0;
  std::vector<double> term(n_ind);
  double loglike = 0;
  if(verbose)
  {
    Rcpp::Rcout << "\tPloidy level: " << m << "\n" ;
    Rcpp::Rcout << "\tNumber of markers: " << n_mrk << "\n" ;
    Rcpp::Rcout << "\tNumber of individuals: " << n_ind << "\n\t" << 1 << " ";
  }
  /* TRANSITION */
  // i: recombination fraction (rf[0]: between first and second)
  // j: state in
  // k: state out
  // T[i*n_gen*n_gen + j*n_gen + k]
  T = log_transition(m, rf);

  /* EMISSION */
  // i: individual
  // j: marker
  // k: genotype state
  // log_emit(i*n_mrk*n_gen + j*n_gen + k)

  //Initializing alpha and beta
  for(int i=0; i < n_ind; i++)
  {
    R_CheckUserInterrupt();
    if(verbose){
      Rcpp::Rcout << ".";
      if((i+1)%50 == 0) Rcpp::Rcout << " " << i << "\n\t" << i+1 << " ";
    }
    for(int k = 0; (unsigned)k < n_gen; k++) {
      alpha[i*n_mrk*n_gen + k] = log_emit(i*n_mrk*n_gen + k); // fist marker: j*n_gen = 0*n_gen
      beta[i*n_mrk*n_gen + (n_mrk-1)*n_gen + k] = 0.0; // last marker: (n_mrk-1)*n_gen
    }

    //forward-backward
    for(k = 1, k1 = n_mrk-2; k < n_mrk; k++, k1--)
    {
      for(int j1 = 0; j1 < n_gen; j1++ )
      {
        alpha[i*n_mrk*n_gen + k*n_gen + j1] = alpha[i*n_mrk*n_gen + (k-1)*n_gen] +
                                              T[(k-1)*n_gen*n_gen + j1];
        beta[i*n_mrk*n_gen + k1*n_gen + j1] = beta[i*n_mrk*n_gen + (k1+1)*n_gen] +
                                              T[k1*n_gen*n_gen + j1*n_gen] +
                                              log_emit(i*n_mrk*n_gen + (k1+1)*n_gen);

        for(int j = 1; j < n_gen; j++ )
        {
          alpha[i*n_mrk*n_gen + k*n_gen + j1] = addlog(alpha[i*n_mrk*n_gen + k*n_gen + j1],
                                                       alpha[i*n_mrk*n_gen +  (k-1)*n_gen + j] +
                                                         T[(k-1)*n_gen*n_gen + j*n_gen + j1]);

          beta[i*n_mrk*n_gen + k1*n_gen + j1] = addlog(beta[i*n_mrk*n_gen + k1*n_gen + j1],
                                                       beta[i*n_mrk*n_gen + (k1+1)*n_gen + j] +
                                                         T[k1*n_gen*n_gen + j1*n_gen + j] +
                                                         log_emit(i*n_mrk*n_gen + (k1+1)*n_gen + j));
        }
        alpha[i*n_mrk*n_gen + k*n_gen + j1] += log_emit(i*n_mrk*n_gen + k*n_gen + j1);
      }
    }
    /* calculate genotype probabilities */
    for(int k=0; k < n_mrk; k++)
    {
      double w = alpha[i*n_mrk*n_gen + k*n_gen] +
        beta[i*n_mrk*n_gen + k*n_gen];
      for(int j=1; (unsigned)j < n_gen; j++)
        w = addlog(w, alpha[i*n_mrk*n_gen + k*n_gen +j] +
          beta[i*n_mrk*n_gen + k*n_gen +j]);
      for(int j=0; (unsigned)j < n_gen; j++){
        probs[count] = exp(alpha[i*n_mrk*n_gen + k*n_gen +j] +
          beta[i*n_mrk*n_gen + k*n_gen +j] - w);
        count++;
      }
    }
    //Termination
    term[i] = alpha[i*n_mrk*n_gen + (n_mrk-1)*n_gen];
    for(int j=1; (unsigned)j < n_gen; j++)
      term[i] = addlog(term[i], alpha[i*n_mrk*n_gen + (n_mrk-1)*n_gen + j]);
    loglike += term[i];
  }//Loop over all individuals
  List z  = List::create(probs, loglike);
  return z ;
}

/* FUNCTION: calc_genprob_no_log
 * -----------------------------------------------------
 */
RcppExport SEXP calc_genprob_no_log(SEXP ploidyR,
                                    SEXP n_mrkR,
                                    SEXP n_indR,
                                    SEXP emitR,
                                    SEXP rfR,
                                    SEXP probsR,
                                    SEXP verboseR)
{
  //convert input to C++ types
  int m = Rcpp::as<int>(ploidyR);
  int n_mrk = Rcpp::as<int>(n_mrkR);
  int n_ind = Rcpp::as<int>(n_indR);
  int n_gen = pow(nChoosek(m, m/2),2);
  Rcpp::NumericVector emit(emitR);
  std::vector<long double> alpha(emit.size());
  std::vector<long double> beta(emit.size());
  std::vector<double> T;
  std::vector<double> rf = Rcpp::as<std::vector<double> >(rfR);
  int verbose = Rcpp::as<int>(verboseR);
  std::vector<long double> probs = Rcpp::as<std::vector<long double> >(probsR);
  int k, k1,count = 0, count2 = 0;
  std::vector<double> term(n_ind);
  double loglike = 0;
  if(verbose)
  {
    Rcpp::Rcout << "\tPloidy level: " << m << "\n" ;
    Rcpp::Rcout << "\tNumber of markers: " << n_mrk << "\n" ;
    Rcpp::Rcout << "\tNumber of individuals: " << n_ind << "\n\t" << 1 << " ";
  }
  /* TRANSITION */
  // i: recombination fraction (rf[0]: between first and second)
  // j: state in
  // k: state out
  // T[i*n_gen*n_gen + j*n_gen + k]
  T = transition(m, rf);

  /* EMISSION */
  // i: individual
  // j: marker
  // k: genotype state
  // emit(i*n_mrk*n_gen + j*n_gen + k)

  //Initializing alpha and beta
  for(int i=0; i < n_ind; i++)
  {
    R_CheckUserInterrupt();
    if(verbose){
      Rcpp::Rcout << ".";
      if((i+1)%50 == 0) Rcpp::Rcout << " " << i << "\n\t" << i+1 << " ";
    }
    for(int k = 0; (unsigned)k < n_gen; k++) {
      alpha[i*n_mrk*n_gen + k] = emit(i*n_mrk*n_gen + k); // fist marker: j*n_gen = 0*n_gen
      beta[i*n_mrk*n_gen + (n_mrk-1)*n_gen + k] = 1.0; // last marker: (n_mrk-1)*n_gen
    }

    //forward-backward
    for(k = 1, k1 = n_mrk-2; k < n_mrk; k++, k1--)
    {
      for(int j1 = 0; j1 < n_gen; j1++ )
      {
        alpha[i*n_mrk*n_gen + k*n_gen + j1] = alpha[i*n_mrk*n_gen + (k-1)*n_gen] *
          T[(k-1)*n_gen*n_gen + j1];
        beta[i*n_mrk*n_gen + k1*n_gen + j1] = beta[i*n_mrk*n_gen + (k1+1)*n_gen] *
          T[k1*n_gen*n_gen + j1*n_gen] *
          emit(i*n_mrk*n_gen + (k1+1)*n_gen);

        for(int j = 1; j < n_gen; j++ )
        {
          alpha[i*n_mrk*n_gen + k*n_gen + j1] = alpha[i*n_mrk*n_gen + k*n_gen + j1] +
                                                       alpha[i*n_mrk*n_gen +  (k-1)*n_gen + j] *
                                                         T[(k-1)*n_gen*n_gen + j*n_gen + j1];

          beta[i*n_mrk*n_gen + k1*n_gen + j1] = beta[i*n_mrk*n_gen + k1*n_gen + j1] +
                                                       beta[i*n_mrk*n_gen + (k1+1)*n_gen + j] *
                                                         T[k1*n_gen*n_gen + j1*n_gen + j] *
                                                         emit(i*n_mrk*n_gen + (k1+1)*n_gen + j);
        }
        alpha[i*n_mrk*n_gen + k*n_gen + j1] *= emit(i*n_mrk*n_gen + k*n_gen + j1);
      }
    }
    /* calculate genotype probabilities */
    for(int k=0; k < n_mrk; k++)
    {
      long double w = alpha[i*n_mrk*n_gen + k*n_gen] *
        beta[i*n_mrk*n_gen + k*n_gen];
      for(int j=1; (unsigned)j < n_gen; j++)
        w = w + alpha[i*n_mrk*n_gen + k*n_gen +j] *
          beta[i*n_mrk*n_gen + k*n_gen +j];
      for(int j=0; (unsigned)j < n_gen; j++){
        probs[count] = alpha[i*n_mrk*n_gen + k*n_gen +j] *
          beta[i*n_mrk*n_gen + k*n_gen +j] / w;
        count++;
      }
    }
    //Termination
    term[i] = alpha[i*n_mrk*n_gen + (n_mrk-1)*n_gen];
    for(int j=1; (unsigned)j < n_gen; j++)
      term[i] = term[i] + alpha[i*n_mrk*n_gen + (n_mrk-1)*n_gen + j];
    loglike += log(term[i]);
  }//Loop over all individuals
  List z  = List::create(probs, loglike);
  return z ;
}
