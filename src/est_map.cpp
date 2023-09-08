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

/* FUNCTION: est_hmm_map
 * -----------------------------------------------------
 */
RcppExport SEXP est_hmm_map(SEXP ploidyR,
                                 SEXP n_mrkR,
                                 SEXP n_indR,
                                 SEXP log_emitR,
                                 SEXP rfR,
                                 SEXP tolR,
                                 SEXP verboseR){
  //convert input to C++ types
  int m = Rcpp::as<int>(ploidyR);
  int n_mrk = Rcpp::as<int>(n_mrkR);
  int n_ind = Rcpp::as<int>(n_indR);
  int n_gen = pow(nChoosek(m, m/2),2);
  Rcpp::NumericVector log_emit(log_emitR);
  std::vector<double> alpha(log_emit.size());
  std::vector<double> beta(log_emit.size());
  vector<double> gamma(n_gen*n_gen);
  std::vector<double> T;
  std::vector<double> rf = Rcpp::as<std::vector<double> >(rfR);
  int verbose = Rcpp::as<int>(verboseR);
  double tol = Rcpp::as<double>(tolR);
  int k, k1, maxit = 1000, flag;
  std::vector<double> term(n_ind), rf_cur(rf.size());
  double loglike = 0;
  if(verbose)
  {
    Rcpp::Rcout << "\tPloidy level: " << m << "\n" ;
    Rcpp::Rcout << "\tNumber of markers: " << n_mrk << "\n" ;
    Rcpp::Rcout << "\tNumber of individuals: " << n_ind << "\n" ;
  }
  /* NUMBER OF RECOMBINATIONS*/
  // i: state in
  // j: state out
  // R[i*n_gen + j]
  std::vector<double> R;
  R = rec_num(m);

  //begin EM algorithm
  for(int it=0; it<maxit; it++)
  {
    //Initializing recombination fraction vector for Baum-Welch
    for(int j=0; j<n_mrk-1; j++)
    {
      rf_cur[j] = rf[j];
      rf[j] = 0.0;
    }

    /* TRANSITION */
    // i: recombination fraction (rf[0]: between first and second)
    // j: state in
    // k: state out
    // T[i*n_gen*n_gen + j*n_gen + k]
    T = log_transition(m, rf_cur);

    /* EMISSION */
    // i: individual
    // j: marker
    // k: genotype state
    // log_emit(i*n_mrk*n_gen + j*n_gen + k)

    //Start loop over all individuals
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

    //Updating recombination fraction
    for(int j = 0; j < n_mrk-1; j++)
    {
      double s = 0.0;
      for(int v = 0; v < n_gen; v++)
      {
        for(int v2 = 0; v2 < n_gen; v2++)
        {
          gamma[v*n_gen + v2] = alpha[i*n_mrk*n_gen + j*n_gen + v] +
                                beta[i*n_mrk*n_gen + (j+1)*n_gen + v2] +
                                T[j*n_gen*n_gen + v*n_gen + v2] +
                                log_emit(i*n_mrk*n_gen + (j+1)*n_gen + v2);
          if(v==0 && v2==0) s = gamma[v*n_gen + v2];
          else s = addlog(s, gamma[v*n_gen + v2]);
        }
      }
      for(int v = 0; v < n_gen; v++)
      {
        for(int v2 = 0; v2 < n_gen; v2++)
        {
          if(s != 0){
            rf[j] +=  R[v*n_gen + v2] * exp(gamma[v*n_gen + v2] - s);
            }
        }
      }
    }
  }//End loop over all individuals

  // rescale
  for(int j=0; j<n_mrk-1; j++)
  {
    rf[j] /= (double)n_ind;
    if(rf[j] < tol/100.0) rf[j] = tol/100.0;
    else if(rf[j] > 0.5-tol/100.0) rf[j] = 0.5-tol/100.0;
  }
  // check convergence
  flag=0;
  for(int j=0; j < n_mrk-1; j++)
  {
    if(fabs(rf[j] - rf_cur[j]) > tol*(rf_cur[j]+tol*100.0))
    {
      flag = 1;
      break;
    }
  }
  if(verbose)
  {
    Rcpp::Rcout << "\t\n Iter: " << it+1 << "\t";
    for(int j = 0; j < n_mrk-1; j++)
    {
      Rcpp::Rcout.precision(3);
      Rcpp::Rcout << std::fixed << rf[j] << " ";
    }
  }
  if(!flag) break;
  }//end of EM algorithm
  if(flag && verbose) Rcpp::Rcout << "Didn't converge!\n";

  //Termination
  for(int i=0; i < n_ind; i++){
  term[i] = alpha[i*n_mrk*n_gen + (n_mrk-1)*n_gen];
  for(int j=1; (unsigned)j < n_gen; j++)
    term[i] = addlog(term[i], alpha[i*n_mrk*n_gen + (n_mrk-1)*n_gen + j]);
  loglike += term[i];
  }
  if(verbose)
    Rcpp::Rcout << "\n";
  List z = List::create(wrap(loglike), rf);
  return(z);
}
/* FUNCTION: est_hmm_map_no_log
 * -----------------------------------------------------
 */
RcppExport SEXP est_hmm_map_no_log(SEXP ploidyR,
                                   SEXP n_mrkR,
                                   SEXP n_indR,
                                   SEXP emitR,
                                   SEXP rfR,
                                   SEXP tolR,
                                   SEXP verboseR){
  //convert input to C++ types
  int m = Rcpp::as<int>(ploidyR);
  int n_mrk = Rcpp::as<int>(n_mrkR);
  int n_ind = Rcpp::as<int>(n_indR);
  int n_gen = pow(nChoosek(m, m/2),2);
  Rcpp::NumericVector emit(emitR);
  std::vector<long double> alpha(emit.size());
  std::vector<long double> beta(emit.size());
  vector<long double> gamma(n_gen*n_gen);
  std::vector<double> T;
  std::vector<long double> rf = Rcpp::as<std::vector<long double> >(rfR);
  int verbose = Rcpp::as<int>(verboseR);
  double tol = Rcpp::as<double>(tolR);
  int k, k1, maxit = 1000, flag;
  std::vector<double> term(n_ind), rf_cur(rf.size());
  double loglike = 0;
  if(verbose)
  {
    Rcpp::Rcout << "\tPloidy level: " << m << "\n" ;
    Rcpp::Rcout << "\tNumber of markers: " << n_mrk << "\n" ;
    Rcpp::Rcout << "\tNumber of individuals: " << n_ind << "\n" ;
  }
  /* NUMBER OF RECOMBINATIONS*/
  // i: state in
  // j: state out
  // R[i*n_gen + j]
  std::vector<double> R;
  R = rec_num(m);

  //begin EM algorithm
  for(int it=0; it<maxit; it++)
  {
    //Initializing recombination fraction vector for Baum-Welch
    for(int j=0; j<n_mrk-1; j++)
    {
      rf_cur[j] = rf[j];
      rf[j] = 0.0;
    }

    /* TRANSITION */
    // i: recombination fraction (rf[0]: between first and second)
    // j: state in
    // k: state out
    // T[i*n_gen*n_gen + j*n_gen + k]
    T = transition(m, rf_cur);

    /* EMISSION */
    // i: individual
    // j: marker
    // k: genotype state
    // emit(i*n_mrk*n_gen + j*n_gen + k)

    //Start loop over all individuals
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

      //Updating recombination fraction
      for(int j = 0; j < n_mrk-1; j++)
      {
        long double s = 0.0;
        for(int v = 0; v < n_gen; v++)
        {
          for(int v2 = 0; v2 < n_gen; v2++)
          {
            gamma[v*n_gen + v2] = alpha[i*n_mrk*n_gen + j*n_gen + v] *
              beta[i*n_mrk*n_gen + (j+1)*n_gen + v2] *
              T[j*n_gen*n_gen + v*n_gen + v2] *
              emit(i*n_mrk*n_gen + (j+1)*n_gen + v2);
            if(v==0 && v2==0) s = gamma[v*n_gen + v2];
            else s += gamma[v*n_gen + v2];
          }
        }
        for(int v = 0; v < n_gen; v++)
        {
          for(int v2 = 0; v2 < n_gen; v2++)
          {
            if(s != 0){
              rf[j] +=  R[v*n_gen + v2] * gamma[v*n_gen + v2]/s;
            }
          }
        }
      }
    }//End loop over all individuals

    // re-scale
    for(int j=0; j<n_mrk-1; j++)
    {
      rf[j] /= (double)n_ind;
      if(rf[j] < tol/100.0) rf[j] = tol/100.0;
      else if(rf[j] > 0.5-tol/100.0) rf[j] = 0.5-tol/100.0;
    }
    // check convergence
    flag=0;
    for(int j=0; j < n_mrk-1; j++)
    {
      if(fabs(rf[j] - rf_cur[j]) > tol*(rf_cur[j]+tol*100.0))
      {
        flag = 1;
        break;
      }
    }
    if(verbose)
    {
      Rcpp::Rcout << "\t\n Iter: " << it+1 << "\t";
      for(int j = 0; j < n_mrk-1; j++)
      {
        Rcpp::Rcout.precision(3);
        Rcpp::Rcout << std::fixed << rf[j] << " ";
      }
    }
    if(!flag) break;
  }//end of EM algorithm
  if(flag && verbose) Rcpp::Rcout << "Didn't converge!\n";

  //Termination
  for(int i=0; i < n_ind; i++){
    term[i] = alpha[i*n_mrk*n_gen + (n_mrk-1)*n_gen];
    for(int j=1; (unsigned)j < n_gen; j++)
      term[i] += alpha[i*n_mrk*n_gen + (n_mrk-1)*n_gen + j];
    loglike += term[i];
  }

  if(verbose)
    Rcpp::Rcout << "\n";
  List z = List::create(wrap(loglike), rf);
  return(z);
}
//end of file
