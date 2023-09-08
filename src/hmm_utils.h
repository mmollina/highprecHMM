/*
 MAPpoly: a package to construct genetic maps in autopolyploids
 Copyright (C) 2014-2020 Marcelo Mollinari

 This file is part of MAPpoly.

 MAPpoly is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 For a copy of the GNU General Public License, please visit
 <http://www.gnu.org/licenses/>.
 */


/*
  File: combinatorial_functions.cpp
  Description: Set of combinatorial functions to be used with mappoly

  Functions Written by Marcelo Mollinari.

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
  First version: Dec 19, 2013
  Last update:   Sep 24, 2020
*/

#include <R.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#define THRESH 200.0

int nChoosek(int n, int k);
int n_rec_given_genk_and_k1(int ploidy, int index1, int index2);
double prob_k1_given_k_lp_lq_m(int m,
                               int lp,
                               int lq,
                               double rf);

std::vector <int>  boolean_lexicographic_k_choose_m_and_collapse(int ploidy,
        std::vector<int>& which_homologous_mk1,
        std::vector<int>& which_homologous_mk2,
        int gen_prog_mk1,
        int gen_prog_mk2);

std::vector <bool> get_boolean_vec_from_lexicographical_index(int ploidy, int index);

double addlog(double a, double b);

double log_prob_k1_given_k_l_m(int m, int l, double rf);

double prob_k1_given_k_l_m(int m, int l, double rf);

std::vector<double> transition(int m, std::vector<double>& rf);

std::vector<double> log_transition(int m, std::vector<double>& rf);

std::vector<double> rec_num(int m);
//end of file
