/*

  Code for the article 'Animal diversity and ecosystem functioning in dynamic food webs'

  Copyright (C) 2016 Christian Guill & Florian D. Schneider

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.


 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>      // for input/output from or to files
//#include <stdlib.h>
#include <iostream>		// for input/output on terminal
#include <sstream>
#include <vector>
#include <cstdlib>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>				// random number generator
#include <gsl/gsl_randist.h>			// random number distributions
#include <gsl/gsl_blas.h>				// linear algebra routines
#include <gsl/gsl_linalg.h>				// linear algebra
#include <gsl/gsl_sort_vector.h>		// vector operations
#include <gsl/gsl_odeiv.h>              // ODE solver

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */

double *web_calc(gsl_rng *r);
static void Generate_RGG_structure(gsl_rng *r, gsl_matrix *Loc1, gsl_matrix *Loc2, gsl_matrix *Loc3, gsl_matrix *Loc4, gsl_matrix *Loc5, gsl_matrix *Loc6, double a, double b, int n, double sigma);
static void output(gsl_matrix *Loc1, gsl_matrix *Loc2, gsl_matrix *Loc3, gsl_matrix *Loc4, gsl_matrix *Loc5, gsl_matrix *Loc6);


