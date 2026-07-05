/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_MATH_EIGEN_H
#define SPARTA_MATH_EIGEN_H

namespace MathEigen {

/** A specialized function which finds the eigenvalues and eigenvectors
 *  of a 3x3 matrix (in double ** format).
 *
 * \param  mat   the 3x3 matrix you wish to diagonalize
 * \param  eval  store the eigenvalues here
 * \param  evec  store the eigenvectors here...
 * \return       0 if eigenvalue calculation converged, 1 if it failed */

int jacobi3(double const *const *mat, double *eval, double **evec);

/** \overload */

int jacobi3(double const mat[3][3], double *eval, double evec[3][3]);

}    // namespace MathEigen

#endif    //#ifndef SPARTA_MATH_EIGEN_H
