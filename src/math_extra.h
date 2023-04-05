/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#ifndef SPARTA_MATH_EXTRA_H
#define SPARTA_MATH_EXTRA_H

#include "spatype.h"
#include "math.h"
#include "stdio.h"
#include "string.h"

namespace MathExtra {

  // 3 vector operations

  inline void norm3(double *v);
  inline void normalize3(const double *v, double *ans);
  inline void snorm3(const double, double *v);
  inline void snormalize3(const double, const double *v, double *ans);
  inline void negate3(double *v);
  inline void scale3(double s, double *v);
  inline void scale3(double s, const double *v, double *ans);
  inline void axpy3(double alpha, const double *x, double *y);
  inline void axpy3(double alpha, const double *x, const double *y,
                    double *ynew);
  inline void add3(const double *v1, const double *v2, double *ans);
  inline void sub3(const double *v1, const double *v2, double *ans);
  inline double len3(const double *v);
  inline double lensq3(const double *v);
  inline double dot3(const double *v1, const double *v2);
  inline void cross3(const double *v1, const double *v2, double *ans);
  inline void reflect3(double *v, const double *norm);

  // 3x3 matrix operations

  inline double det3(const double mat[3][3]);
  inline void diag_times3(const double *diagonal, const double mat[3][3],
                          double ans[3][3]);
  inline void plus3(const double m[3][3], const double m2[3][3],
                    double ans[3][3]);
  inline void times3(const double m[3][3], const double m2[3][3],
                     double ans[3][3]);
  inline void transpose_times3(const double mat1[3][3],
                               const double mat2[3][3],
                               double ans[3][3]);
  inline void times3_transpose(const double mat1[3][3],
                               const double mat2[3][3],
                               double ans[3][3]);
  inline void invert3(const double mat[3][3], double ans[3][3]);
  inline void matvec(const double mat[3][3], const double*vec, double *ans);
  inline void matvec(const double *ex, const double *ey, const double *ez,
                     const double *vec, double *ans);
  inline void transpose_matvec(const double mat[3][3], const double*vec,
                               double *ans);
  inline void transpose_matvec(const double *ex, const double *ey,
                               const double *ez, const double *v,
                               double *ans);
  inline void transpose_diag3(const double mat[3][3], const double*vec,
                              double ans[3][3]);
  inline void vecmat(const double *v, const double m[3][3], double *ans);
  inline void scalar_times3(const double f, double m[3][3]);

  // quaternion operations

  inline void axisangle_to_quat(const double *v, const double angle,
                                double *quat);
  inline void quat_to_mat(const double *quat, double mat[3][3]);

  // misc methods

  int bounds(char *, int, int &, int &);
  char *num2str(SPARTA_NS::bigint, char *);
}

/* ----------------------------------------------------------------------
   normalize a vector in place
------------------------------------------------------------------------- */

void MathExtra::norm3(double *v)
{
  double scale = 1.0/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] *= scale;
  v[1] *= scale;
  v[2] *= scale;
}

/* ----------------------------------------------------------------------
   normalize a vector, return in ans
------------------------------------------------------------------------- */

void MathExtra::normalize3(const double *v, double *ans)
{
  double scale = 1.0/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  ans[0] = v[0]*scale;
  ans[1] = v[1]*scale;
  ans[2] = v[2]*scale;
}

/* ----------------------------------------------------------------------
   scale a vector to length in place
------------------------------------------------------------------------- */

void MathExtra::snorm3(const double length, double *v)
{
  double scale = length/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] *= scale;
  v[1] *= scale;
  v[2] *= scale;
}

/* ----------------------------------------------------------------------
   scale a vector to length
------------------------------------------------------------------------- */

void MathExtra::snormalize3(const double length, const double *v, double *ans)
{
  double scale = length/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  ans[0] = v[0]*scale;
  ans[1] = v[1]*scale;
  ans[2] = v[2]*scale;
}

/* ----------------------------------------------------------------------
   negate vector v in place
------------------------------------------------------------------------- */

void MathExtra::negate3(double *v)
{
  v[0] = -v[0];
  v[1] = -v[1];
  v[2] = -v[2];
}

/* ----------------------------------------------------------------------
   scale vector v by s in place
------------------------------------------------------------------------- */

void MathExtra::scale3(double s, double *v)
{
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}

/* ----------------------------------------------------------------------
   scale vector v by s, return in ans
------------------------------------------------------------------------- */

void MathExtra::scale3(double s, const double *v, double *ans)
{
  ans[0] = s*v[0];
  ans[1] = s*v[1];
  ans[2] = s*v[2];
}

/* ----------------------------------------------------------------------
   axpy: y = alpha*x + y
   y is replaced by result
------------------------------------------------------------------------- */

void MathExtra::axpy3(double alpha, const double *x, double *y)
{
  y[0] += alpha*x[0];
  y[1] += alpha*x[1];
  y[2] += alpha*x[2];
}

/* ----------------------------------------------------------------------
   axpy: ynew = alpha*x + y
------------------------------------------------------------------------- */

void MathExtra::axpy3(double alpha, const double *x, const double *y,
                      double *ynew)
{
  ynew[0] += alpha*x[0] + y[0];
  ynew[1] += alpha*x[1] + y[1];
  ynew[2] += alpha*x[2] + y[2];
}

/* ----------------------------------------------------------------------
   ans = v1 + v2
------------------------------------------------------------------------- */

void MathExtra::add3(const double *v1, const double *v2, double *ans)
{
  ans[0] = v1[0] + v2[0];
  ans[1] = v1[1] + v2[1];
  ans[2] = v1[2] + v2[2];
}

/* ----------------------------------------------------------------------
   ans = v1 - v2
------------------------------------------------------------------------- */

void MathExtra::sub3(const double *v1, const double *v2, double *ans)
{
  ans[0] = v1[0] - v2[0];
  ans[1] = v1[1] - v2[1];
  ans[2] = v1[2] - v2[2];
}

/* ----------------------------------------------------------------------
   length of vector v
------------------------------------------------------------------------- */

double MathExtra::len3(const double *v)
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/* ----------------------------------------------------------------------
   squared length of vector v, or dot product of v with itself
------------------------------------------------------------------------- */

double MathExtra::lensq3(const double *v)
{
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

/* ----------------------------------------------------------------------
   dot product of 2 vectors
------------------------------------------------------------------------- */

double MathExtra::dot3(const double *v1, const double *v2)
{
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

/* ----------------------------------------------------------------------
   cross product of 2 vectors
------------------------------------------------------------------------- */

void MathExtra::cross3(const double *v1, const double *v2, double *ans)
{
  ans[0] = v1[1]*v2[2] - v1[2]*v2[1];
  ans[1] = v1[2]*v2[0] - v1[0]*v2[2];
  ans[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

/* ----------------------------------------------------------------------
   reflect vector v around unit normal n
   return updated v of same length = v - 2(v dot n)n
------------------------------------------------------------------------- */

void MathExtra::reflect3(double *v, const double *n)
{
  double dot = dot3(v,n);
  v[0] -= 2.0*dot*n[0];
  v[1] -= 2.0*dot*n[1];
  v[2] -= 2.0*dot*n[2];
}

/* ----------------------------------------------------------------------
   determinant of a matrix
------------------------------------------------------------------------- */

double MathExtra::det3(const double m[3][3])
{
  double ans = m[0][0]*m[1][1]*m[2][2] - m[0][0]*m[1][2]*m[2][1] -
    m[1][0]*m[0][1]*m[2][2] + m[1][0]*m[0][2]*m[2][1] +
    m[2][0]*m[0][1]*m[1][2] - m[2][0]*m[0][2]*m[1][1];
  return ans;
}

/* ----------------------------------------------------------------------
   diagonal matrix times a full matrix
------------------------------------------------------------------------- */

void MathExtra::diag_times3(const double *d, const double m[3][3],
                            double ans[3][3])
{
  ans[0][0] = d[0]*m[0][0];
  ans[0][1] = d[0]*m[0][1];
  ans[0][2] = d[0]*m[0][2];
  ans[1][0] = d[1]*m[1][0];
  ans[1][1] = d[1]*m[1][1];
  ans[1][2] = d[1]*m[1][2];
  ans[2][0] = d[2]*m[2][0];
  ans[2][1] = d[2]*m[2][1];
  ans[2][2] = d[2]*m[2][2];
}

/* ----------------------------------------------------------------------
   add two matrices
------------------------------------------------------------------------- */

void MathExtra::plus3(const double m[3][3], const double m2[3][3],
                      double ans[3][3])
{
  ans[0][0] = m[0][0]+m2[0][0];
  ans[0][1] = m[0][1]+m2[0][1];
  ans[0][2] = m[0][2]+m2[0][2];
  ans[1][0] = m[1][0]+m2[1][0];
  ans[1][1] = m[1][1]+m2[1][1];
  ans[1][2] = m[1][2]+m2[1][2];
  ans[2][0] = m[2][0]+m2[2][0];
  ans[2][1] = m[2][1]+m2[2][1];
  ans[2][2] = m[2][2]+m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply mat1 times mat2
------------------------------------------------------------------------- */

void MathExtra::times3(const double m[3][3], const double m2[3][3],
                       double ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0] + m[0][1]*m2[1][0] + m[0][2]*m2[2][0];
  ans[0][1] = m[0][0]*m2[0][1] + m[0][1]*m2[1][1] + m[0][2]*m2[2][1];
  ans[0][2] = m[0][0]*m2[0][2] + m[0][1]*m2[1][2] + m[0][2]*m2[2][2];
  ans[1][0] = m[1][0]*m2[0][0] + m[1][1]*m2[1][0] + m[1][2]*m2[2][0];
  ans[1][1] = m[1][0]*m2[0][1] + m[1][1]*m2[1][1] + m[1][2]*m2[2][1];
  ans[1][2] = m[1][0]*m2[0][2] + m[1][1]*m2[1][2] + m[1][2]*m2[2][2];
  ans[2][0] = m[2][0]*m2[0][0] + m[2][1]*m2[1][0] + m[2][2]*m2[2][0];
  ans[2][1] = m[2][0]*m2[0][1] + m[2][1]*m2[1][1] + m[2][2]*m2[2][1];
  ans[2][2] = m[2][0]*m2[0][2] + m[2][1]*m2[1][2] + m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply the transpose of mat1 times mat2
------------------------------------------------------------------------- */

void MathExtra::transpose_times3(const double m[3][3], const double m2[3][3],
                                 double ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0] + m[1][0]*m2[1][0] + m[2][0]*m2[2][0];
  ans[0][1] = m[0][0]*m2[0][1] + m[1][0]*m2[1][1] + m[2][0]*m2[2][1];
  ans[0][2] = m[0][0]*m2[0][2] + m[1][0]*m2[1][2] + m[2][0]*m2[2][2];
  ans[1][0] = m[0][1]*m2[0][0] + m[1][1]*m2[1][0] + m[2][1]*m2[2][0];
  ans[1][1] = m[0][1]*m2[0][1] + m[1][1]*m2[1][1] + m[2][1]*m2[2][1];
  ans[1][2] = m[0][1]*m2[0][2] + m[1][1]*m2[1][2] + m[2][1]*m2[2][2];
  ans[2][0] = m[0][2]*m2[0][0] + m[1][2]*m2[1][0] + m[2][2]*m2[2][0];
  ans[2][1] = m[0][2]*m2[0][1] + m[1][2]*m2[1][1] + m[2][2]*m2[2][1];
  ans[2][2] = m[0][2]*m2[0][2] + m[1][2]*m2[1][2] + m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply mat1 times transpose of mat2
------------------------------------------------------------------------- */

void MathExtra::times3_transpose(const double m[3][3], const double m2[3][3],
                                 double ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0] + m[0][1]*m2[0][1] + m[0][2]*m2[0][2];
  ans[0][1] = m[0][0]*m2[1][0] + m[0][1]*m2[1][1] + m[0][2]*m2[1][2];
  ans[0][2] = m[0][0]*m2[2][0] + m[0][1]*m2[2][1] + m[0][2]*m2[2][2];
  ans[1][0] = m[1][0]*m2[0][0] + m[1][1]*m2[0][1] + m[1][2]*m2[0][2];
  ans[1][1] = m[1][0]*m2[1][0] + m[1][1]*m2[1][1] + m[1][2]*m2[1][2];
  ans[1][2] = m[1][0]*m2[2][0] + m[1][1]*m2[2][1] + m[1][2]*m2[2][2];
  ans[2][0] = m[2][0]*m2[0][0] + m[2][1]*m2[0][1] + m[2][2]*m2[0][2];
  ans[2][1] = m[2][0]*m2[1][0] + m[2][1]*m2[1][1] + m[2][2]*m2[1][2];
  ans[2][2] = m[2][0]*m2[2][0] + m[2][1]*m2[2][1] + m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   invert a matrix
   does NOT checks for singular or badly scaled matrix
------------------------------------------------------------------------- */

void MathExtra::invert3(const double m[3][3], double ans[3][3])
{
  double den = m[0][0]*m[1][1]*m[2][2]-m[0][0]*m[1][2]*m[2][1];
  den += -m[1][0]*m[0][1]*m[2][2]+m[1][0]*m[0][2]*m[2][1];
  den += m[2][0]*m[0][1]*m[1][2]-m[2][0]*m[0][2]*m[1][1];

  ans[0][0] = (m[1][1]*m[2][2]-m[1][2]*m[2][1]) / den;
  ans[0][1] = -(m[0][1]*m[2][2]-m[0][2]*m[2][1]) / den;
  ans[0][2] = (m[0][1]*m[1][2]-m[0][2]*m[1][1]) / den;
  ans[1][0] = -(m[1][0]*m[2][2]-m[1][2]*m[2][0]) / den;
  ans[1][1] = (m[0][0]*m[2][2]-m[0][2]*m[2][0]) / den;
  ans[1][2] = -(m[0][0]*m[1][2]-m[0][2]*m[1][0]) / den;
  ans[2][0] = (m[1][0]*m[2][1]-m[1][1]*m[2][0]) / den;
  ans[2][1] = -(m[0][0]*m[2][1]-m[0][1]*m[2][0]) / den;
  ans[2][2] = (m[0][0]*m[1][1]-m[0][1]*m[1][0]) / den;
}

/* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */

void MathExtra::matvec(const double m[3][3], const double *v, double *ans)
{
  ans[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
  ans[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
  ans[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
}

/* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */

void MathExtra::matvec(const double *ex, const double *ey, const double *ez,
                       const double *v, double *ans)
{
  ans[0] = ex[0]*v[0] + ey[0]*v[1] + ez[0]*v[2];
  ans[1] = ex[1]*v[0] + ey[1]*v[1] + ez[1]*v[2];
  ans[2] = ex[2]*v[0] + ey[2]*v[1] + ez[2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

void MathExtra::transpose_matvec(const double m[3][3], const double *v,
                                 double *ans)
{
  ans[0] = m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2];
  ans[1] = m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2];
  ans[2] = m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

void MathExtra::transpose_matvec(const double *ex, const double *ey,
                                 const double *ez, const double *v,
                                 double *ans)
{
  ans[0] = ex[0]*v[0] + ex[1]*v[1] + ex[2]*v[2];
  ans[1] = ey[0]*v[0] + ey[1]*v[1] + ey[2]*v[2];
  ans[2] = ez[0]*v[0] + ez[1]*v[1] + ez[2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times diagonal matrix
------------------------------------------------------------------------- */

void MathExtra::transpose_diag3(const double m[3][3], const double *d,
                                double ans[3][3])
{
  ans[0][0] = m[0][0]*d[0];
  ans[0][1] = m[1][0]*d[1];
  ans[0][2] = m[2][0]*d[2];
  ans[1][0] = m[0][1]*d[0];
  ans[1][1] = m[1][1]*d[1];
  ans[1][2] = m[2][1]*d[2];
  ans[2][0] = m[0][2]*d[0];
  ans[2][1] = m[1][2]*d[1];
  ans[2][2] = m[2][2]*d[2];
}

/* ----------------------------------------------------------------------
   row vector times matrix
------------------------------------------------------------------------- */

void MathExtra::vecmat(const double *v, const double m[3][3], double *ans)
{
  ans[0] = v[0]*m[0][0] + v[1]*m[1][0] + v[2]*m[2][0];
  ans[1] = v[0]*m[0][1] + v[1]*m[1][1] + v[2]*m[2][1];
  ans[2] = v[0]*m[0][2] + v[1]*m[1][2] + v[2]*m[2][2];
}

/* ----------------------------------------------------------------------
   matrix times scalar, in place
------------------------------------------------------------------------- */

void MathExtra::scalar_times3(const double f, double m[3][3])
{
  m[0][0] *= f; m[0][1] *= f; m[0][2] *= f;
  m[1][0] *= f; m[1][1] *= f; m[1][2] *= f;
  m[2][0] *= f; m[2][1] *= f; m[2][2] *= f;
}

/* ----------------------------------------------------------------------
   compute quaternion from axis-angle rotation
   v MUST be a unit vector
------------------------------------------------------------------------- */

void MathExtra::axisangle_to_quat(const double *v, const double angle,
                                  double *quat)
{
  double halfa = 0.5*angle;
  double sina = sin(halfa);
  quat[0] = cos(halfa);
  quat[1] = v[0]*sina;
  quat[2] = v[1]*sina;
  quat[3] = v[2]*sina;
}

/* ----------------------------------------------------------------------
   compute rotation matrix from quaternion
   quat = [w i j k]
------------------------------------------------------------------------- */

void MathExtra::quat_to_mat(const double *quat, double mat[3][3])
{
  double w2 = quat[0]*quat[0];
  double i2 = quat[1]*quat[1];
  double j2 = quat[2]*quat[2];
  double k2 = quat[3]*quat[3];
  double twoij = 2.0*quat[1]*quat[2];
  double twoik = 2.0*quat[1]*quat[3];
  double twojk = 2.0*quat[2]*quat[3];
  double twoiw = 2.0*quat[1]*quat[0];
  double twojw = 2.0*quat[2]*quat[0];
  double twokw = 2.0*quat[3]*quat[0];

  mat[0][0] = w2+i2-j2-k2;
  mat[0][1] = twoij-twokw;
  mat[0][2] = twojw+twoik;

  mat[1][0] = twoij+twokw;
  mat[1][1] = w2-i2+j2-k2;
  mat[1][2] = twojk-twoiw;

  mat[2][0] = twoik-twojw;
  mat[2][1] = twojk+twoiw;
  mat[2][2] = w2-i2-j2+k2;
}

#endif
