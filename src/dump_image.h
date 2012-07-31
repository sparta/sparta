/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   www.sandia.gov/sparta.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(image,DumpImage)

#else

#ifndef SPARTA_DUMP_IMAGE_H
#define SPARTA_DUMP_IMAGE_H

#include "dump_molecule.h"

namespace SPARTA_NS {

class DumpImage : public DumpMolecule {
 public:
  DumpImage(class SPARTA *, int, char**);
  ~DumpImage();
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);

 private:
  int filetype;
  int acolor,adiam;                // what determines color/diam of atoms
  double adiamvalue;               // atom diameter value
  int atomflag;                    // 0/1 for draw atoms
  char *thetastr,*phistr;          // variables for view theta,phi
  int thetavar,phivar;             // index to theta,phi vars
  int cflag;                       // static/dynamic box center
  double cx,cy,cz;                 // fractional box center
  char *cxstr,*cystr,*czstr;       // variables for box center
  int cxvar,cyvar,czvar;           // index to box center vars
  char *upxstr,*upystr,*upzstr;    // view up vector variables
  int upxvar,upyvar,upzvar;        // index to up vector vars
  char *zoomstr,*perspstr;         // view zoom and perspective variables
  int zoomvar,perspvar;            // index to zoom,persp vars
  int boxflag,axesflag;            // 0/1 for draw box and axes
  double boxdiam,axeslen,axesdiam; // params for drawing box and axes

  int viewflag;                    // overall view is static or dynamic

  double *diamtype,*diamelement;           // per-type diameters
  double **colortype,**colorelement;       // per-type colors

  class Image *image;              // class that renders each image

  double **bufcopy;                // buffer for communicating bond/atom info
  int maxbufcopy;

  void init_style();
  int modify_param(int, char **);
  void write();

  void box_center();
  void view_params();
  void box_bounds();

  void create_image();
};

}

#endif
#endif
