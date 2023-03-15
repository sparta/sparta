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

#include "math.h"
#include "ctype.h"
#include "stdlib.h"
#include "string.h"
#include "dump_image.h"
#include "image.h"
#include "domain.h"
#include "comm.h"
#include "region.h"
#include "particle.h"
#include "grid.h"
#include "surf.h"
#include "input.h"
#include "variable.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "math_extra.h"
#include "math_const.h"
#include "error.h"
#include "memory.h"

#ifdef SPARTA_JPEG
#include "jpeglib.h"
#endif

//#define RCB_DEBUG 1     // un-comment to include RCB proc boxes in image

using namespace SPARTA_NS;
using namespace MathConst;

#define BIG 1.0e20
#define INVOKED_PER_GRID 16
#define INVOKED_PER_SURF 32

enum{PPM,JPG,PNG};        // multiple files
enum{NUMERIC,TYPE,PROC,ATTRIBUTE,ONE};
enum{STATIC,DYNAMIC};
enum{COMPUTE,FIX,VARIABLE};
enum{PARTICLE,GRID,SURF,XPLANE,YPLANE,ZPLANE};

/* ---------------------------------------------------------------------- */

DumpImage::DumpImage(SPARTA *sparta, int narg, char **arg) :
  DumpParticle(sparta, narg, arg)
{
  if (binary || multiproc) error->all(FLERR,"Invalid dump image filename");

  // set filetype based on filename suffix

  int n = strlen(filename);
  if (strlen(filename) > 4 && strcmp(&filename[n-4],".jpg") == 0)
    filetype = JPG;
  else if (strlen(filename) > 4 && strcmp(&filename[n-4],".JPG") == 0)
    filetype = JPG;
  else if (strlen(filename) > 5 && strcmp(&filename[n-5],".jpeg") == 0)
    filetype = JPG;
  else if (strlen(filename) > 5 && strcmp(&filename[n-5],".JPEG") == 0)
    filetype = JPG;
  else if (strlen(filename) > 4 && strcmp(&filename[n-4],".png") == 0)
    filetype = PNG;
  else if (strlen(filename) > 4 && strcmp(&filename[n-4],".PNG") == 0)
    filetype = PNG;
  else filetype = PPM;

#ifndef SPARTA_JPEG
  if (filetype == JPG)
    error->all(FLERR,"Support for writing images in JPEG format not included");
#endif
#ifndef SPARTA_PNG
  if (filetype == PNG)
    error->all(FLERR,"Support for writing images in PNG format not included");
#endif

  // particle color,diameter settings

  if (nfield != 2) error->all(FLERR,"Illegal dump image command");

  pcolor = ATTRIBUTE;
  if (strcmp(arg[5],"type") == 0) pcolor = TYPE;
  else if (strcmp(arg[5],"proc") == 0) pcolor = PROC;

  pdiam = ATTRIBUTE;
  if (strcmp(arg[6],"type") == 0) pdiam = TYPE;

  // create Image class with 6 colormaps
  // colormaps for particles, grid, surf, grid xyz planes
  // change defaults for 2d

  image = new Image(sparta,6);
  if (domain->dimension == 2) {
    image->theta = 0.0;
    image->phi = 0.0;
    image->up[0] = 0.0; image->up[1] = 1.0; image->up[2] = 0.0;
  }

  boxcolor = image->color2rgb("yellow");
  surfcolorone = image->color2rgb("gray");
  glinecolor = image->color2rgb("white");
  slinecolor = image->color2rgb("white");

  // set defaults for optional args

  particleflag = 1;
  gridflag = 0;
  gridxflag = gridyflag = gridzflag = 0;
  surfflag = 0;
  thetastr = phistr = NULL;
  cflag = STATIC;
  cx = cy = cz = 0.5;
  cxstr = cystr = czstr = NULL;
  upxstr = upystr = upzstr = NULL;
  zoomstr = NULL;
  perspstr = NULL;
  boxflag = 1;
  boxdiam = 0.02;
  glineflag = 0;
  slineflag = 0;
  axesflag = 0;

  idgrid = idgridx = idgridy = idgridz = idsurf = NULL;

  // parse optional args

  int iarg = ioptional;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pdiam") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      pdiam = NUMERIC;
      pdiamvalue = atof(arg[iarg+1]);
      if (pdiamvalue <= 0.0) error->all(FLERR,"Illegal dump image command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"particle") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"yes") == 0) particleflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) particleflag = 0;
      else error->all(FLERR,"Illegal dump image command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"grid") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      gridflag = 1;
      if (strcmp(arg[iarg+1],"proc") == 0) gcolor = PROC;
      else {
        gcolor = ATTRIBUTE;
        if (strncmp(arg[iarg+1],"c_",2) && strncmp(arg[iarg+1],"f_",2) &&
            strncmp(arg[iarg+1],"v_",2))
          error->all(FLERR,"Illegal dump image command");
        if (arg[iarg+1][0] == 'c') gridwhich = COMPUTE;
        else if (arg[iarg+1][0] == 'f') gridwhich = FIX;
        else if (arg[iarg+1][0] == 'v') gridwhich = VARIABLE;

        int n = strlen(arg[iarg+1]);
        char *suffix = new char[n];
        strcpy(suffix,&arg[iarg+1][2]);

        char *ptr = strchr(suffix,'[');
        if (ptr) {
          if (suffix[strlen(suffix)-1] != ']')
            error->all(FLERR,"Illegal fix ave/grid command");
          gridcol = atoi(ptr+1);
          *ptr = '\0';
        } else gridcol = 0;
        n = strlen(suffix) + 1;
        idgrid = new char[n];
        strcpy(idgrid,suffix);
        delete [] suffix;
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"gridx") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      gridxflag = 1;
      gridxcoord = atof(arg[iarg+1]);
      if (strcmp(arg[iarg+2],"proc") == 0) gxcolor = PROC;
      else {
        gxcolor = ATTRIBUTE;
        if (strncmp(arg[iarg+2],"c_",2) && strncmp(arg[iarg+2],"f_",2) &&
            strncmp(arg[iarg+2],"v_",2))
          error->all(FLERR,"Illegal dump image command");
        if (arg[iarg+2][0] == 'c') gridxwhich = COMPUTE;
        else if (arg[iarg+2][0] == 'f') gridxwhich = FIX;
        else if (arg[iarg+2][0] == 'v') gridxwhich = VARIABLE;

        int n = strlen(arg[iarg+2]);
        char *suffix = new char[n];
        strcpy(suffix,&arg[iarg+2][2]);

        char *ptr = strchr(suffix,'[');
        if (ptr) {
          if (suffix[strlen(suffix)-1] != ']')
            error->all(FLERR,"Illegal fix ave/grid command");
          gridxcol = atoi(ptr+1);
          *ptr = '\0';
        } else gridxcol = 0;
        n = strlen(suffix) + 1;
        idgridx = new char[n];
        strcpy(idgridx,suffix);
        delete [] suffix;
      }
      iarg += 3;

    } else if (strcmp(arg[iarg],"gridy") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      gridyflag = 1;
      gridycoord = atof(arg[iarg+1]);
      if (strcmp(arg[iarg+2],"proc") == 0) gycolor = PROC;
      else {
        gycolor = ATTRIBUTE;
        if (strncmp(arg[iarg+2],"c_",2) && strncmp(arg[iarg+2],"f_",2) &&
            strncmp(arg[iarg+2],"v_",2))
          error->all(FLERR,"Illegal dump image command");
        if (arg[iarg+2][0] == 'c') gridywhich = COMPUTE;
        else if (arg[iarg+2][0] == 'f') gridywhich = FIX;
        else if (arg[iarg+2][0] == 'v') gridywhich = VARIABLE;

        int n = strlen(arg[iarg+2]);
        char *suffix = new char[n];
        strcpy(suffix,&arg[iarg+2][2]);

        char *ptr = strchr(suffix,'[');
        if (ptr) {
          if (suffix[strlen(suffix)-1] != ']')
            error->all(FLERR,"Illegal fix ave/grid command");
          gridycol = atoi(ptr+1);
          *ptr = '\0';
        } else gridycol = 0;
        n = strlen(suffix) + 1;
        idgridy = new char[n];
        strcpy(idgridy,suffix);
        delete [] suffix;
      }
      iarg += 3;

    } else if (strcmp(arg[iarg],"gridz") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      gridzflag = 1;
      gridzcoord = atof(arg[iarg+1]);
      if (strcmp(arg[iarg+2],"proc") == 0) gzcolor = PROC;
      else {
        gzcolor = ATTRIBUTE;
        if (strncmp(arg[iarg+2],"c_",2) && strncmp(arg[iarg+2],"f_",2) &&
            strncmp(arg[iarg+2],"v_",2))
          error->all(FLERR,"Illegal dump image command");
        if (arg[iarg+2][0] == 'c') gridzwhich = COMPUTE;
        else if (arg[iarg+2][0] == 'f') gridzwhich = FIX;
        else if (arg[iarg+2][0] == 'v') gridzwhich = VARIABLE;

        int n = strlen(arg[iarg+2]);
        char *suffix = new char[n];
        strcpy(suffix,&arg[iarg+2][2]);

        char *ptr = strchr(suffix,'[');
        if (ptr) {
          if (suffix[strlen(suffix)-1] != ']')
            error->all(FLERR,"Illegal fix ave/grid command");
          gridzcol = atoi(ptr+1);
          *ptr = '\0';
        } else gridzcol = 0;
        n = strlen(suffix) + 1;
        idgridz = new char[n];
        strcpy(idgridz,suffix);
        delete [] suffix;
      }
      iarg += 3;

    } else if (strcmp(arg[iarg],"surf") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      surfflag = 1;
      if (strcmp(arg[iarg+1],"one") == 0) scolor = ONE;
      else if (strcmp(arg[iarg+1],"proc") == 0) scolor = PROC;
      else {
        if (surf->implicit)
          error->all(FLERR,"Cannot use dump image surf with "
                     "compute/fix/variable for implicit surfs");
        scolor = ATTRIBUTE;
        if (strncmp(arg[iarg+1],"c_",2) && strncmp(arg[iarg+1],"f_",2) &&
            strncmp(arg[iarg+1],"v_",2))
          error->all(FLERR,"Illegal dump image command");
        if (arg[iarg+1][0] == 'c') surfwhich = COMPUTE;
        else if (arg[iarg+1][0] == 'f') surfwhich = FIX;
        else if (arg[iarg+1][0] == 'v') surfwhich = VARIABLE;

        int n = strlen(arg[iarg+1]);
        char *suffix = new char[n];
        strcpy(suffix,&arg[iarg+1][2]);

        char *ptr = strchr(suffix,'[');
        if (ptr) {
          if (suffix[strlen(suffix)-1] != ']')
            error->all(FLERR,"Illegal fix ave/grid command");
          surfcol = atoi(ptr+1);
          *ptr = '\0';
        } else surfcol = 0;
        n = strlen(suffix) + 1;
        idsurf = new char[n];
        strcpy(idsurf,suffix);
        delete [] suffix;
      }
      sdiamvalue = atof(arg[iarg+2]);
      iarg += 3;

    } else if (strcmp(arg[iarg],"size") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      int width = atoi(arg[iarg+1]);
      int height = atoi(arg[iarg+2]);
      if (width <= 0 || height <= 0)
        error->all(FLERR,"Illegal dump image command");
      image->width = width;
      image->height = height;
      iarg += 3;

    } else if (strcmp(arg[iarg],"view") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        thetastr = new char[n];
        strcpy(thetastr,&arg[iarg+1][2]);
      } else {
        double theta = atof(arg[iarg+1]);
        if (theta < 0.0 || theta > 180.0)
          error->all(FLERR,"Invalid dump image theta value");
        theta *= MY_PI/180.0;
        image->theta = theta;
      }
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        phistr = new char[n];
        strcpy(phistr,&arg[iarg+2][2]);
      } else {
        double phi = atof(arg[iarg+2]);
        phi *= MY_PI/180.0;
        image->phi = phi;
      }
      iarg += 3;

    } else if (strcmp(arg[iarg],"center") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"s") == 0) cflag = STATIC;
      else if (strcmp(arg[iarg+1],"d") == 0) cflag = DYNAMIC;
      else error->all(FLERR,"Illegal dump image command");
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        cxstr = new char[n];
        strcpy(cxstr,&arg[iarg+2][2]);
        cflag = DYNAMIC;
      } else cx = atof(arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        cystr = new char[n];
        strcpy(cystr,&arg[iarg+3][2]);
        cflag = DYNAMIC;
      } else cy = atof(arg[iarg+3]);
      if (strstr(arg[iarg+4],"v_") == arg[iarg+4]) {
        int n = strlen(&arg[iarg+4][2]) + 1;
        czstr = new char[n];
        strcpy(czstr,&arg[iarg+4][2]);
        cflag = DYNAMIC;
      } else cz = atof(arg[iarg+4]);
      iarg += 5;

    } else if (strcmp(arg[iarg],"up") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal dump image command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        upxstr = new char[n];
        strcpy(upxstr,&arg[iarg+1][2]);
      } else image->up[0] = atof(arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        upystr = new char[n];
        strcpy(upystr,&arg[iarg+2][2]);
      } else image->up[1] = atof(arg[iarg+1]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        upzstr = new char[n];
        strcpy(upzstr,&arg[iarg+3][2]);
      } else image->up[2] = atof(arg[iarg+3]);
      iarg += 4;

    } else if (strcmp(arg[iarg],"zoom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        zoomstr = new char[n];
        strcpy(zoomstr,&arg[iarg+1][2]);
      } else {
        double zoom = atof(arg[iarg+1]);
        if (zoom <= 0.0) error->all(FLERR,"Illegal dump image command");
        image->zoom = zoom;
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"persp") == 0) {
      error->all(FLERR,"Dump image persp option is not yet supported");
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        perspstr = new char[n];
        strcpy(perspstr,&arg[iarg+1][2]);
      } else {
        double persp = atof(arg[iarg+1]);
        if (persp < 0.0) error->all(FLERR,"Illegal dump image command");
        image->persp = persp;
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"box") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"yes") == 0) boxflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) boxflag = 0;
      else error->all(FLERR,"Illegal dump image command");
      boxdiam = atof(arg[iarg+2]);
      if (boxdiam < 0.0) error->all(FLERR,"Illegal dump image command");
      iarg += 3;

    } else if (strcmp(arg[iarg],"gline") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"yes") == 0) glineflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) glineflag = 0;
      else error->all(FLERR,"Illegal dump image command");
      glinediam = atof(arg[iarg+2]);
      if (glinediam < 0.0) error->all(FLERR,"Illegal dump image command");
      iarg += 3;

    } else if (strcmp(arg[iarg],"sline") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"yes") == 0) slineflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) slineflag = 0;
      else error->all(FLERR,"Illegal dump image command");
      slinediam = atof(arg[iarg+2]);
      if (slinediam < 0.0) error->all(FLERR,"Illegal dump image command");
      iarg += 3;

    } else if (strcmp(arg[iarg],"axes") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"yes") == 0) axesflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) axesflag = 0;
      else error->all(FLERR,"Illegal dump image command");
      axeslen = atof(arg[iarg+2]);
      axesdiam = atof(arg[iarg+3]);
      if (axeslen < 0.0 || axesdiam < 0.0)
        error->all(FLERR,"Illegal dump image command");
      iarg += 4;

    } else if (strcmp(arg[iarg],"shiny") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      double shiny = atof(arg[iarg+1]);
      if (shiny < 0.0 || shiny > 1.0)
        error->all(FLERR,"Illegal dump image command");
      image->shiny = shiny;
      iarg += 2;

    } else if (strcmp(arg[iarg],"ssao") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"yes") == 0) image->ssao = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) image->ssao = 0;
      else error->all(FLERR,"Illegal dump image command");
      int seed = atoi(arg[iarg+2]);
      if (seed <= 0) error->all(FLERR,"Illegal dump image command");
      image->seed = seed;
      double ssaoint = atof(arg[iarg+3]);
      if (ssaoint < 0.0 || ssaoint > 1.0)
        error->all(FLERR,"Illegal dump image command");
      image->ssaoint = ssaoint;
      iarg += 4;

    } else error->all(FLERR,"Illegal dump image command");
  }

  // error check

  if (gridflag && (gridxflag || gridyflag || gridzflag))
      error->all(FLERR,"Dump image cannot use grid and gridx/gridy/gridz");

  // allocate image buffer now that image size is known

  image->buffers();

  // additional defaults for dump_modify options

  pcolortype = new double*[ntypes+1];
  pdiamtype = new double[ntypes+1];

  for (int i = 1; i <= ntypes; i++) {
    pdiamtype[i] = 1.0;
    if (i % 6 == 1) pcolortype[i] = image->color2rgb("red");
    else if (i % 6 == 2) pcolortype[i] = image->color2rgb("green");
    else if (i % 6 == 3) pcolortype[i] = image->color2rgb("blue");
    else if (i % 6 == 4) pcolortype[i] = image->color2rgb("yellow");
    else if (i % 6 == 5) pcolortype[i] = image->color2rgb("aqua");
    else if (i % 6 == 0) pcolortype[i] = image->color2rgb("purple");
  }

  if (me % 6 == 0) pcolorproc = image->color2rgb("red");
  else if (me % 6 == 1) pcolorproc = image->color2rgb("green");
  else if (me % 6 == 2) pcolorproc = image->color2rgb("blue");
  else if (me % 6 == 3) pcolorproc = image->color2rgb("yellow");
  else if (me % 6 == 4) pcolorproc = image->color2rgb("aqua");
  else if (me % 6 == 5) pcolorproc = image->color2rgb("purple");

  if (me % 6 == 0) gcolorproc = image->color2rgb("red");
  else if (me % 6 == 1) gcolorproc = image->color2rgb("green");
  else if (me % 6 == 2) gcolorproc = image->color2rgb("blue");
  else if (me % 6 == 3) gcolorproc = image->color2rgb("yellow");
  else if (me % 6 == 4) gcolorproc = image->color2rgb("aqua");
  else if (me % 6 == 5) gcolorproc = image->color2rgb("purple");

  if (me % 6 == 0) scolorproc = image->color2rgb("red");
  else if (me % 6 == 1) scolorproc = image->color2rgb("green");
  else if (me % 6 == 2) scolorproc = image->color2rgb("blue");
  else if (me % 6 == 3) scolorproc = image->color2rgb("yellow");
  else if (me % 6 == 4) scolorproc = image->color2rgb("aqua");
  else if (me % 6 == 5) scolorproc = image->color2rgb("purple");

  // viewflag = DYNAMIC if any view parameter is dynamic

  viewflag = STATIC;
  if (thetastr || phistr || cflag == DYNAMIC ||
      upxstr || upystr || upzstr || zoomstr || perspstr) viewflag = DYNAMIC;

  if (cflag == STATIC) box_center();
  if (viewflag == STATIC) view_params();
}

/* ---------------------------------------------------------------------- */

DumpImage::~DumpImage()
{
  delete image;

  delete [] idgrid;
  delete [] idgridx;
  delete [] idgridy;
  delete [] idgridz;
  delete [] idsurf;

  delete [] pcolortype;
  delete [] pdiamtype;
}

/* ---------------------------------------------------------------------- */

void DumpImage::init_style()
{
  if (multifile == 0)
    error->all(FLERR,"Dump image requires one snapshot per file");

  DumpParticle::init_style();

  // check variables

  if (thetastr) {
    thetavar = input->variable->find(thetastr);
    if (thetavar < 0)
      error->all(FLERR,"Variable name for dump image theta does not exist");
    if (!input->variable->equal_style(thetavar))
      error->all(FLERR,"Variable for dump image theta is invalid style");
  }
  if (phistr) {
    phivar = input->variable->find(phistr);
    if (phivar < 0)
      error->all(FLERR,"Variable name for dump image phi does not exist");
    if (!input->variable->equal_style(phivar))
      error->all(FLERR,"Variable for dump image phi is invalid style");
  }
  if (cxstr) {
    cxvar = input->variable->find(cxstr);
    if (cxvar < 0)
      error->all(FLERR,"Variable name for dump image center does not exist");
    if (!input->variable->equal_style(cxvar))
      error->all(FLERR,"Variable for dump image center is invalid style");
  }
  if (cystr) {
    cyvar = input->variable->find(cystr);
    if (cyvar < 0)
      error->all(FLERR,"Variable name for dump image center does not exist");
    if (!input->variable->equal_style(cyvar))
      error->all(FLERR,"Variable for dump image center is invalid style");
  }
  if (czstr) {
    czvar = input->variable->find(czstr);
    if (czvar < 0)
      error->all(FLERR,"Variable name for dump image center does not exist");
    if (!input->variable->equal_style(czvar))
      error->all(FLERR,"Variable for dump image center is invalid style");
  }
  if (upxstr) {
    upxvar = input->variable->find(upxstr);
    if (upxvar < 0)
      error->all(FLERR,"Variable name for dump image center does not exist");
    if (!input->variable->equal_style(upxvar))
      error->all(FLERR,"Variable for dump image center is invalid style");
  }
  if (upystr) {
    upyvar = input->variable->find(upystr);
    if (upyvar < 0)
      error->all(FLERR,"Variable name for dump image center does not exist");
    if (!input->variable->equal_style(upyvar))
      error->all(FLERR,"Variable for dump image center is invalid style");
  }
  if (upzstr) {
    upzvar = input->variable->find(upzstr);
    if (upzvar < 0)
      error->all(FLERR,"Variable name for dump image center does not exist");
    if (!input->variable->equal_style(upzvar))
      error->all(FLERR,"Variable for dump image center is invalid style");
  }
  if (zoomstr) {
    zoomvar = input->variable->find(zoomstr);
    if (zoomvar < 0)
      error->all(FLERR,"Variable name for dump image zoom does not exist");
    if (!input->variable->equal_style(zoomvar))
      error->all(FLERR,"Variable for dump image zoom is invalid style");
  }
  if (perspstr) {
    perspvar = input->variable->find(perspstr);
    if (perspvar < 0)
      error->all(FLERR,"Variable name for dump image persp does not exist");
    if (!input->variable->equal_style(perspvar))
      error->all(FLERR,"Variable for dump image persp is invalid style");
  }

  // setup grid/surf references to computes, fixes, variables
  // NOTE: need to add lookup for variables

  if (gridflag && gcolor == ATTRIBUTE) {
    if (gridwhich == COMPUTE) {
      gridindex = modify->find_compute(idgrid);
      if (gridindex < 0)
        error->all(FLERR,"Could not find dump image compute ID");
      Compute *compute = modify->compute[gridindex];
      if (!compute->per_grid_flag)
        error->all(FLERR,"Dump image compute is not a per-grid compute");
      if (gridcol == 0 && compute->size_per_grid_cols != 0)
        error->all(FLERR,"Dump image compute does not produce a vector");
      if (gridcol > 0 && compute->size_per_grid_cols < gridcol)
        error->all(FLERR,"Dump image compute does not have requested column");
    } else if (gridwhich == FIX) {
      gridindex = modify->find_fix(idgrid);
      if (gridindex < 0)
        error->all(FLERR,"Could not find dump image fix ID");
      Fix *fix = modify->fix[gridindex];
      if (!fix->per_grid_flag)
        error->all(FLERR,"Dump image fix does not produce per-grid values");
      if (gridcol == 0 && fix->size_per_grid_cols != 0)
        error->all(FLERR,"Dump image fix does not produce a vector");
      if (gridcol > 0 && fix->size_per_grid_cols < gridcol)
        error->all(FLERR,"Dump image fix does not have requested column");
      if (nevery % fix->per_grid_freq)
        error->all(FLERR,"Dump image and fix not computed at compatible times");
    }
  }

  if (gridxflag && gxcolor == ATTRIBUTE) {
    if (gridxwhich == COMPUTE) {
      gridxindex = modify->find_compute(idgridx);
      if (gridxindex < 0)
        error->all(FLERR,"Could not find dump image compute ID");
      Compute *compute = modify->compute[gridxindex];
      if (!compute->per_grid_flag)
        error->all(FLERR,"Dump image compute is not a per-grid compute");
      if (gridxcol == 0 && compute->size_per_grid_cols != 0)
        error->all(FLERR,"Dump image compute does not produce a vector");
      if (gridxcol > 0 && compute->size_per_grid_cols < gridxcol)
        error->all(FLERR,"Dump image compute does not have requested column");
    } else if (gridxwhich == FIX) {
      gridxindex = modify->find_fix(idgridx);
      if (gridxindex < 0)
        error->all(FLERR,"Could not find dump image fix ID");
      Fix *fix = modify->fix[gridxindex];
      if (!fix->per_grid_flag)
        error->all(FLERR,"Dump image fix does not produce per-grid values");
      if (gridxcol == 0 && fix->size_per_grid_cols != 0)
        error->all(FLERR,"Dump image fix does not produce a vector");
      if (gridxcol > 0 && fix->size_per_grid_cols < gridxcol)
        error->all(FLERR,"Dump image fix does not have requested column");
      if (nevery % fix->per_grid_freq)
        error->all(FLERR,"Dump image and fix not computed at compatible times");
    }
  }

  if (gridyflag && gycolor == ATTRIBUTE) {
    if (gridywhich == COMPUTE) {
      gridyindex = modify->find_compute(idgridy);
      if (gridyindex < 0)
        error->all(FLERR,"Could not find dump image compute ID");
      Compute *compute = modify->compute[gridyindex];
      if (!compute->per_grid_flag)
        error->all(FLERR,"Dump image compute is not a per-grid compute");
      if (gridycol == 0 && compute->size_per_grid_cols != 0)
        error->all(FLERR,"Dump image compute does not produce a vector");
      if (gridycol > 0 && compute->size_per_grid_cols < gridycol)
        error->all(FLERR,"Dump image compute does not have requested column");
    } else if (gridywhich == FIX) {
      gridyindex = modify->find_fix(idgridy);
      if (gridyindex < 0)
        error->all(FLERR,"Could not find dump image fix ID");
      Fix *fix = modify->fix[gridyindex];
      if (!fix->per_grid_flag)
        error->all(FLERR,"Dump image fix does not produce per-grid values");
      if (gridycol == 0 && fix->size_per_grid_cols != 0)
        error->all(FLERR,"Dump image fix does not produce a vector");
      if (gridycol > 0 && fix->size_per_grid_cols < gridycol)
        error->all(FLERR,"Dump image fix does not have requested column");
      if (nevery % fix->per_grid_freq)
        error->all(FLERR,"Dump image and fix not computed at compatible times");
    }
  }

  if (gridzflag && gzcolor == ATTRIBUTE) {
    if (gridzwhich == COMPUTE) {
      gridzindex = modify->find_compute(idgridz);
      if (gridzindex < 0)
        error->all(FLERR,"Could not find dump image compute ID");
      Compute *compute = modify->compute[gridzindex];
      if (!compute->per_grid_flag)
        error->all(FLERR,"Dump image compute is not a per-grid compute");
      if (gridzcol == 0 && compute->size_per_grid_cols != 0)
        error->all(FLERR,"Dump image compute does not produce a vector");
      if (gridzcol > 0 && compute->size_per_grid_cols < gridzcol)
        error->all(FLERR,"Dump image compute does not have requested column");
    } else if (gridzwhich == FIX) {
      gridzindex = modify->find_fix(idgridz);
      if (gridzindex < 0)
        error->all(FLERR,"Could not find dump image fix ID");
      Fix *fix = modify->fix[gridzindex];
      if (!fix->per_grid_flag)
        error->all(FLERR,"Dump image fix does not produce per-grid values");
      if (gridzcol == 0 && fix->size_per_grid_cols != 0)
        error->all(FLERR,"Dump image fix does not produce a vector");
      if (gridzcol > 0 && fix->size_per_grid_cols < gridzcol)
        error->all(FLERR,"Dump image fix does not have requested column");
      if (nevery % fix->per_grid_freq)
        error->all(FLERR,"Dump image and fix not computed at compatible times");
    }
  }

  if (surfflag && scolor == ATTRIBUTE) {
    if (surfwhich == COMPUTE) {
      surfindex = modify->find_compute(idsurf);
      if (surfindex < 0)
        error->all(FLERR,"Could not find dump image compute ID");
      Compute *compute = modify->compute[surfindex];
      if (!compute->per_surf_flag)
        error->all(FLERR,"Dump image compute is not a per-surf compute");
      if (surfcol == 0 && compute->size_per_surf_cols != 0)
        error->all(FLERR,"Dump image compute does not produce a vector");
      if (surfcol > 0 && compute->size_per_surf_cols < surfcol)
        error->all(FLERR,"Dump image compute does not have requested column");
    } else if (surfwhich == FIX) {
      surfindex = modify->find_fix(idsurf);
      if (surfindex < 0)
        error->all(FLERR,"Could not find dump image fix ID");
      Fix *fix = modify->fix[surfindex];
      if (!fix->per_surf_flag)
        error->all(FLERR,"Dump image fix does not produce per-surf values");
      if (surfcol == 0 && fix->size_per_surf_cols != 0)
        error->all(FLERR,"Dump image fix does not produce a vector");
      if (surfcol > 0 && fix->size_per_surf_cols < surfcol)
        error->all(FLERR,"Dump image fix does not have requested column");
      if (nevery % fix->per_surf_freq)
        error->all(FLERR,"Dump image and fix not computed at compatible times");
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpImage::write()
{
  // open new file

  openfile();

  // reset box center and view parameters if dynamic

  if (cflag == DYNAMIC) box_center();
  if (viewflag == DYNAMIC) view_params();

  // nme = # of particles this proc will contribute to dump

  nme = count();

  if (nme > maxbuf) {
    maxbuf = nme;
    memory->destroy(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
  }

  // pack buf with color,diameter

  pack();

  // set minmax color ranges if using dynamic color maps

  if (particleflag && pcolor == ATTRIBUTE && image->map_dynamic(PARTICLE)) {
    double two[2],twoall[2];
    double lo = BIG;
    double hi = -BIG;
    int m = 0;
    for (int i = 0; i < nchoose; i++) {
      lo = MIN(lo,buf[m]);
      hi = MAX(hi,buf[m]);
      m += size_one;
    }
    two[0] = -lo;
    two[1] = hi;
    MPI_Allreduce(two,twoall,2,MPI_DOUBLE,MPI_MAX,world);
    int flag = image->map_minmax(PARTICLE,-twoall[0],twoall[1]);
    if (flag) error->all(FLERR,"Invalid color map min/max values");
  }

  Compute *c;
  Fix *f;
  int ppgflag;

  if (gridflag && gcolor == ATTRIBUTE && image->map_dynamic(GRID)) {
    double value,two[2],twoall[2];
    double lo = BIG;
    double hi = -BIG;

    if (gridwhich == COMPUTE) {
      c = modify->compute[gridindex];
      if (!(c->invoked_flag & INVOKED_PER_GRID)) {
        c->compute_per_grid();
        c->invoked_flag |= INVOKED_PER_GRID;
      }
      ppgflag = 0;
      if (c->post_process_grid_flag) {
        c->post_process_grid(gridcol,1,NULL,NULL,NULL,1);
        ppgflag = 1;
      } else if (c->post_process_isurf_grid_flag) c->post_process_isurf_grid();
    } else if (gridwhich == FIX) {
      f = modify->fix[gridindex];
    }

    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;

    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      if (gridwhich == COMPUTE) {
        if (gridcol == 0 || ppgflag) value = c->vector_grid[icell];
        else value = c->array_grid[icell][gridcol-1];
      } else if (gridwhich == FIX) {
        if (gridcol == 0) value = f->vector_grid[icell];
        else value = f->array_grid[icell][gridcol-1];
      }
      lo = MIN(lo,value);
      hi = MAX(hi,value);
    }
    two[0] = -lo;
    two[1] = hi;
    MPI_Allreduce(two,twoall,2,MPI_DOUBLE,MPI_MAX,world);
    int flag = image->map_minmax(GRID,-twoall[0],twoall[1]);
    if (flag) error->all(FLERR,"Invalid color map min/max values");
  }

  if (gridxflag && gxcolor == ATTRIBUTE && image->map_dynamic(XPLANE)) {
    double value,two[2],twoall[2];
    double lo = BIG;
    double hi = -BIG;

    if (gridxwhich == COMPUTE) {
      c = modify->compute[gridxindex];
      if (!(c->invoked_flag & INVOKED_PER_GRID)) {
        c->compute_per_grid();
        c->invoked_flag |= INVOKED_PER_GRID;
      }
      ppgflag = 0;
      if (c->post_process_grid_flag) {
        c->post_process_grid(gridxcol,1,NULL,NULL,NULL,1);
        ppgflag = 1;
      } else if (c->post_process_isurf_grid_flag) c->post_process_isurf_grid();
    } else if (gridxwhich == FIX) {
      f = modify->fix[gridxindex];
    }

    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;

    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      if (gridxwhich == COMPUTE) {
        if (gridxcol == 0 || ppgflag) value = c->vector_grid[icell];
        else value = c->array_grid[icell][gridxcol-1];
      } else if (gridxwhich == FIX) {
        if (gridxcol == 0) value = f->vector_grid[icell];
        else value = f->array_grid[icell][gridxcol-1];
      }
      lo = MIN(lo,value);
      hi = MAX(hi,value);
    }
    two[0] = -lo;
    two[1] = hi;
    MPI_Allreduce(two,twoall,2,MPI_DOUBLE,MPI_MAX,world);
    int flag = image->map_minmax(XPLANE,-twoall[0],twoall[1]);
    if (flag) error->all(FLERR,"Invalid color map min/max values");
  }

  if (gridyflag && gycolor == ATTRIBUTE && image->map_dynamic(YPLANE)) {
    double value,two[2],twoall[2];
    double lo = BIG;
    double hi = -BIG;

    if (gridywhich == COMPUTE) {
      c = modify->compute[gridyindex];
      if (!(c->invoked_flag & INVOKED_PER_GRID)) {
        c->compute_per_grid();
        c->invoked_flag |= INVOKED_PER_GRID;
      }
      ppgflag = 0;
      if (c->post_process_grid_flag) {
        c->post_process_grid(gridycol,1,NULL,NULL,NULL,1);
        ppgflag = 1;
      } else if (c->post_process_isurf_grid_flag) c->post_process_isurf_grid();
    } else if (gridywhich == FIX) {
      f = modify->fix[gridyindex];
    }

    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;

    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      if (gridywhich == COMPUTE) {
        if (gridycol == 0 || ppgflag) value = c->vector_grid[icell];
        else value = c->array_grid[icell][gridycol-1];
      } else if (gridywhich == FIX) {
        if (gridycol == 0) value = f->vector_grid[icell];
        else value = f->array_grid[icell][gridycol-1];
      }
      lo = MIN(lo,value);
      hi = MAX(hi,value);
    }
    two[0] = -lo;
    two[1] = hi;
    MPI_Allreduce(two,twoall,2,MPI_DOUBLE,MPI_MAX,world);
    int flag = image->map_minmax(YPLANE,-twoall[0],twoall[1]);
    if (flag) error->all(FLERR,"Invalid color map min/max values");
  }

  if (gridzflag && gzcolor == ATTRIBUTE && image->map_dynamic(ZPLANE)) {
    double value,two[2],twoall[2];
    double lo = BIG;
    double hi = -BIG;

    if (gridzwhich == COMPUTE) {
      c = modify->compute[gridzindex];
      if (!(c->invoked_flag & INVOKED_PER_GRID)) {
        c->compute_per_grid();
        c->invoked_flag |= INVOKED_PER_GRID;
      }
      ppgflag = 0;
      if (c->post_process_grid_flag) {
        c->post_process_grid(gridzcol,1,NULL,NULL,NULL,1);
        ppgflag = 1;
      } else if (c->post_process_isurf_grid_flag) c->post_process_isurf_grid();
    } else if (gridzwhich == FIX) {
      f = modify->fix[gridzindex];
    }

    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;

    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      if (gridzwhich == COMPUTE) {
        if (gridzcol == 0 || ppgflag) value = c->vector_grid[icell];
        else value = c->array_grid[icell][gridzcol-1];
      } else if (gridzwhich == FIX) {
        if (gridzcol == 0) value = f->vector_grid[icell];
        else value = f->array_grid[icell][gridzcol-1];
      }
      lo = MIN(lo,value);
      hi = MAX(hi,value);
    }
    two[0] = -lo;
    two[1] = hi;
    MPI_Allreduce(two,twoall,2,MPI_DOUBLE,MPI_MAX,world);
    int flag = image->map_minmax(ZPLANE,-twoall[0],twoall[1]);
    if (flag) error->all(FLERR,"Invalid color map min/max values");
  }

  if (surfflag && scolor == ATTRIBUTE && image->map_dynamic(SURF)) {
    double value,two[2],twoall[2];
    double lo = BIG;
    double hi = -BIG;

    if (surfwhich == COMPUTE) {
      c = modify->compute[surfindex];
      if (!(c->invoked_flag & INVOKED_PER_GRID)) {
        c->compute_per_grid();
        c->invoked_flag |= INVOKED_PER_GRID;
      }
      c->post_process_surf();
    } else if (surfwhich == FIX) {
      f = modify->fix[surfindex];
    }

    // only computes/fixes for explicit surfs allowed for this setting

    int nsurf = surf->nown;

    for (int isurf = 0; isurf < nsurf; isurf++) {
      if (surfwhich == COMPUTE) {
        if (surfcol == 0) value = c->vector_surf[isurf];
        else value = c->array_surf[isurf][surfcol-1];
      } else if (surfwhich == FIX) {
        if (surfcol == 0) value = f->vector_surf[isurf];
        else value = f->array_surf[isurf][surfcol-1];
      }
      lo = MIN(lo,value);
      hi = MAX(hi,value);
    }
    two[0] = -lo;
    two[1] = hi;
    MPI_Allreduce(two,twoall,2,MPI_DOUBLE,MPI_MAX,world);
    int flag = image->map_minmax(SURF,-twoall[0],twoall[1]);
    if (flag) error->all(FLERR,"Invalid color map min/max values");
  }

  // create my portion of image for my particles, grid cells, surfs
  // then merge per-proc images

  image->clear();
  create_image();
  image->merge();

  // write image file

  if (me == 0) {
    if (filetype == JPG) image->write_JPG(fp);
    else if (filetype == PNG) image->write_PNG(fp);
    else image->write_PPM(fp);
    if (multifile) {
      fclose(fp);
      fp = NULL;
    }
  }
}

/* ----------------------------------------------------------------------
   simulation box bounds
------------------------------------------------------------------------- */

void DumpImage::box_bounds()
{
  boxxlo = domain->boxlo[0];
  boxxhi = domain->boxhi[0];
  boxylo = domain->boxlo[1];
  boxyhi = domain->boxhi[1];
  boxzlo = domain->boxlo[2];
  boxzhi = domain->boxhi[2];
}

/* ----------------------------------------------------------------------
   reset view parameters
   called once from constructor if view is STATIC
   called every snapshot from write() if view is DYNAMIC
------------------------------------------------------------------------- */

void DumpImage::box_center()
{
  box_bounds();

  if (cxstr) cx = input->variable->compute_equal(cxvar);
  if (cystr) cy = input->variable->compute_equal(cyvar);
  if (czstr) cz = input->variable->compute_equal(czvar);

  image->xctr = boxxlo + cx*(boxxhi-boxxlo);
  image->yctr = boxylo + cy*(boxyhi-boxylo);
  image->zctr = boxzlo + cz*(boxzhi-boxzlo);
}

/* ----------------------------------------------------------------------
   reset view parameters in Image class
   called once from constructor if view is STATIC
   called every snapshot from write() if view is DYNAMIC
------------------------------------------------------------------------- */

void DumpImage::view_params()
{
  // view direction theta and phi

  if (thetastr) {
    double theta = input->variable->compute_equal(thetavar);
    if (theta < 0.0 || theta > 180.0)
      error->all(FLERR,"Invalid dump image theta value");
    theta *= MY_PI/180.0;
    image->theta = theta;
  }

  if (phistr) {
    double phi = input->variable->compute_equal(phivar);
    phi *= MY_PI/180.0;
    image->phi = phi;
  }

  // up vector

  if (upxstr) image->up[0] = input->variable->compute_equal(upxvar);
  if (upystr) image->up[1] = input->variable->compute_equal(upyvar);
  if (upzstr) image->up[2] = input->variable->compute_equal(upzvar);

  // zoom and perspective

  if (zoomstr) image->zoom = input->variable->compute_equal(zoomvar);
  if (image->zoom <= 0.0) error->all(FLERR,"Invalid dump image zoom value");
  if (perspstr) image->persp = input->variable->compute_equal(perspvar);
  if (image->persp < 0.0) error->all(FLERR,"Invalid dump image persp value");

  // current simulation box bounds

  box_bounds();

  // remainder of view setup is internal to Image class

  image->view_params(boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi);
}

/* ----------------------------------------------------------------------
   create image for particles on this proc
   every pixel has depth
------------------------------------------------------------------------- */

void DumpImage::create_image()
{
  int i,j,m,itype;
  double diameter;
  double *color;

  // render my partiless
  // region is used as constraint by parent class

  if (particleflag) {
    Particle::OnePart *particles = particle->particles;

    m = 0;
    for (i = 0; i < nchoose; i++) {
      j = clist[i];

      if (pcolor == TYPE) {
        itype = static_cast<int> (buf[m]);
        color = pcolortype[itype];
      } else if (pcolor == PROC) {
        color = pcolorproc;
      } else if (pcolor == ATTRIBUTE) {
        color = image->map_value2color(PARTICLE,buf[m]);
      }

      if (pdiam == NUMERIC) {
        diameter = pdiamvalue;
      } else if (pdiam == TYPE) {
        itype = static_cast<int> (buf[m+1]);
        diameter = pdiamtype[itype];
      } else if (pdiam == ATTRIBUTE) {
        diameter = buf[m+1];
      }

      image->draw_sphere(particles[j].x,color,diameter);
      m += size_one;
    }
  }

  // render my grid cells
  // for 2d, draw as rectangle, so particles and grid lines show up
  // for 3d, draw as brick, so particles will be hidden inside
  // use region as constraint

  if (gridflag) {
    double value;
    double x[3],diam[3];
    double *lo,*hi;

    Compute *c;
    Fix *f;
    int ppgflag;

    if (gcolor == ATTRIBUTE && gridwhich == COMPUTE) {
      c = modify->compute[gridindex];
      if (!(c->invoked_flag & INVOKED_PER_GRID)) {
        c->compute_per_grid();
        c->invoked_flag |= INVOKED_PER_GRID;
      }
      ppgflag = 0;
      if (c->post_process_grid_flag) {
        c->post_process_grid(gridcol,1,NULL,NULL,NULL,1);
        ppgflag = 1;
      } else if (c->post_process_isurf_grid_flag) c->post_process_isurf_grid();
    } else if (gcolor == ATTRIBUTE && gridwhich == FIX) {
      f = modify->fix[gridindex];
    }

    Region *region = NULL;
    if (iregion >= 0) region = domain->regions[iregion];

    int dimension = domain->dimension;
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;

    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;

      if (gcolor == PROC) {
        color = gcolorproc;
      } else if (gcolor == ATTRIBUTE) {
        if (gridwhich == COMPUTE) {
          if (gridcol == 0 || ppgflag) value = c->vector_grid[icell];
          else value = c->array_grid[icell][gridcol-1];
        } else if (gridwhich == FIX) {
          if (gridcol == 0) value = f->vector_grid[icell];
          else value = f->array_grid[icell][gridcol-1];
        }
        color = image->map_value2color(GRID,value);
      }

      lo = cells[icell].lo;
      hi = cells[icell].hi;

      diam[0] = hi[0]-lo[0];
      diam[1] = hi[1]-lo[1];
      diam[2] = hi[2]-lo[2];
      x[0] = 0.5*(lo[0]+hi[0]);
      x[1] = 0.5*(lo[1]+hi[1]);
      x[2] = 0.5*(lo[2]+hi[2]);
      if (dimension == 2) diam[2] = x[2] = 0.0;

      if (region && !region->match(x)) continue;

      image->draw_brick(x,color,diam);
    }
  }

  // render my grid plane(s)
  // do not use region as constraint
  // NOTE: could add this if allow for specifying different regions

  if (gridxflag || gridyflag || gridzflag) {
    double value;
    double x[3],diam[3];
    double *lo,*hi;

    Compute *cx,*cy,*cz;
    Fix *fx,*fy,*fz;
    int ppgflagx,ppgflagy,ppgflagz;

    if (gridxflag && gxcolor == ATTRIBUTE) {
      if (gridxwhich == COMPUTE) {
        cx = modify->compute[gridxindex];
        if (!(cx->invoked_flag & INVOKED_PER_GRID)) {
          cx->compute_per_grid();
          cx->invoked_flag |= INVOKED_PER_GRID;
        }
        ppgflagx = 0;
        if (cx->post_process_grid_flag) {
          cx->post_process_grid(gridxcol,1,NULL,NULL,NULL,1);
          ppgflagx = 1;
        } else if (cx->post_process_isurf_grid_flag)
          cx->post_process_isurf_grid();
      } else if (gridxwhich == FIX) {
        fx = modify->fix[gridxindex];
      }
    }

    if (gridyflag && gycolor == ATTRIBUTE) {
      if (gridywhich == COMPUTE) {
        cy = modify->compute[gridyindex];
        if (!(cy->invoked_flag & INVOKED_PER_GRID)) {
          cy->compute_per_grid();
          cy->invoked_flag |= INVOKED_PER_GRID;
        }
        ppgflagy = 0;
        if (cy->post_process_grid_flag) {
          cy->post_process_grid(gridycol,1,NULL,NULL,NULL,1);
          ppgflagy = 1;
        } else if (cy->post_process_isurf_grid_flag)
          cy->post_process_isurf_grid();
      } else if (gridywhich == FIX) {
        fy = modify->fix[gridyindex];
      }
    }

    if (gridzflag && gzcolor == ATTRIBUTE) {
      if (gridzwhich == COMPUTE) {
        cz = modify->compute[gridzindex];
        if (!(cz->invoked_flag & INVOKED_PER_GRID)) {
          cz->compute_per_grid();
          cz->invoked_flag |= INVOKED_PER_GRID;
        }
        ppgflagz = 0;
        if (cz->post_process_grid_flag) {
          cz->post_process_grid(gridzcol,1,NULL,NULL,NULL,1);
          ppgflagz = 1;
        } else if (cz->post_process_isurf_grid_flag)
          cz->post_process_isurf_grid();
      } else if (gridzwhich == FIX) {
        fz = modify->fix[gridzindex];
      }
    }

    Region *region = NULL;

    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;

    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;

      // draw as rectangle, so grid lines show up

      lo = cells[icell].lo;
      hi = cells[icell].hi;

      if (gridxflag) {
        if (gridxcoord >= lo[0] && gridxcoord < hi[0]) {
          if (gxcolor == PROC) {
            color = gcolorproc;
          } else if (gxcolor == ATTRIBUTE) {
            if (gridxwhich == COMPUTE) {
              if (gridxcol == 0 || ppgflagx) value = cx->vector_grid[icell];
              else value = cx->array_grid[icell][gridxcol-1];
            } else if (gridxwhich == FIX) {
              if (gridxcol == 0) value = fx->vector_grid[icell];
              else value = fx->array_grid[icell][gridxcol-1];
            }
            color = image->map_value2color(XPLANE,value);
          }
          diam[0] = 0.0;
          diam[1] = hi[1]-lo[1];
          diam[2] = hi[2]-lo[2];
          x[0] = gridxcoord;
          x[1] = 0.5*(lo[1]+hi[1]);
          x[2] = 0.5*(lo[2]+hi[2]);
          if (!region || region->match(x)) image->draw_brick(x,color,diam);
        }
      }

      if (gridyflag) {
        if (gridycoord >= lo[1] && gridycoord < hi[1]) {
          if (gycolor == PROC) {
            color = gcolorproc;
          } else if (gycolor == ATTRIBUTE) {
            if (gridywhich == COMPUTE) {
              if (gridycol == 0) value = cy->vector_grid[icell];
              else value = cy->array_grid[icell][gridycol-1];
            } else if (gridywhich == FIX) {
              if (gridycol == 0) value = fy->vector_grid[icell];
              else value = fy->array_grid[icell][gridycol-1];
            }
            color = image->map_value2color(YPLANE,value);
          }
          diam[0] = hi[0]-lo[0];
          diam[1] = 0.0;
          diam[2] = hi[2]-lo[2];
          x[0] = 0.5*(lo[0]+hi[0]);
          x[1] = gridycoord;
          x[2] = 0.5*(lo[2]+hi[2]);
          if (!region || region->match(x)) image->draw_brick(x,color,diam);
        }
      }

      if (gridzflag) {
        if (gridzcoord >= lo[2] && gridzcoord < hi[2]) {
          if (gzcolor == PROC) {
            color = gcolorproc;
          } else if (gzcolor == ATTRIBUTE) {
            if (gridzwhich == COMPUTE) {
              if (gridzcol == 0) value = cz->vector_grid[icell];
              else value = cz->array_grid[icell][gridzcol-1];
            } else if (gridzwhich == FIX) {
              if (gridzcol == 0) value = fz->vector_grid[icell];
              else value = fz->array_grid[icell][gridzcol-1];
            }
            color = image->map_value2color(ZPLANE,value);
          }
          diam[0] = hi[0]-lo[0];
          diam[1] = hi[1]-lo[1];
          diam[2] = 0.0;
          x[0] = 0.5*(lo[0]+hi[0]);
          x[1] = 0.5*(lo[1]+hi[1]);
          x[2] = gridzcoord;
          if (!region || region->match(x)) image->draw_brick(x,color,diam);
        }
      }
    }
  }

  // render outline of my grid plane(s)
  // do not use region as constraint
  // NOTE: could add this if allow for specifying different regions

  if (glineflag && (gridxflag || gridyflag || gridzflag)) {
    double x[3];
    double *lo,*hi;

    Region *region = NULL;
    //if (iregion >= 0) region = domain->regions[iregion];

    diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
    if (domain->dimension == 3) diameter = MIN(diameter,boxzhi-boxzlo);
    diameter *= glinediam;

    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;

    double box[8][3];

    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      lo = cells[icell].lo;
      hi = cells[icell].hi;

      if (gridxflag) {
        if (gridxcoord >= lo[0] && gridxcoord < hi[0]) {
          x[0] = gridxcoord;
          x[1] = 0.5*(lo[1]+hi[1]);
          x[2] = 0.5*(lo[2]+hi[2]);
          if (!region || region->match(x)) {
            box[0][0] = gridxcoord; box[0][1] = lo[1]; box[0][2] = lo[2];
            box[1][0] = gridxcoord; box[1][1] = hi[1]; box[1][2] = lo[2];
            box[2][0] = gridxcoord; box[2][1] = lo[1]; box[2][2] = hi[2];
            box[3][0] = gridxcoord; box[3][1] = hi[1]; box[3][2] = hi[2];
            image->draw_box2d(box,glinecolor,diameter);
          }
        }
      }
      if (gridyflag) {
        if (gridycoord >= lo[1] && gridycoord < hi[1]) {
          x[0] = 0.5*(lo[0]+hi[0]);
          x[1] = gridycoord;
          x[2] = 0.5*(lo[2]+hi[2]);
          if (!region || region->match(x)) {
            box[0][0] = lo[0]; box[0][1] = gridycoord; box[0][2] = lo[2];
            box[1][0] = hi[0]; box[1][1] = gridycoord; box[1][2] = lo[2];
            box[2][0] = lo[0]; box[2][1] = gridycoord; box[2][2] = hi[2];
            box[3][0] = hi[0]; box[3][1] = gridycoord; box[3][2] = hi[2];
            image->draw_box2d(box,glinecolor,diameter);
          }
        }
      }
      if (gridzflag) {
        if (gridzcoord >= lo[2] && gridzcoord < hi[2]) {
          x[0] = 0.5*(lo[0]+hi[0]);
          x[1] = 0.5*(lo[1]+hi[1]);
          x[2] = gridzcoord;
          if (!region || region->match(x)) {
            box[0][0] = lo[0]; box[0][1] = lo[1]; box[0][2] = gridzcoord;
            box[1][0] = hi[0]; box[1][1] = lo[1]; box[1][2] = gridzcoord;
            box[2][0] = lo[0]; box[2][1] = hi[1]; box[2][2] = gridzcoord;
            box[3][0] = hi[0]; box[3][1] = hi[1]; box[3][2] = gridzcoord;
            image->draw_box2d(box,glinecolor,diameter);
          }
        }
      }
    }
  }

  // render outline of my grid cells
  // use region as constraint

  else if (glineflag) {
    double x[3];
    double *lo,*hi;

    Region *region = NULL;
    if (iregion >= 0) region = domain->regions[iregion];

    diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
    if (domain->dimension == 3) diameter = MIN(diameter,boxzhi-boxzlo);
    diameter *= glinediam;

    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;

    double box[8][3];

    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      lo = cells[icell].lo;
      hi = cells[icell].hi;

      x[0] = 0.5*(lo[0]+hi[0]);
      x[1] = 0.5*(lo[1]+hi[1]);
      x[2] = 0.5*(lo[2]+hi[2]);
      if (region && !region->match(x)) continue;

      box[0][0] = lo[0]; box[0][1] = lo[1]; box[0][2] = lo[2];
      box[1][0] = hi[0]; box[1][1] = lo[1]; box[1][2] = lo[2];
      box[2][0] = lo[0]; box[2][1] = hi[1]; box[2][2] = lo[2];
      box[3][0] = hi[0]; box[3][1] = hi[1]; box[3][2] = lo[2];
      box[4][0] = lo[0]; box[4][1] = lo[1]; box[4][2] = hi[2];
      box[5][0] = hi[0]; box[5][1] = lo[1]; box[5][2] = hi[2];
      box[6][0] = lo[0]; box[6][1] = hi[1]; box[6][2] = hi[2];
      box[7][0] = hi[0]; box[7][1] = hi[1]; box[7][2] = hi[2];

      image->draw_box(box,glinecolor,diameter);
    }
  }

  // render my surf elements
  // do not use region as constraint

  if (surfflag) {
    double value;

    Compute *c;
    Fix *f;

    if (scolor == ATTRIBUTE && surfwhich == COMPUTE) {
      c = modify->compute[surfindex];
      if (!(c->invoked_flag & INVOKED_PER_GRID)) {
        c->compute_per_grid();
        c->invoked_flag |= INVOKED_PER_GRID;
      }
      c->post_process_surf();
    } else if (scolor == ATTRIBUTE && surfwhich == FIX) {
      f = modify->fix[surfindex];
    }

    diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
    if (domain->dimension == 3) diameter = MIN(diameter,boxzhi-boxzlo);
    diameter *= sdiamvalue;

    int me = comm->me;
    int nprocs = comm->nprocs;
    int dim = domain->dimension;
    int nsurf,strided;
    Surf::Line *lines;
    Surf::Tri *tris;

    if (!surf->distributed) {
      nsurf = surf->nown;
      strided = 1;
      lines = surf->lines;
      tris = surf->tris;
    } else if (!surf->implicit) {
      nsurf = surf->nown;
      strided = 0;
      lines = surf->mylines;
      tris = surf->mytris;
    } else {
      nsurf = surf->nlocal;
      strided = 0;
      lines = surf->lines;
      tris = surf->tris;
    }

    for (int isurf = 0; isurf < nsurf; isurf++) {
      if (scolor == ONE) {
        color = surfcolorone;
      } else if (scolor == PROC) {
        color = scolorproc;
      } else if (scolor == ATTRIBUTE) {
        if (surfwhich == COMPUTE) {
          if (surfcol == 0) value = c->vector_surf[isurf];
          else value = c->array_surf[isurf][surfcol-1];
        } else if (surfwhich == FIX) {
          if (surfcol == 0) value = f->vector_surf[isurf];
          else value = f->array_surf[isurf][surfcol-1];
        }
        color = image->map_value2color(SURF,value);
      }

      if (strided) m = me + isurf*nprocs;
      else m = isurf;

      if (dim == 2)
        image->draw_line(lines[m].p1,lines[m].p2,color,diameter);
      else
        image->draw_triangle(tris[m].p1,tris[m].p2,tris[m].p3,color);
    }
  }

  // render outline of my surf elements
  // do not use region as constraint

  if (slineflag) {
    diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
    if (domain->dimension == 3) diameter = MIN(diameter,boxzhi-boxzlo);
    diameter *= slinediam;

    int me = comm->me;
    int nprocs = comm->nprocs;
    int dim = domain->dimension;
    int nsurf,strided;
    Surf::Line *lines;
    Surf::Tri *tris;

    if (!surf->distributed) {
      nsurf = surf->nown;
      strided = 1;
      lines = surf->lines;
      tris = surf->tris;
    } else if (!surf->implicit) {
      nsurf = surf->nown;
      strided = 0;
      lines = surf->mylines;
      tris = surf->mytris;
    } else {
      nsurf = surf->nlocal;
      strided = 0;
      lines = surf->lines;
      tris = surf->tris;
    }

    if (dim == 2) {
      for (int isurf = 0; isurf < nsurf; isurf++) {
        if (strided) m = me + isurf*nprocs;
        else m = isurf;
        image->draw_line(lines[m].p1,lines[m].p2,slinecolor,diameter);
      }
    } else {
      for (int isurf = 0; isurf < nsurf; isurf++) {
        if (strided) m = me + isurf*nprocs;
        else m = isurf;
        image->draw_line(tris[m].p1,tris[m].p2,slinecolor,diameter);
        image->draw_line(tris[m].p2,tris[m].p3,slinecolor,diameter);
        image->draw_line(tris[m].p3,tris[m].p1,slinecolor,diameter);
      }
    }
  }

  // render outline of simulation box, orthogonal or triclinic

  if (boxflag) {
    diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
    if (domain->dimension == 3) diameter = MIN(diameter,boxzhi-boxzlo);
    diameter *= boxdiam;

    double box[8][3];
    box[0][0] = boxxlo; box[0][1] = boxylo; box[0][2] = boxzlo;
    box[1][0] = boxxhi; box[1][1] = boxylo; box[1][2] = boxzlo;
    box[2][0] = boxxlo; box[2][1] = boxyhi; box[2][2] = boxzlo;
    box[3][0] = boxxhi; box[3][1] = boxyhi; box[3][2] = boxzlo;
    box[4][0] = boxxlo; box[4][1] = boxylo; box[4][2] = boxzhi;
    box[5][0] = boxxhi; box[5][1] = boxylo; box[5][2] = boxzhi;
    box[6][0] = boxxlo; box[6][1] = boxyhi; box[6][2] = boxzhi;
    box[7][0] = boxxhi; box[7][1] = boxyhi; box[7][2] = boxzhi;

    image->draw_box(box,boxcolor,diameter);
  }

  // render XYZ axes in red/green/blue
  // offset by 10% of box size and scale by axeslen

  if (axesflag) {
    diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
    if (domain->dimension == 3) diameter = MIN(diameter,boxzhi-boxzlo);
    diameter *= axesdiam;

    double axes[4][3];
    axes[0][0] = boxxlo; axes[0][1] = boxylo; axes[0][2] = boxzlo;
    axes[1][0] = boxxhi; axes[1][1] = boxylo; axes[1][2] = boxzlo;
    axes[2][0] = boxxlo; axes[2][1] = boxyhi; axes[2][2] = boxzlo;
    axes[3][0] = boxxlo; axes[3][1] = boxylo; axes[3][2] = boxzhi;

    double offset = MAX(boxxhi-boxxlo,boxyhi-boxylo);
    if (domain->dimension == 3) offset = MAX(offset,boxzhi-boxzlo);
    offset *= 0.1;
    axes[0][0] -= offset; axes[0][1] -= offset; axes[0][2] -= offset;
    axes[1][0] -= offset; axes[1][1] -= offset; axes[1][2] -= offset;
    axes[2][0] -= offset; axes[2][1] -= offset; axes[2][2] -= offset;
    axes[3][0] -= offset; axes[3][1] -= offset; axes[3][2] -= offset;

    axes[1][0] = axes[0][0] + axeslen*(axes[1][0]-axes[0][0]);
    axes[1][1] = axes[0][1] + axeslen*(axes[1][1]-axes[0][1]);
    axes[1][2] = axes[0][2] + axeslen*(axes[1][2]-axes[0][2]);
    axes[2][0] = axes[0][0] + axeslen*(axes[2][0]-axes[0][0]);
    axes[2][1] = axes[0][1] + axeslen*(axes[2][1]-axes[0][1]);
    axes[2][2] = axes[0][2] + axeslen*(axes[2][2]-axes[0][2]);
    axes[3][0] = axes[0][0] + axeslen*(axes[3][0]-axes[0][0]);
    axes[3][1] = axes[0][1] + axeslen*(axes[3][1]-axes[0][1]);
    axes[3][2] = axes[0][2] + axeslen*(axes[3][2]-axes[0][2]);

    image->draw_axes(axes,diameter);
  }

  // DEBUG - render each proc's RCB box
  // need to add logic for 3d to this

#ifdef RCB_DEBUG

  diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
  if (domain->dimension == 3) diameter = MIN(diameter,boxzhi-boxzlo);
  diameter *= 0.5*boxdiam;

  double *rcblo = update->rcblo;
  double *rcbhi = update->rcbhi;

  double box[4][3];
  box[0][0] = rcblo[0]; box[0][1] = rcblo[1]; box[0][2] = boxzhi;
  box[1][0] = rcbhi[0]; box[1][1] = rcblo[1]; box[1][2] = boxzhi;
  box[2][0] = rcblo[0]; box[2][1] = rcbhi[1]; box[2][2] = boxzhi;
  box[3][0] = rcbhi[0]; box[3][1] = rcbhi[1]; box[3][2] = boxzhi;
  image->draw_box2d(box,boxcolor,diameter);

#endif

}

/* ---------------------------------------------------------------------- */

int DumpImage::modify_param(int narg, char **arg)
{
  int n = DumpParticle::modify_param(narg,arg);
  if (n) return n;

  if (strcmp(arg[0],"pcolor") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");
    int err,nlo,nhi;
    if (pcolor == TYPE)
      err = MathExtra::bounds(arg[1],particle->nspecies,nlo,nhi);
    else if (pcolor == PROC)
      err = MathExtra::bounds(arg[1],nprocs,nlo,nhi);
    else error->all(FLERR,"Illegal dump_modify command");
    if (err) error->all(FLERR,"Illegal dump_modify command");

    // ptrs = list of ncount colornames separated by '/'

    int ncount = 1;
    char *nextptr;
    char *ptr = arg[2];
    while ((nextptr = strchr(ptr,'/'))) {
      ptr = nextptr + 1;
      ncount++;
    }
    char **ptrs = new char*[ncount+1];
    ncount = 0;
    ptrs[ncount++] = strtok(arg[2],"/");
    while ((ptrs[ncount++] = strtok(NULL,"/")));
    ncount--;

    // assign each of ncount colors in round-robin fashion to types or procs
    // for PROC case, map 0-Nprocs-1 proc ID to 1 to Nprocs colors

    if (pcolor == TYPE) {
      int m = 0;
      for (int i = nlo; i <= nhi; i++) {
        pcolortype[i] = image->color2rgb(ptrs[m%ncount]);
        if (pcolortype[i] == NULL)
          error->all(FLERR,"Invalid color in dump_modify command");
        m++;
      }
    } else if (pcolor == PROC) {
      if (me+1 >= nlo && me+1 <= nhi) {
        int m = (me+1-nlo) % ncount;
        pcolorproc = image->color2rgb(ptrs[m]);
        if (pcolorproc == NULL)
          error->all(FLERR,"Invalid color in dump_modify command");
      }
    }

    delete [] ptrs;
    return 3;
  }

  if (strcmp(arg[0],"pdiam") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");
    int err,nlo,nhi;
    err = MathExtra::bounds(arg[1],particle->nspecies,nlo,nhi);
    if (err) error->all(FLERR,"Illegal dump_modify command");
    double diam = atof(arg[2]);
    if (diam <= 0.0) error->all(FLERR,"Illegal dump_modify command");
    for (int i = nlo; i <= nhi; i++) pdiamtype[i] = diam;
    return 3;
  }

  if (strcmp(arg[0],"backcolor") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    double *color = image->color2rgb(arg[1]);
    if (color == NULL) error->all(FLERR,"Invalid color in dump_modify command");
    image->background[0] = static_cast<int> (color[0]*255.0);
    image->background[1] = static_cast<int> (color[1]*255.0);
    image->background[2] = static_cast<int> (color[2]*255.0);
    return 2;
  }

  if (strcmp(arg[0],"boxcolor") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    boxcolor = image->color2rgb(arg[1]);
    if (boxcolor == NULL)
      error->all(FLERR,"Invalid color in dump_modify command");
    return 2;
  }

  if (strcmp(arg[0],"cmap") == 0) {
    if (narg < 7) error->all(FLERR,"Illegal dump_modify command");
    int which;
    if (strcmp(arg[1],"particle") == 0) which = PARTICLE;
    else if (strcmp(arg[1],"grid") == 0) which = GRID;
    else if (strcmp(arg[1],"surf") == 0) which = SURF;
    else if (strcmp(arg[1],"gridx") == 0) which = XPLANE;
    else if (strcmp(arg[1],"gridy") == 0) which = YPLANE;
    else if (strcmp(arg[1],"gridz") == 0) which = ZPLANE;
    else error->all(FLERR,"Illegal dump_modify command");
    if (strlen(arg[4]) != 2) error->all(FLERR,"Illegal dump_modify command");
    int factor = 2;
    if (arg[4][0] == 's') factor = 1;
    int nentry = atoi(arg[6]);
    if (nentry < 1) error->all(FLERR,"Illegal dump_modify command");
    int n = 7 + factor*nentry;
    if (narg < n) error->all(FLERR,"Illegal dump_modify command");
    int flag = image->map_reset(which,n-2,&arg[2]);
    if (flag) error->all(FLERR,"Illegal dump_modify command");
    return n;
  }

  if (strcmp(arg[0],"color") == 0) {
    if (narg < 5) error->all(FLERR,"Illegal dump_modify command");
    int flag = image->addcolor(arg[1],atof(arg[2]),atof(arg[3]),atof(arg[4]));
    if (flag) error->all(FLERR,"Illegal dump_modify command");
    return 5;
  }

  if (strcmp(arg[0],"gcolor") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");
    int err,nlo,nhi;
    if (gcolor == PROC) err = MathExtra::bounds(arg[1],nprocs,nlo,nhi);
    else error->all(FLERR,"Illegal dump_modify command");
    if (err) error->all(FLERR,"Illegal dump_modify command");

    // ptrs = list of ncount colornames separated by '/'

    int ncount = 1;
    char *nextptr;
    char *ptr = arg[2];
    while ((nextptr = strchr(ptr,'/'))) {
      ptr = nextptr + 1;
      ncount++;
    }
    char **ptrs = new char*[ncount+1];
    ncount = 0;
    ptrs[ncount++] = strtok(arg[2],"/");
    while ((ptrs[ncount++] = strtok(NULL,"/")));
    ncount--;

    // assign each of ncount colors in round-robin fashion to procs
    // for PROC case, assign Ith color to I-1 value in colorproc
    // this is so can use bounds() above from 1 to Nprocs inclusive

    if (gcolor == PROC) {
      if (me+1 >= nlo && me+1 <= nhi) {
        int m = (me+1-nlo) % ncount;
        gcolorproc = image->color2rgb(ptrs[m]);
        if (gcolorproc == NULL)
          error->all(FLERR,"Invalid color in dump_modify command");
      }
    }

    delete [] ptrs;
    return 3;
  }

  if (strcmp(arg[0],"glinecolor") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    glinecolor = image->color2rgb(arg[1]);
    if (glinecolor == NULL)
      error->all(FLERR,"Invalid color in dump_modify command");
    return 2;
  }

  if (strcmp(arg[0],"scolor") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");
    int err,nlo,nhi;
    if (scolor == ONE) err = 0;
    else if (scolor == PROC) err = MathExtra::bounds(arg[1],nprocs,nlo,nhi);
    else error->all(FLERR,"Illegal dump_modify command");
    if (err) error->all(FLERR,"Illegal dump_modify command");

    // ptrs = list of ncount colornames separated by '/'

    int ncount = 1;
    char *nextptr;
    char *ptr = arg[2];
    while ((nextptr = strchr(ptr,'/'))) {
      ptr = nextptr + 1;
      ncount++;
    }
    char **ptrs = new char*[ncount+1];
    ncount = 0;
    ptrs[ncount++] = strtok(arg[2],"/");
    while ((ptrs[ncount++] = strtok(NULL,"/")));
    ncount--;

    // assign each of ncount colors in round-robin fashion to procs
    // for PROC case, assign Ith color to I-1 value in colorproc
    // this is so can use bounds() above from 1 to Nprocs inclusive

    if (scolor == ONE) {
      if (ncount > 1) error->all(FLERR,"Illegal dump_modify command");
      surfcolorone = image->color2rgb(arg[2]);
      if (surfcolorone == NULL)
        error->all(FLERR,"Invalid color in dump_modify command");
    } else if (scolor == PROC) {
      if (me+1 >= nlo && me+1 <= nhi) {
        int m = (me+1-nlo) % ncount;
        scolorproc = image->color2rgb(ptrs[m]);
        if (scolorproc == NULL)
          error->all(FLERR,"Invalid color in dump_modify command");
      }
    }

    delete [] ptrs;
    return 3;
  }

  if (strcmp(arg[0],"slinecolor") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    slinecolor = image->color2rgb(arg[1]);
    if (slinecolor == NULL)
      error->all(FLERR,"Invalid color in dump_modify command");
    return 2;
  }

  return 0;
}
