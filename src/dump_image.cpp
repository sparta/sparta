/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "ctype.h"
#include "stdlib.h"
#include "string.h"
#include "dump_image.h"
#include "image.h"
#include "domain.h"
#include "particle.h"
#include "input.h"
#include "variable.h"
#include "math_extra.h"
#include "math_const.h"
#include "error.h"
#include "memory.h"

#ifdef DSMC_JPEG
#include "jpeglib.h"
#endif

using namespace DSMC_NS;
using namespace MathConst;

enum{PPM,JPG};
enum{NUMERIC,ATOM,TYPE,ELEMENT,ATTRIBUTE};
enum{STATIC,DYNAMIC};
enum{NO,YES};

/* ---------------------------------------------------------------------- */

DumpImage::DumpImage(DSMC *dsmc, int narg, char **arg) : 
  DumpMolecule(dsmc, narg, arg)
{
  if (binary || multiproc) error->all(FLERR,"Invalid dump image filename");

  // set filetype based on filename suffix

  int n = strlen(filename);
  if (strlen(filename) > 4 && strcmp(&filename[n-4],".jpg") == 0)
    filetype = JPG;
  else if (strlen(filename) > 5 && strcmp(&filename[n-5],".jpeg") == 0)
    filetype = JPG;
  else filetype = PPM;

#ifndef DSMC_JPEG
  if (filetype == JPG) error->all(FLERR,"Cannot dump JPG file");
#endif

  // atom color,diameter settings

  if (nfield != 2) error->all(FLERR,"Illegal dump image command");

  acolor = ATTRIBUTE;
  if (strcmp(arg[5],"type") == 0) acolor = TYPE;
  else if (strcmp(arg[5],"element") == 0) acolor = ELEMENT;

  adiam = ATTRIBUTE;
  if (strcmp(arg[6],"type") == 0) adiam = TYPE;
  else if (strcmp(arg[6],"element") == 0) adiam = ELEMENT;

  // create Image class

  image = new Image(dsmc);

  // set defaults for optional args

  atomflag = YES;
  thetastr = phistr = NULL;
  cflag = STATIC;
  cx = cy = cz = 0.5;
  cxstr = cystr = czstr = NULL;

  if (domain->dimension == 3) {
    image->up[0] = 0.0; image->up[1] = 0.0; image->up[2] = 1.0;
  } else {
    image->up[0] = 0.0; image->up[1] = 1.0; image->up[2] = 0.0;
  }

  upxstr = upystr = upzstr = NULL;
  zoomstr = NULL;
  perspstr = NULL;
  boxflag = YES;
  boxdiam = 0.02;
  axesflag = NO;

  // parse optional args

  int iarg = ioptional;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"adiam") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      adiam = NUMERIC;
      adiamvalue = atof(arg[iarg+1]);
      if (adiamvalue <= 0.0) error->all(FLERR,"Illegal dump image command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"atom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"yes") == 0) atomflag = YES;
      else if (strcmp(arg[iarg+1],"no") == 0) atomflag = NO;
      else error->all(FLERR,"Illegal dump image command");
      iarg += 2;

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
      if (strcmp(arg[iarg+1],"yes") == 0) boxflag = YES;
      else if (strcmp(arg[iarg+1],"no") == 0) boxflag = NO;
      else error->all(FLERR,"Illegal dump image command");
      boxdiam = atof(arg[iarg+2]);
      if (boxdiam < 0.0) error->all(FLERR,"Illegal dump image command");
      iarg += 3;

    } else if (strcmp(arg[iarg],"axes") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"yes") == 0) axesflag = YES;
      else if (strcmp(arg[iarg+1],"no") == 0) axesflag = NO;
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
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"yes") == 0) image->ssao = YES;
      else if (strcmp(arg[iarg+1],"no") == 0) image->ssao = NO;
      else error->all(FLERR,"Illegal dump image command");
      double ssaoint = atof(arg[iarg+2]);
      if (ssaoint < 0.0 || ssaoint > 1.0)
	error->all(FLERR,"Illegal dump image command");
      image->ssaoint = ssaoint;
      iarg += 3;

    } else error->all(FLERR,"Illegal dump image command");
  }

  // allocate image buffer now that image size is known

  image->buffers();

  // additional defaults for dump_modify options

  diamtype = new double[ntypes+1];
  diamelement = new double[ntypes+1];
  colortype = new double*[ntypes+1];
  colorelement = new double*[ntypes+1];

  for (int i = 1; i <= ntypes; i++) {
    diamtype[i] = 1.0;
    if (i % 6 == 1) colortype[i] = image->color2rgb("red");
    else if (i % 6 == 2) colortype[i] = image->color2rgb("green");
    else if (i % 6 == 3) colortype[i] = image->color2rgb("blue");
    else if (i % 6 == 4) colortype[i] = image->color2rgb("yellow");
    else if (i % 6 == 5) colortype[i] = image->color2rgb("aqua");
    else if (i % 6 == 0) colortype[i] = image->color2rgb("cyan");
  }

  // viewflag = DYNAMIC if any view parameter is dynamic

  viewflag = STATIC;
  if (thetastr || phistr || cflag == DYNAMIC || 
      upxstr || upystr || upzstr || zoomstr || perspstr) viewflag = DYNAMIC;

  if (cflag == STATIC) box_center();
  if (viewflag == STATIC) view_params();

  // local data

  maxbufcopy = 0;
  bufcopy = NULL;
}

/* ---------------------------------------------------------------------- */

DumpImage::~DumpImage()
{
  delete image;

  delete [] diamtype;
  delete [] diamelement;
  delete [] colortype;
  delete [] colorelement;
}

/* ---------------------------------------------------------------------- */

void DumpImage::init_style()
{
  if (multifile == 0) 
    error->all(FLERR,"Dump image requires one snapshot per file");

  DumpMolecule::init_style();

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

  // set up type -> element mapping

  if (atomflag && acolor == ELEMENT) {
    for (int i = 1; i <= ntypes; i++) {
      colorelement[i] = image->element2color(typenames[i]);
      if (colorelement[i] == NULL)
	error->all(FLERR,"Invalid dump image element name");
    }
  }

  if (atomflag && adiam == ELEMENT) {
    for (int i = 1; i <= ntypes; i++) {
      diamelement[i] = image->element2diam(typenames[i]);
      if (diamelement[i] == 0.0)
	error->all(FLERR,"Invalid dump image element name");
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

  // nme = # of atoms this proc will contribute to dump
  // pack buf with x,y,z,color,diameter
  // set minmax color range if using color map
  // create my portion of image for my particles
  
  nme = count();

  if (nme > maxbuf) {
    maxbuf = nme;
    memory->destroy(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
  }

  pack();
  if (acolor == ATTRIBUTE) image->color_minmax(nchoose,buf,size_one);

  // create image on each proc, then merge them

  image->clear();
  create_image();
  image->merge();

  // write image file

  if (me == 0) {
    if (filetype == JPG) image->write_JPG(fp);
    else image->write_PPM(fp);
    fclose(fp);
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
   create image for atoms on this proc
   every pixel has depth 
------------------------------------------------------------------------- */

void DumpImage::create_image()
{
  int i,j,m,itype,atom1,atom2;
  double diameter,delx,dely,delz;
  double *color,*color1,*color2;
  double xmid[3];

  // render my atoms

  if (atomflag) {
    Particle::OnePart *particles = particle->particles;

    m = 0;
    for (i = 0; i < nchoose; i++) {
      j = clist[i];
      
      if (acolor == TYPE) {
	itype = static_cast<int> (buf[m]);
	color = colortype[itype];
      } else if (acolor == ELEMENT) {
	itype = static_cast<int> (buf[m]);
	color = colorelement[itype];
      } else if (acolor == ATTRIBUTE) {
	color = image->value2color(buf[m]);
      }

      if (adiam == NUMERIC) {
	diameter = adiamvalue;
      } else if (adiam == TYPE) {
	itype = static_cast<int> (buf[m+1]);
	diameter = diamtype[itype];
      } else if (adiam == ELEMENT) {
	itype = static_cast<int> (buf[m+1]);
	diameter = diamelement[itype];
      } else if (adiam == ATTRIBUTE) {
	diameter = buf[m+1];
      }

      image->draw_sphere(particles[j].x,color,diameter);
      m += size_one;
    }
  }

  // render outline of simulation box, orthogonal or triclinic

  if (boxflag) {
    double diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
    if (domain->dimension == 3) diameter = MIN(diameter,boxzhi-boxzlo);
    diameter *= boxdiam;

    double (*boxcorners)[3];
    double box[8][3];
    box[0][0] = boxxlo; box[0][1] = boxylo; box[0][2] = boxzlo;
    box[1][0] = boxxhi; box[1][1] = boxylo; box[1][2] = boxzlo;
    box[2][0] = boxxlo; box[2][1] = boxyhi; box[2][2] = boxzlo;
    box[3][0] = boxxhi; box[3][1] = boxyhi; box[3][2] = boxzlo;
    box[4][0] = boxxlo; box[4][1] = boxylo; box[4][2] = boxzhi;
    box[5][0] = boxxhi; box[5][1] = boxylo; box[5][2] = boxzhi;
    box[6][0] = boxxlo; box[6][1] = boxyhi; box[6][2] = boxzhi;
    box[7][0] = boxxhi; box[7][1] = boxyhi; box[7][2] = boxzhi;
    boxcorners = box;

    image->draw_box(box,diameter);
  }

  // render XYZ axes in red/green/blue
  // offset by 10% of box size and scale by axeslen

  if (axesflag) {
    double diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
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
}

/* ---------------------------------------------------------------------- */

int DumpImage::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;

  if (comm_forward == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = choose[j];
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = choose[j];
      buf[m++] = bufcopy[j][0];
      buf[m++] = bufcopy[j][1];
    }
  }

  return comm_forward;
}

/* ---------------------------------------------------------------------- */

void DumpImage::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  if (comm_forward == 1)
    for (i = first; i < last; i++) choose[i] = static_cast<int> (buf[m++]);
  else {
    for (i = first; i < last; i++) {
      choose[i] = static_cast<int> (buf[m++]);
      bufcopy[i][0] = buf[m++];
      bufcopy[i][1] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int DumpImage::modify_param(int narg, char **arg)
{
  int n = DumpMolecule::modify_param(narg,arg);
  if (n) return n;

  if (strcmp(arg[0],"acolor") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");
    int nlo,nhi;
    MathExtra::bounds(arg[1],particle->nspecies,nlo,nhi);

    // ptrs = list of ncount colornames separated by '/'

    int ncount = 1;
    char *nextptr;
    char *ptr = arg[2];
    while (nextptr = strchr(ptr,'/')) {
      ptr = nextptr + 1;
      ncount++;
    }
    char **ptrs = new char*[ncount+1];
    ncount = 0;
    ptrs[ncount++] = strtok(arg[2],"/");
    while (ptrs[ncount++] = strtok(NULL,"/"));
    ncount--;

    // assign each of ncount colors in round-robin fashion to types

    int m = 0;
    for (int i = nlo; i <= nhi; i++) {
      colortype[i] = image->color2rgb(ptrs[m%ncount]);
      if (colortype[i] == NULL)
	error->all(FLERR,"Invalid color in dump_modify command");
      m++;
    }

    delete [] ptrs;
    return 3;
  }

  if (strcmp(arg[0],"adiam") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");
    int nlo,nhi;
    MathExtra::bounds(arg[1],particle->nspecies,nlo,nhi);
    double diam = atof(arg[2]);
    if (diam <= 0.0) error->all(FLERR,"Illegal dump_modify command");
    for (int i = nlo; i <= nhi; i++) diamtype[i] = diam;
    return 3;
  }

  if (strcmp(arg[0],"amap") == 0) {
    if (narg < 6) error->all(FLERR,"Illegal dump_modify command");
    if (strlen(arg[3]) != 2) error->all(FLERR,"Illegal dump_modify command");
    int factor = 2;
    if (arg[3][0] == 's') factor = 1;
    int nentry = atoi(arg[5]);
    if (nentry < 1) error->all(FLERR,"Illegal dump_modify command");
    int n = 6 + factor*nentry;
    if (narg < n) error->all(FLERR,"Illegal dump_modify command");
    int flag = image->colormap(n-1,&arg[1]);
    if (flag) error->all(FLERR,"Illegal dump_modify command");
    return n;
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
    image->boxcolor = image->color2rgb(arg[1]);
    if (image->boxcolor == NULL) 
      error->all(FLERR,"Invalid color in dump_modify command");
    return 2;
  }

  if (strcmp(arg[0],"color") == 0) {
    if (narg < 5) error->all(FLERR,"Illegal dump_modify command");
    int flag = image->addcolor(arg[1],atof(arg[2]),atof(arg[3]),atof(arg[4]));
    if (flag) error->all(FLERR,"Illegal dump_modify command");
    return 5;
  }

  return 0;
}
