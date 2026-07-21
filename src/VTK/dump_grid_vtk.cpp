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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "dump_grid_vtk.h"
#include "grid.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

#include <vtkVersion.h>
#include <vtkCellType.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkLongLongArray.h>
#include <vtkStringArray.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

using namespace SPARTA_NS;

enum{INT,DOUBLE,BIGINT,UINT,BIGUINT,STRING};    // same as Dump
enum{VTK,VTU,PVTU};                             // VTK file formats (grid)

/* ---------------------------------------------------------------------- */

DumpGridVTK::DumpGridVTK(SPARTA *sparta, int narg, char **arg) :
  DumpGrid(sparta, narg, arg)
{
  buffer_allow = 0;
  buffer_flag = 0;

  n_calls_ = 0;
  filecurrent = NULL;
  parallelfilecurrent = NULL;

  // VTK cell geometry: voxel (3d) or pixel (2d)

  if (dimension == 3) { celltype = VTK_VOXEL; ncorner = 8; ngeom = 6; }
  else                { celltype = VTK_PIXEL; ncorner = 4; ngeom = 4; }

  // determine file format from extension: .vtu / .vtk (unstructured only)

  vtk_file_format = VTK;
  char *suffix = filename + strlen(filename) - strlen(".vtu");
  if (suffix > filename && strcmp(suffix,".vtu") == 0)
    vtk_file_format = multiproc ? PVTU : VTU;
  else {
    suffix = filename + strlen(filename) - strlen(".vtp");
    if (suffix > filename && strcmp(suffix,".vtp") == 0)
      error->all(FLERR,"Dump grid/vtk does not support .vtp; grid cells are "
                 "volumetric, use .vtu");
  }

  if (vtk_file_format == VTK && multiproc)
    error->all(FLERR,"Dump grid/vtk legacy .vtk format does not support "
               "'%' in filename; use .vtu");

  if (!multifile)
    error->all(FLERR,"Dump grid/vtk requires '*' in filename for "
               "per-timestep output");

  // build cell-data field list from the user attributes, then append
  // the box-bound geometry columns to the packed buffer

  ndata = size_one;
  setup_vtk_fields();
  augment_geometry_fields();

  if (filewriter) reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

DumpGridVTK::~DumpGridVTK()
{
  delete [] filecurrent;
  delete [] parallelfilecurrent;

  // free the vformat entries for the appended geometry columns
  // (the DumpGrid destructor only frees [0,nfield) = the user columns)

  for (int i = ndata; i < size_one; i++) {
    delete [] vformat[i];
    vformat[i] = NULL;
  }
}

/* ----------------------------------------------------------------------
   build one cell-data array descriptor per user data column
------------------------------------------------------------------------- */

void DumpGridVTK::setup_vtk_fields()
{
  std::vector<std::string> names;
  char *copy = new char[strlen(columns)+1];
  strcpy(copy,columns);
  char *tok = strtok(copy," ");
  while (tok) { names.push_back(tok); tok = strtok(NULL," "); }
  delete [] copy;

  if ((int) names.size() != ndata)
    error->all(FLERR,"Dump grid/vtk internal field count mismatch");

  fields.clear();
  for (int i = 0; i < ndata; i++) {
    VTKField f;
    f.col = i;
    f.type = vtype[i];
    f.name = names[i];
    fields.push_back(f);
  }
}

/* ----------------------------------------------------------------------
   append ngeom geometry columns (the cell box bounds) to the end of the
   packed buffer so the geometry travels through the MPI gather and the
   filewriter can rebuild cell corners in buf2arrays().

   MAINTENANCE NOTE: this reaches into and reallocates arrays that are
   allocated and freed by the base classes (Dump / DumpGrid):
     size_one          - grown to ndata+ngeom (drives buf width + pack loop)
     vtype             - grown; freed by ~DumpGrid (delete [] vtype)
     format_default    - rebuilt; freed by ~Dump
     format_column_user- regrown to size_one; freed by ~Dump over [0,size_one)
     vformat           - regrown to size_one; ~DumpGrid frees only [0,nfield),
                         so ~DumpGridVTK frees the [ndata,size_one) tail
   The geometry pack functions are NOT stored in pack_choice (that stays
   size ndata); pack() invokes them directly.  If a base dump changes how it
   allocates or frees any of the above, this and ~DumpGridVTK must be updated
   in lock-step.
------------------------------------------------------------------------- */

void DumpGridVTK::augment_geometry_fields()
{
  int newsize = ndata + ngeom;

  // grow vtype for the appended geometry columns (all DOUBLE); the geometry
  // packers themselves are invoked directly in pack(), not via pack_choice

  int *newvtype = new int[newsize];
  for (int i = 0; i < ndata; i++) newvtype[i] = vtype[i];
  for (int i = ndata; i < newsize; i++) newvtype[i] = DOUBLE;
  delete [] vtype; vtype = newvtype;
  size_one = newsize;

  // resize the format bookkeeping arrays that Dump::init() indexes by size_one

  delete [] vformat;
  vformat = new char*[newsize];
  for (int i = 0; i < newsize; i++) vformat[i] = NULL;

  delete [] format_column_user;
  format_column_user = new char*[newsize];
  for (int i = 0; i < newsize; i++) format_column_user[i] = NULL;

  delete [] format_default;
  format_default = new char[4*newsize+1];
  format_default[0] = '\0';
  for (int i = 0; i < newsize; i++) {
    if (vtype[i] == DOUBLE) strcat(format_default,"%g ");
    else if (vtype[i] == INT) strcat(format_default,"%d ");
    else if (vtype[i] == BIGINT) strcat(format_default,BIGINT_FORMAT " ");
    else if (vtype[i] == UINT) strcat(format_default,"%u ");
    else if (vtype[i] == BIGUINT) strcat(format_default,BIGUINT_FORMAT " ");
    else strcat(format_default,"%g ");
  }
}

/* ---------------------------------------------------------------------- */

void DumpGridVTK::init_style()
{
  // reuse DumpGrid to resolve compute/fix/variable/custom and build cpart;
  // multifile is enforced in the constructor, so no dump file is opened here

  DumpGrid::init_style();
}

/* ---------------------------------------------------------------------- */

// the VTK writers manage their own files; see dump_particle_vtk.cpp
void DumpGridVTK::openfile() {}

void DumpGridVTK::write_header(bigint) { n_calls_ = 0; }

/* ----------------------------------------------------------------------
   pack the user data columns via the inherited pack functions, then pack
   the appended box-bound geometry columns by calling the inherited DumpGrid
   pack methods directly (order must match buf2arrays)
------------------------------------------------------------------------- */

void DumpGridVTK::pack()
{
  for (int n = 0; n < ndata; n++) (this->*pack_choice[n])(n);

  int g = ndata;
  pack_xlo(g++);
  pack_ylo(g++);
  if (dimension == 3) pack_zlo(g++);
  pack_xhi(g++);
  pack_yhi(g++);
  if (dimension == 3) pack_zhi(g++);
}

/* ---------------------------------------------------------------------- */

void DumpGridVTK::write_data(int n, double *mybuf)
{
  if (vtk_file_format == VTU || vtk_file_format == PVTU) write_vtu(n,mybuf);
  else write_vtk(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpGridVTK::buf2arrays(int n, double *mybuf)
{
  for (int i = 0; i < n; i++) {
    int base = i*size_one;
    double *g = &mybuf[base+ndata];

    vtkIdType ids[8];
    if (dimension == 3) {
      double xlo=g[0],ylo=g[1],zlo=g[2],xhi=g[3],yhi=g[4],zhi=g[5];
      ids[0] = points->InsertNextPoint(xlo,ylo,zlo);
      ids[1] = points->InsertNextPoint(xhi,ylo,zlo);
      ids[2] = points->InsertNextPoint(xlo,yhi,zlo);
      ids[3] = points->InsertNextPoint(xhi,yhi,zlo);
      ids[4] = points->InsertNextPoint(xlo,ylo,zhi);
      ids[5] = points->InsertNextPoint(xhi,ylo,zhi);
      ids[6] = points->InsertNextPoint(xlo,yhi,zhi);
      ids[7] = points->InsertNextPoint(xhi,yhi,zhi);
    } else {
      double xlo=g[0],ylo=g[1],xhi=g[2],yhi=g[3];
      ids[0] = points->InsertNextPoint(xlo,ylo,0.0);
      ids[1] = points->InsertNextPoint(xhi,ylo,0.0);
      ids[2] = points->InsertNextPoint(xlo,yhi,0.0);
      ids[3] = points->InsertNextPoint(xhi,yhi,0.0);
    }
    ugrid->InsertNextCell(celltype,ncorner,ids);

    for (int f = 0; f < (int) fields.size(); f++) {
      int c = base + fields[f].col;
      vtkAbstractArray *paa = myarrays[f];
      if (fields[f].type == DOUBLE)
        ((vtkDoubleArray *) paa)->InsertNextValue(mybuf[c]);
      else if (fields[f].type == STRING) {
        // idstr: buf holds the numeric cell ID (ubuf-encoded); convert to the
        // hierarchical string form, mirroring DumpGrid::write_text
        char str[32];
        grid->id_num2str((cellint) ubuf(mybuf[c]).i,str);
        ((vtkStringArray *) paa)->InsertNextValue(str);
      } else
        ((vtkLongLongArray *) paa)->
          InsertNextValue((long long) ubuf(mybuf[c]).i);
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpGridVTK::reset_vtk_data_containers()
{
  points = vtkSmartPointer<vtkPoints>::New();
  ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  ugrid->Allocate();

  myarrays.clear();
  for (int f = 0; f < (int) fields.size(); f++) {
    vtkSmartPointer<vtkAbstractArray> arr;
    if (fields[f].type == DOUBLE)
      arr = vtkSmartPointer<vtkDoubleArray>::New();
    else if (fields[f].type == STRING)
      arr = vtkSmartPointer<vtkStringArray>::New();
    else
      arr = vtkSmartPointer<vtkLongLongArray>::New();
    arr->SetName(fields[f].name.c_str());
    myarrays.push_back(arr);
  }
}

/* ---------------------------------------------------------------------- */

void DumpGridVTK::write_vtk(int n, double *mybuf)
{
  ++n_calls_;
  buf2arrays(n,mybuf);
  if (n_calls_ < nclusterprocs) return;

  setFileCurrent();

  ugrid->SetPoints(points);
  for (int f = 0; f < (int) myarrays.size(); f++)
    ugrid->GetCellData()->AddArray(myarrays[f]);

  vtkSmartPointer<vtkUnstructuredGridWriter> writer =
    vtkSmartPointer<vtkUnstructuredGridWriter>::New();
  writer->SetHeader("Generated by SPARTA");
  if (binary) writer->SetFileTypeToBinary();
  else        writer->SetFileTypeToASCII();
#if VTK_MAJOR_VERSION < 6
  writer->SetInput(ugrid);
#else
  writer->SetInputData(ugrid);
#endif
  writer->SetFileName(filecurrent);
  writer->Write();

  reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

void DumpGridVTK::write_vtu(int n, double *mybuf)
{
  ++n_calls_;
  buf2arrays(n,mybuf);
  if (n_calls_ < nclusterprocs) return;

  setFileCurrent();

  ugrid->SetPoints(points);
  for (int f = 0; f < (int) myarrays.size(); f++)
    ugrid->GetCellData()->AddArray(myarrays[f]);

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  if (binary) writer->SetDataModeToBinary();
  else        writer->SetDataModeToAscii();
#if VTK_MAJOR_VERSION < 6
  writer->SetInput(ugrid);
#else
  writer->SetInputData(ugrid);
#endif
  writer->SetFileName(filecurrent);
  writer->Write();

  if (me == 0 && multiproc) write_pvtk();

  reset_vtk_data_containers();
}

/* ----------------------------------------------------------------------
   filename helpers (see dump_particle_vtk.cpp for the same pattern)
------------------------------------------------------------------------- */

static char *subst_percent(const char *name, int id)
{
  const char *ptr = strchr(name,'%');
  if (!ptr) { char *s = new char[strlen(name)+1]; strcpy(s,name); return s; }
  char *s = new char[strlen(name)+16];
  int pre = (int)(ptr-name);
  strncpy(s,name,pre); s[pre] = '\0';
  sprintf(s+pre,"%d%s",id,ptr+1);
  return s;
}

static char *subst_star(const char *name, bigint ntimestep, int padflag)
{
  const char *ptr = strchr(name,'*');
  if (!ptr) { char *s = new char[strlen(name)+1]; strcpy(s,name); return s; }
  char *s = new char[strlen(name)+16];
  int pre = (int)(ptr-name);
  strncpy(s,name,pre); s[pre] = '\0';
  if (padflag == 0)
    sprintf(s+pre,BIGINT_FORMAT "%s",ntimestep,ptr+1);
  else {
    char bif[8],pad[16];
    strcpy(bif,BIGINT_FORMAT);
    sprintf(pad,"%%0%d%s%%s",padflag,&bif[1]);
    sprintf(s+pre,pad,ntimestep,ptr+1);
  }
  return s;
}

/* ---------------------------------------------------------------------- */

void DumpGridVTK::setFileCurrent()
{
  delete [] filecurrent;
  filecurrent = NULL;

  int id = me;
  if (multiproc > 1)
    id = (me + nclusterprocs == nprocs) ? multiproc-1 : me/nclusterprocs;

  char *tmp = subst_percent(filename,id);
  filecurrent = subst_star(tmp,update->ntimestep,padflag);
  delete [] tmp;

  if (multiproc && me == 0) {
    delete [] parallelfilecurrent;
    parallelfilecurrent = NULL;

    char *ptr = strchr(filename,'%');
    char *base = new char[strlen(filename)+16];
    int pre = (int)(ptr-filename);
    strncpy(base,filename,pre); base[pre] = '\0';
    sprintf(base+pre,"%s",ptr+1);
    char *ext = strrchr(base,'.');
    ext++;
    ext[0] = 'p'; ext[1] = 'v'; ext[2] = 't'; ext[3] = 'u'; ext[4] = '\0';

    parallelfilecurrent = subst_star(base,update->ntimestep,padflag);
    delete [] base;
  }
}

/* ---------------------------------------------------------------------- */

std::string DumpGridVTK::pvtk_piece_filename(int id)
{
  char *tmp = subst_percent(filename,id);
  char *full = subst_star(tmp,update->ntimestep,padflag);
  delete [] tmp;
  std::string fname = full;
  delete [] full;
  size_t pos = fname.find_last_of("/\\");
  if (pos != std::string::npos) fname = fname.substr(pos+1);
  return fname;
}

/* ---------------------------------------------------------------------- */

void DumpGridVTK::write_pvtk()
{
  FILE *fp = fopen(parallelfilecurrent,"w");
  if (!fp) error->one(FLERR,"Cannot open dump grid/vtk parallel file");

  int one = 1;
  const char *byte_order =
    (*((char *) &one)) ? "LittleEndian" : "BigEndian";

  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" "
          "byte_order=\"%s\">\n",byte_order);
  fprintf(fp,"  <PUnstructuredGrid GhostLevel=\"0\">\n");

  fprintf(fp,"    <PCellData>\n");
  for (int f = 0; f < (int) fields.size(); f++) {
    const char *type = "Int64";
    if (fields[f].type == DOUBLE) type = "Float64";
    else if (fields[f].type == STRING) type = "String";
    fprintf(fp,"      <PDataArray type=\"%s\" Name=\"%s\"/>\n",
            type,fields[f].name.c_str());
  }
  fprintf(fp,"    </PCellData>\n");

  fprintf(fp,"    <PPoints>\n");
  fprintf(fp,"      <PDataArray type=\"Float32\" Name=\"Points\" "
          "NumberOfComponents=\"3\"/>\n");
  fprintf(fp,"    </PPoints>\n");

  int npieces = (multiproc > 1) ? multiproc : nprocs;
  for (int i = 0; i < npieces; i++)
    fprintf(fp,"    <Piece Source=\"%s\"/>\n",pvtk_piece_filename(i).c_str());

  fprintf(fp,"  </PUnstructuredGrid>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}

/* ---------------------------------------------------------------------- */

int DumpGridVTK::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"binary") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"yes") == 0) binary = 1;
    else if (strcmp(arg[1],"no") == 0) binary = 0;
    else error->all(FLERR,"Illegal dump_modify command");
    return 2;
  }
  return 0;
}
