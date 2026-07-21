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
#include "dump_surf_vtk.h"
#include "surf.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

#include <vtkVersion.h>
#include <vtkCellType.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkLongLongArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

using namespace SPARTA_NS;

enum{INT,DOUBLE,BIGINT,UINT,BIGUINT,STRING};    // same as Dump
enum{VTK,VTP,VTU,PVTP,PVTU};                    // VTK file formats

/* ---------------------------------------------------------------------- */

DumpSurfVTK::DumpSurfVTK(SPARTA *sparta, int narg, char **arg) :
  DumpSurf(sparta, narg, arg)
{
  buffer_allow = 0;
  buffer_flag = 0;

  n_calls_ = 0;
  filecurrent = NULL;
  parallelfilecurrent = NULL;

  // surf element geometry: line (2d) or triangle (3d)

  if (dimension == 3) { ncorner = 3; ngeom = 9; }
  else                { ncorner = 2; ngeom = 4; }

  // determine file format from extension: .vtp, .vtu, or legacy .vtk

  vtk_file_format = VTK;
  char *suffix = filename + strlen(filename) - strlen(".vtp");
  if (suffix > filename && strcmp(suffix,".vtp") == 0) {
    vtk_file_format = multiproc ? PVTP : VTP;
  } else {
    suffix = filename + strlen(filename) - strlen(".vtu");
    if (suffix > filename && strcmp(suffix,".vtu") == 0)
      vtk_file_format = multiproc ? PVTU : VTU;
  }

  if (vtk_file_format == VTK && multiproc)
    error->all(FLERR,"Dump surf/vtk legacy .vtk format does not support "
               "'%' in filename; use .vtu or .vtp");

  if (!multifile)
    error->all(FLERR,"Dump surf/vtk requires '*' in filename for "
               "per-timestep output");

  ndata = size_one;
  setup_vtk_fields();
  augment_geometry_fields();

  if (filewriter) reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

DumpSurfVTK::~DumpSurfVTK()
{
  delete [] filecurrent;
  delete [] parallelfilecurrent;

  // free vformat entries for the appended geometry columns
  // (the DumpSurf destructor only frees [0,nfield) = the user columns)

  for (int i = ndata; i < size_one; i++) {
    delete [] vformat[i];
    vformat[i] = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurfVTK::setup_vtk_fields()
{
  std::vector<std::string> names;
  char *copy = new char[strlen(columns)+1];
  strcpy(copy,columns);
  char *tok = strtok(copy," ");
  while (tok) { names.push_back(tok); tok = strtok(NULL," "); }
  delete [] copy;

  if ((int) names.size() != ndata)
    error->all(FLERR,"Dump surf/vtk internal field count mismatch");

  fields.clear();
  for (int i = 0; i < ndata; i++) {
    if (vtype[i] == STRING)
      error->all(FLERR,"Dump surf/vtk does not support string attributes");
    VTKField f;
    f.col = i;
    f.type = vtype[i];
    f.name = names[i];
    fields.push_back(f);
  }
}

/* ----------------------------------------------------------------------
   append ngeom vertex columns to the end of the packed buffer so the
   element geometry travels through the MPI gather and the filewriter can
   rebuild line/triangle vertices in buf2arrays().

   MAINTENANCE NOTE: like DumpGridVTK::augment_geometry_fields(), this grows
   base-owned arrays (size_one, vtype, format_default, format_column_user,
   vformat) and relies on the base destructors' free ranges; ~DumpSurfVTK
   frees the vformat [ndata,size_one) tail.  The geometry pack functions are
   invoked directly in pack(), not stored in pack_choice.  Keep this in
   lock-step with the base dumps' allocation/free behavior.
------------------------------------------------------------------------- */

void DumpSurfVTK::augment_geometry_fields()
{
  int newsize = ndata + ngeom;

  // grow vtype for the appended vertex columns (all DOUBLE); the geometry
  // packers themselves are invoked directly in pack(), not via pack_choice

  int *newvtype = new int[newsize];
  for (int i = 0; i < ndata; i++) newvtype[i] = vtype[i];
  for (int i = ndata; i < newsize; i++) newvtype[i] = DOUBLE;
  delete [] vtype; vtype = newvtype;
  size_one = newsize;

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

void DumpSurfVTK::init_style()
{
  // reuse DumpSurf to resolve compute/fix/variable/custom and build the
  // owned-element lists; multifile is enforced, so no dump file is opened

  DumpSurf::init_style();
}

/* ---------------------------------------------------------------------- */

// the VTK writers manage their own files; see dump_particle_vtk.cpp
void DumpSurfVTK::openfile() {}

void DumpSurfVTK::write_header(bigint) { n_calls_ = 0; }

/* ----------------------------------------------------------------------
   pack the user data columns via the inherited pack functions, then pack
   the appended vertex geometry columns by calling the inherited DumpSurf
   pack methods directly (order must match buf2arrays)
------------------------------------------------------------------------- */

void DumpSurfVTK::pack()
{
  for (int n = 0; n < ndata; n++) (this->*pack_choice[n])(n);

  int g = ndata;
  if (dimension == 3) {
    pack_v1x(g++); pack_v1y(g++); pack_v1z(g++);
    pack_v2x(g++); pack_v2y(g++); pack_v2z(g++);
    pack_v3x(g++); pack_v3y(g++); pack_v3z(g++);
  } else {
    pack_v1x(g++); pack_v1y(g++);
    pack_v2x(g++); pack_v2y(g++);
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurfVTK::write_data(int n, double *mybuf)
{
  if (vtk_file_format == VTP || vtk_file_format == PVTP) write_vtp(n,mybuf);
  else if (vtk_file_format == VTU || vtk_file_format == PVTU) write_vtu(n,mybuf);
  else write_vtk(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpSurfVTK::buf2arrays(int n, double *mybuf)
{
  for (int i = 0; i < n; i++) {
    int base = i*size_one;
    double *g = &mybuf[base+ndata];

    vtkIdType ids[3];
    if (dimension == 3) {
      ids[0] = points->InsertNextPoint(g[0],g[1],g[2]);
      ids[1] = points->InsertNextPoint(g[3],g[4],g[5]);
      ids[2] = points->InsertNextPoint(g[6],g[7],g[8]);
    } else {
      ids[0] = points->InsertNextPoint(g[0],g[1],0.0);
      ids[1] = points->InsertNextPoint(g[2],g[3],0.0);
    }
    cellArray->InsertNextCell(ncorner,ids);

    for (int f = 0; f < (int) fields.size(); f++) {
      int c = base + fields[f].col;
      vtkAbstractArray *paa = myarrays[f];
      if (fields[f].type == DOUBLE)
        ((vtkDoubleArray *) paa)->InsertNextValue(mybuf[c]);
      else
        ((vtkLongLongArray *) paa)->
          InsertNextValue((long long) ubuf(mybuf[c]).i);
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurfVTK::reset_vtk_data_containers()
{
  points = vtkSmartPointer<vtkPoints>::New();
  cellArray = vtkSmartPointer<vtkCellArray>::New();

  myarrays.clear();
  for (int f = 0; f < (int) fields.size(); f++) {
    vtkSmartPointer<vtkAbstractArray> arr;
    if (fields[f].type == DOUBLE)
      arr = vtkSmartPointer<vtkDoubleArray>::New();
    else
      arr = vtkSmartPointer<vtkLongLongArray>::New();
    arr->SetName(fields[f].name.c_str());
    myarrays.push_back(arr);
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurfVTK::write_vtk(int n, double *mybuf)
{
  ++n_calls_;
  buf2arrays(n,mybuf);
  if (n_calls_ < nclusterprocs) return;

  setFileCurrent();

  vtkSmartPointer<vtkUnstructuredGrid> ugrid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();
  ugrid->SetPoints(points);
  ugrid->SetCells(dimension == 3 ? VTK_TRIANGLE : VTK_LINE,cellArray);
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

void DumpSurfVTK::write_vtp(int n, double *mybuf)
{
  ++n_calls_;
  buf2arrays(n,mybuf);
  if (n_calls_ < nclusterprocs) return;

  setFileCurrent();

  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(points);
  if (dimension == 3) polyData->SetPolys(cellArray);
  else                polyData->SetLines(cellArray);
  for (int f = 0; f < (int) myarrays.size(); f++)
    polyData->GetCellData()->AddArray(myarrays[f]);

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  if (binary) writer->SetDataModeToBinary();
  else        writer->SetDataModeToAscii();
#if VTK_MAJOR_VERSION < 6
  writer->SetInput(polyData);
#else
  writer->SetInputData(polyData);
#endif
  writer->SetFileName(filecurrent);
  writer->Write();

  if (me == 0 && multiproc) write_pvtk(PVTP);

  reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

void DumpSurfVTK::write_vtu(int n, double *mybuf)
{
  ++n_calls_;
  buf2arrays(n,mybuf);
  if (n_calls_ < nclusterprocs) return;

  setFileCurrent();

  vtkSmartPointer<vtkUnstructuredGrid> ugrid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();
  ugrid->SetPoints(points);
  ugrid->SetCells(dimension == 3 ? VTK_TRIANGLE : VTK_LINE,cellArray);
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

  if (me == 0 && multiproc) write_pvtk(PVTU);

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

void DumpSurfVTK::setFileCurrent()
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
    ext[0] = 'p'; ext[1] = 'v'; ext[2] = 't';
    ext[3] = (vtk_file_format == PVTP) ? 'p' : 'u';
    ext[4] = '\0';

    parallelfilecurrent = subst_star(base,update->ntimestep,padflag);
    delete [] base;
  }
}

/* ---------------------------------------------------------------------- */

std::string DumpSurfVTK::pvtk_piece_filename(int id)
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

void DumpSurfVTK::write_pvtk(int fileformat)
{
  FILE *fp = fopen(parallelfilecurrent,"w");
  if (!fp) error->one(FLERR,"Cannot open dump surf/vtk parallel file");

  const char *gridtype =
    (fileformat == PVTP) ? "PPolyData" : "PUnstructuredGrid";
  int one = 1;
  const char *byte_order =
    (*((char *) &one)) ? "LittleEndian" : "BigEndian";

  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"%s\" version=\"1.0\" byte_order=\"%s\">\n",
          gridtype,byte_order);
  fprintf(fp,"  <%s GhostLevel=\"0\">\n",gridtype);

  fprintf(fp,"    <PCellData>\n");
  for (int f = 0; f < (int) fields.size(); f++) {
    const char *type = (fields[f].type == DOUBLE) ? "Float64" : "Int64";
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

  fprintf(fp,"  </%s>\n",gridtype);
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}

/* ---------------------------------------------------------------------- */

int DumpSurfVTK::modify_param(int narg, char **arg)
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
