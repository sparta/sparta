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

#include "string.h"
#include "stdlib.h"
#include "compute_fft_grid_kokkos.h"
#include "update.h"
#include "domain.h"
#include "grid_kokkos.h"
#include "modify.h"
#include "fix.h"
#include "input.h"
#include "kokkos.h"
#include "variable.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"
#include "comm.h"

#ifdef SPARTA_MAP
#include <map>
#elif defined SPARTA_UNORDERED_MAP
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

using namespace SPARTA_NS;

enum{COMPUTE,FIX,VARIABLE};
enum{FORWARD=1,BACWARD=-1};

#define INVOKED_PER_GRID 16

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

/* ---------------------------------------------------------------------- */

ComputeFFTGridKokkos::ComputeFFTGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeFFTGrid(sparta, narg, arg)
{
  kokkos_flag = 1;

  fft2dKK = NULL;
  fft3dKK = NULL;
  irregular1KK = irregular2KK = NULL;

#if defined (SPARTA_KOKKOS_GPU)
  #if defined(FFT_KOKKOS_KISS)
    if (comm->me == 0)
      error->warning(FLERR,"Using default KISS FFT with Kokkos GPU backends may give suboptimal performance");
  #endif
#endif
}

/* ---------------------------------------------------------------------- */

ComputeFFTGridKokkos::~ComputeFFTGridKokkos()
{
  if (copymode) return;

  if (dimension == 3) delete fft3dKK;
  else delete fft2dKK;

  memoryKK->destroy_kokkos(k_fftwork, fftwork);
  fftwork = NULL;

  memoryKK->destroy_kokkos(k_ingrid, ingrid);
  ingrid = NULL;

  memoryKK->destroy_kokkos(k_vector_grid, vector_grid);
  memoryKK->destroy_kokkos(k_array_grid, array_grid);
  vector_grid = NULL;
  array_grid = NULL;

  delete irregular1KK;
  delete irregular2KK;
}

/* ---------------------------------------------------------------------- */

void ComputeFFTGridKokkos::post_constructor()
{
  // partition for FFTs
  // allocate bufs for grid and FFT decomps
  // NOTE: could avoid allocating inbuf in some cases, depends on values

  copymode = 1;

  fft_create();

  MemKK::realloc_kokkos(d_fft, "fft/grid:fft", 2*nfft);
  d_fft_char = DAT::t_char_1d((char *)d_fft.data(),d_fft.size()*sizeof(FFT_SCALAR));

  memoryKK->create_kokkos(k_fftwork, fftwork, nfft, "fft/grid:fftwork");
  d_fftwork = k_fftwork.d_view;
  d_fftwork_char = DAT::t_char_1d((char *)d_fftwork.data(),d_fftwork.size()*sizeof(double));

  reallocate();

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

void ComputeFFTGridKokkos::compute_per_grid()
{
  if (sparta->kokkos->prewrap) {
    ComputeFFTGrid::compute_per_grid();
  } else {
    compute_per_grid_kokkos();
    if (ncol == 1) {
      k_vector_grid.modify_device();
      k_vector_grid.sync_host();
    } else {
      k_array_grid.modify_device();
      k_array_grid.sync_host();
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeFFTGridKokkos::compute_per_grid_kokkos()
{
  copymode = 1;

  invoked_per_grid = update->ntimestep;

  // check that grid has not adapted
  // NOTE: also need to check it has not been re-balanced?

  if (grid->maxlevel != 1)
    error->all(FLERR,"Compute fft/grid requires uniform one-level grid");

  // if sumflag set, zero output vector/array, but not K-space indices
  // so can sum each value's result into it

  if (sumflag) {
    if (ncol == 1)
      Kokkos::deep_copy(d_vector_grid,0.0);
    else
      Kokkos::deep_copy(d_array_grid,0.0);
  }

  // process values, one FFT per value

  for (int m = 0; m < nvalues; m++) {
    const int vidx = value2index[m];
    const int aidx = argindex[m];

    // invoke compute if not previously invoked
    // for per-grid compute, invoke post_process_grid() if necessary

    if (which[m] == COMPUTE) {
      Compute *c = modify->compute[vidx];

      if (!c->kokkos_flag)
        error->all(FLERR,"Cannot (yet) use non-Kokkos computes with compute fft/grid/kk");

      KokkosBase* cKKBase = dynamic_cast<KokkosBase*>(c);

      if (!(c->invoked_flag & INVOKED_PER_GRID)) {
        cKKBase->compute_per_grid_kokkos();
        c->invoked_flag |= INVOKED_PER_GRID;
      }

      if (c->post_process_grid_flag) {
        DAT::t_float_2d_lr d_tmp1;
        DAT::t_float_1d_strided d_tmp2;
        cKKBase->post_process_grid_kokkos(aidx,1,d_tmp1,NULL,d_tmp2);
      }

      if (aidx == 0 || c->post_process_grid_flag) {
        d_ingrid = cKKBase->d_vector_grid;
      } else {
        auto d_carray = cKKBase->d_array_grid;
        auto d_ingrid = k_ingrid.d_view;
        const int n = grid->nlocal;
        const int aidxm1 = aidx - 1;

        Kokkos::parallel_for(n, SPARTA_LAMBDA(int i) {
          d_ingrid[i] = d_carray(i,aidxm1);
        });
      }

    // access fix fields, check if fix frequency is a match

    } else if (which[m] == FIX) {
      Fix *fix = modify->fix[vidx];

      KokkosBase* fixKKBase = dynamic_cast<KokkosBase*>(fix);

      if (!fixKKBase || !fix->per_grid_flag)
        error->all(FLERR,"Unsupported fix used by compute fft/grid/kk");

      if (update->ntimestep % modify->fix[vidx]->per_grid_freq)
        error->all(FLERR,"Fix used in compute fft/grid not "
                   "computed at compatible time");

      if (aidx == 0) {
        d_ingrid = fixKKBase->d_vector_grid;
      } else {
        auto d_farray = fixKKBase->d_array_grid;
        auto d_ingrid = k_ingrid.d_view;
        const int n = grid->nlocal;
        const int aidxm1 = aidx - 1;

        Kokkos::parallel_for(n, SPARTA_LAMBDA(int i) {
          d_ingrid[i] = d_farray(i,aidxm1);
        });
      }

    // evaluate particle-style or grid-style variable

    } else if (which[m] == VARIABLE) {
      input->variable->compute_grid(vidx,ingrid,1,0);
      k_ingrid.modify_host();
      k_ingrid.sync_device();
    }

    // ------------------------------
    // perform FFT on inbufptr values
    // ------------------------------

    // irregular comm to move grid values from SPARTA owners -> FFT owners

    auto d_ingrid_char = DAT::t_char_1d((char *)d_ingrid.data(),d_ingrid.size()*sizeof(double));

    irregular1KK->exchange_uniform(d_ingrid_char,sizeof(double),(char *)d_fftwork_char.data(),
                                   d_fftwork_char);

    // convert SPARTA grid value to FFT complex
    // use map1 to morph recvbuf to local FFT layout

    Kokkos::parallel_for(nfft, SPARTA_CLASS_LAMBDA(int i) {
      const int n = 2*d_map1(i);
      d_fft[n] = d_fftwork[i];
      d_fft[n+1] = ZEROF;
    });

    // perform forward FFT

    if (dimension == 3) fft3dKK->compute(d_fft,d_fft,FORWARD);
    else fft2dKK->compute(d_fft,d_fft,FORWARD);

    // reverse irregular comm to move results from FFT grid -> SPARTA grid
    // if conjugate set:
    //   convert complex FFT datums back to real via c times c*
    //   comm single floating point value per grid cell
    // else: comm entire complex value
    // copy or sum received values into output vec or array via map2

    if (conjugate) {

      Kokkos::parallel_for(nfft, SPARTA_CLASS_LAMBDA(int i) {
        const FFT_SCALAR real = d_fft[2*i];
        const FFT_SCALAR imag = d_fft[2*i+1];
        d_fftwork[i] = real*real + imag*imag;
      });

      irregular2KK->exchange_uniform(d_fftwork_char,sizeof(double),
                                   (char *) d_gridwork_char.data(),d_gridwork_char);
      if (sumflag) {
        if (ncol == 1) {
          const int n = grid->nlocal;
          Kokkos::parallel_for(n, SPARTA_CLASS_LAMBDA(int i) {
            d_vector_grid[d_map2[i]] += d_gridwork[i];
          });
        } else {
          const int icol = startcol;
          const int n = grid->nlocal;
          Kokkos::parallel_for(n, SPARTA_CLASS_LAMBDA(int i) {
            d_array_grid(d_map2[i],icol) += d_gridwork[i];
          });
        }
      } else {
        if (ncol == 1) {
          const int n = grid->nlocal;
          Kokkos::parallel_for(n, SPARTA_CLASS_LAMBDA(int i) {
            d_vector_grid[d_map2[i]] = d_gridwork[i];
          });
        } else {
          const int icol = m + startcol;
          const int n = grid->nlocal;
          Kokkos::parallel_for(n, SPARTA_CLASS_LAMBDA(int i) {
            d_array_grid(d_map2[i],icol) = d_gridwork[i];
          });
        }
      }
    } else {
      irregular2KK->exchange_uniform(d_fft_char, 2*sizeof(FFT_SCALAR),
                                      (char *) d_gridworkcomplex_char.data(), d_gridworkcomplex_char);

      if (sumflag) {
        const int icol = startcol;
        const int n = grid->nlocal;
        Kokkos::parallel_for(n, SPARTA_CLASS_LAMBDA(int i) {
          d_array_grid(d_map2[i],icol) += d_gridworkcomplex[2*i];
          d_array_grid(d_map2[i],icol+1) += d_gridworkcomplex[2*i+1];
        });
      } else {
        const int icol = 2*m + startcol;
        const int n = grid->nlocal;
        Kokkos::parallel_for(n, SPARTA_CLASS_LAMBDA(int i) {
          d_array_grid(d_map2[i],icol) = d_gridworkcomplex[2*i];
          d_array_grid(d_map2[i],icol+1) = d_gridworkcomplex[2*i+1];
        });
      }
    }
  }

  // scale results if requested

  if (scalefactor == 1.0) return;

  if (ncol == 1) {
    Kokkos::parallel_for(nglocal, SPARTA_CLASS_LAMBDA(int i) {
      d_vector_grid[i] *= scalefactor;
    });
  } else {
    Kokkos::parallel_for(nglocal, SPARTA_CLASS_LAMBDA(int i) {
      for (int m = startcol; m < ncol; m++)
        d_array_grid(i,m) *= scalefactor;
    });
  }

  copymode = 0;
}

/* ----------------------------------------------------------------------
   reallocate irregular comm patterns if local grid storage changes
   called by init() and whenever grid is rebalanced
------------------------------------------------------------------------- */

void ComputeFFTGridKokkos::reallocate()
{
  copymode = 1;

  delete irregular1KK;
  delete irregular2KK;

  irregular_create();

  memoryKK->destroy_kokkos(k_ingrid, ingrid);
  memoryKK->destroy_kokkos(k_vector_grid, vector_grid);
  memoryKK->destroy_kokkos(k_array_grid, array_grid);

  nglocal = grid->nlocal;

  memoryKK->create_kokkos(k_ingrid,ingrid,nglocal,"fft/grid:ingrid");
  d_ingrid = k_ingrid.d_view;

  if (startcol || conjugate) {
    MemKK::realloc_kokkos(d_gridwork,"fft/grid:gridwork",nglocal);
    d_gridwork_char = DAT::t_char_1d((char *)d_gridwork.data(),d_gridwork.size()*sizeof(double));
  }

  if (!conjugate) {
    MemKK::realloc_kokkos(d_gridworkcomplex,"fft/grid:gridworkcomplex",2*nglocal);
    d_gridworkcomplex_char = DAT::t_char_1d((char *)d_gridworkcomplex.data(),d_gridworkcomplex.size()*sizeof(FFT_SCALAR));
  }

  if (ncol == 1) {
    memoryKK->create_kokkos(k_vector_grid,vector_grid,nglocal,"fft/grid:vector_grid");
    d_vector_grid = k_vector_grid.d_view;
  } else {
    memoryKK->create_kokkos(k_array_grid,array_grid,nglocal,ncol,"fft/grid:array_grid");
    d_array_grid = k_array_grid.d_view;
  }

  // one-time setup of vector of K-space vector magnitudes if requested
  // compute values in K-space, irregular comm to grid decomposition
  // kx,ky,kz = indices of FFT grid cell in K-space
  // convert to distance from (0,0,0) cell using PBC
  // klen = length of K-space vector

  if (!startcol) return;

  int i,j,k;
  double ikx,iky,ikz;
  double klen;

  int nxhalf = nx/2;
  int nyhalf = ny/2;
  int nzhalf = nz/2;
  if (dimension == 2) nzhalf = 1;

  int icol = 0;

  for (int m = 0; m < 4; m++) {
    if (m == 0 && !kx) continue;
    if (m == 1 && !ky) continue;
    if (m == 2 && !kz) continue;
    if (m == 3 && !kmag) continue;

    int n = 0;
    for (k = nzlo; k <= nzhi; k++) {
      if (k < nzhalf) ikz = k;
      else ikz = nz - k;
      for (j = nylo; j <= nyhi; j++) {
        if (j < nyhalf) iky = j;
        else iky = ny - j;
        for (i = nxlo; i <= nxhi; i++) {
          if (i < nxhalf) ikx = i;
          else ikx = nx - i;

          if (m == 0) klen = ikx;
          else if (m == 1) klen = iky;
          else if (m == 2) klen = ikz;
          else if (m == 3) klen = sqrt(ikx*ikx + iky*iky + ikz*ikz);

          fftwork[n++] = klen;
        }
      }
    }

    k_fftwork.modify_host();
    k_fftwork.sync_device();

    irregular2KK->exchange_uniform(d_fftwork_char,sizeof(double),
                                 (char *) d_gridwork_char.data(),d_gridwork_char);

    copymode = 1;

    Kokkos::parallel_for(nglocal, SPARTA_CLASS_LAMBDA(int i) {
      d_array_grid(d_map2[i],icol) = d_gridwork[i];
    });

    icol++;
  }

  copymode = 0;
}

/* ----------------------------------------------------------------------
   one-time creation and partitioning of the FFT grid
------------------------------------------------------------------------- */

void ComputeFFTGridKokkos::fft_create()
{
  // create FFT grid and partition it
  // proc partitioning = Npx by Npy by Npz
  // this proc's index within partitioning = xme,yme,zme
  // sub-domain this proc owns (inclusive): (Nxlo,Nxhi) (Nylo,Nyhi) (Nzlo,Nzhi)
  // nfft = # of FFT points on this proc

  nx = grid->unx;
  ny = grid->uny;
  nz = grid->unz;

  // warn if any grid dimension is not factorable by 2,3,5

  int flag = 0;
  if (!factorable(nx)) flag = 1;
  if (!factorable(ny)) flag = 1;
  if (dimension == 3 && !factorable(nz)) flag = 1;
  if (flag && me == 0)
    error->warning(FLERR,"Compute fft/grid FFT grid is not "
                   "factorable by 2,3,5");

  if (dimension == 3) {
    npx = 1;
    if (nz >= nprocs) {
      npy = 1;
      npz = nprocs;
    } else procs2grid2d(nprocs,ny,nz,npy,npz);
  } else {
    npx = 1;
    npy = nprocs;
    npz = 1;
  }

  int xme = me % npx;
  int yme = (me/npx) % npy;
  int zme = me / (npx*npy);

  nxlo = xme*nx/npx;
  nxhi = (xme+1)*nx/npx - 1;
  nylo = yme*ny/npy;
  nyhi = (yme+1)*ny/npy - 1;
  nzlo = zme*nz/npz;
  nzhi = (zme+1)*nz/npz - 1;

  nxfft = nxhi - nxlo + 1;
  nyfft = nyhi - nylo + 1;
  nzfft = nzhi - nzlo + 1;

  nfft = nxfft * nyfft * nzfft;

  //printf("FFT %d: nxyz %d %d %d np xyz %d %d %d: "
  //       "x %d %d y %d %d z %d %d: %d\n",
  //       me,nx,ny,nz,npx,npy,npz,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,nfft);

  bigint nfft2 = nfft;
  nfft2 *= 2;
  if (nfft2 > MAXSMALLINT)
    error->all(FLERR,"Compute fft/grid FFT is too large per-processor");

  // create FFT plan

  int collective_flag = 0; // not yet supported in Kokkos version
  int gpu_aware_flag = sparta->kokkos->gpu_aware_flag;

  int tmp;
  if (dimension == 3) {
    fft3dKK = new FFT3dKokkos<DeviceType>(sparta,world,nx,ny,nz,
                              nxlo,nxhi,nylo,nyhi,nzlo,nzhi,
                              nxlo,nxhi,nylo,nyhi,nzlo,nzhi,
                              0,0,&tmp,collective_flag,gpu_aware_flag);
  } else {
    fft2dKK = new FFT2dKokkos<DeviceType>(sparta,world,nx,ny,
                              nxlo,nxhi,nylo,nyhi,nxlo,nxhi,nylo,nyhi,
                              0,0,&tmp,collective_flag,gpu_aware_flag);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeFFTGridKokkos::irregular_create()
{
  copymode = 1;

  int *proclist1,*proclist2,*proclist3;
  DAT::tdual_int_1d k_proclist1,k_proclist2,k_proclist3;
  DAT::t_char_1d d_sbuf1,d_rbuf1,d_sbuf2,d_rbuf2;

  // plan for moving data from SPARTA grid -> FFT grid
  // send my cell IDs in SPARTA grid->cells order
  //   to proc who owns each in FFT partition
  // create map1 = how to reorder them on receiving FFT proc

  irregular1KK = new IrregularKokkos(sparta);

  GridKokkos* gridKK = (GridKokkos*) grid;
  int nglocal = grid->nlocal;

  t_cell_1d d_cells;
  if (sparta->kokkos->prewrap) {
    d_cells = t_cell_1d("grid:cells",nglocal);
    auto h_cells = Kokkos::create_mirror_view(d_cells);
    for (int i = 0; i < nglocal; i++)
      h_cells[i] = grid->cells[i];
    Kokkos::deep_copy(d_cells,h_cells);
  } else {
    d_cells = gridKK->k_cells.d_view;
    gridKK->sync(Device,CELL_MASK);
  }

  memoryKK->create_kokkos(k_proclist1,proclist1,nglocal,"fft/grid:proclist1");
  auto d_proclist1 = k_proclist1.d_view;

  // use cell ID to determine which proc owns it in FFT partitioning
  // while loops are for case where some procs own nothing in FFT partition,
  //   due to more procs that rows/columns of FFT partition
  //   in this case formula (iy/ny * npy) may undershoot correct ipy proc
  //   increment until find ipy consistent with
  //   ipy*ny/npy to (ipy+1)*ny/npy bounds of that proc's FFT partition

  Kokkos::parallel_for(nglocal, SPARTA_CLASS_LAMBDA(int i) {
    const cellint gid = d_cells[i].id;
    const int iy = ((gid-1) / nx) % ny;
    const int iz = (gid-1) / (nx*ny);

    int ipy = static_cast<int> (1.0*iy/ny * npy);
    while (1) {
      if (iy >= ipy*ny/npy && iy < (ipy+1)*ny/npy) break;
      ipy++;
    }
    int ipz = static_cast<int> (1.0*iz/nz * npz);
    while (1) {
      if (iz >= ipz*nz/npz && iz < (ipz+1)*nz/npz) break;
      ipz++;
    }

    d_proclist1[i] = ipz*npy + ipy;
  });

  k_proclist1.modify_device();
  k_proclist1.sync_host();

  int nrecv = irregular1KK->create_data_uniform(nglocal,proclist1);

  if (nrecv != nfft)
    error->one(FLERR,"Compute fft/grid FFT mapping is inconsistent");

  MemKK::realloc_kokkos(d_sbuf1,"fft/grid:sbuf1",(bigint)nglocal*sizeof(cellint));
  MemKK::realloc_kokkos(d_rbuf1,"fft/grid:rbuf1",(bigint)nfft*sizeof(cellint));

  DAT::t_cellint_1d d_idsend = DAT::t_cellint_1d((cellint *)d_sbuf1.data(),nglocal);

  Kokkos::parallel_for(nglocal, SPARTA_LAMBDA(int i) {
    d_idsend[i] = d_cells[i].id;
  });

  irregular1KK->exchange_uniform(d_sbuf1,sizeof(cellint),d_rbuf1.data(),d_rbuf1);

  MemKK::realloc_kokkos(d_map1,"fft/grid:map1",nfft);

  DAT::t_cellint_1d d_idrecv = DAT::t_cellint_1d((cellint *)d_rbuf1.data(),nfft);

  Kokkos::parallel_for(nfft, SPARTA_CLASS_LAMBDA(int i) {
    const cellint gid = d_idrecv[i];
    const int ix = (gid-1) % nx;
    const int iy = ((gid-1) / nx) % ny;
    const int iz = (gid-1) / (nx*ny);
    d_map1[i] = (iz-nzlo)*nxfft*nyfft + (iy-nylo)*nxfft + (ix-nxlo);
  });

  // plan for moving data from FFT grid -> SPARTA grid
  // send my cell IDs in FFT grid order
  //   back to proc who owns each in SPARTA grid->cells
  // proclist3 generated from irregular1->reverse()
  //   must permute proclist3 -> proclist2 via map1
  //   same way that received grid data will be permuted by map1 to FFT data
  // create map2 = how to reorder received FFT data on SPARTA grid proc

  irregular2KK = new IrregularKokkos(sparta);

  memoryKK->create_kokkos(k_proclist3,proclist3,nfft,"fft/grid:proclist3");
  auto d_proclist3 = k_proclist3.d_view;
  irregular1KK->reverse(nrecv,proclist3);

  k_proclist3.modify_host();
  k_proclist3.sync_device();

  memoryKK->create_kokkos(k_proclist2,proclist2,nfft,"fft/grid:proclist2");
  auto d_proclist2 = k_proclist2.d_view;

  Kokkos::parallel_for(nfft, SPARTA_CLASS_LAMBDA(int i) {
    d_proclist2[d_map1[i]] = d_proclist3[i];
  });

  k_proclist2.modify_device();
  k_proclist2.sync_host();

  nrecv = irregular2KK->create_data_uniform(nfft,proclist2);
  if (nrecv != nglocal)
    error->one(FLERR,"Compute fft/grid FFT mapping is inconsistent");

  MemKK::realloc_kokkos(d_sbuf2,"fft/grid:sbuf2",(bigint)nfft*sizeof(cellint));
  MemKK::realloc_kokkos(d_rbuf2,"fft/grid:rbuf2",(bigint)nglocal*sizeof(cellint));

  d_idsend = DAT::t_cellint_1d((cellint *)d_sbuf2.data(),nfft);

  Kokkos::parallel_for(nfft, SPARTA_CLASS_LAMBDA(int i) {
    d_idsend[d_map1[i]] = d_idrecv[i];
  });

  irregular2KK->exchange_uniform(d_sbuf2,sizeof(cellint),d_rbuf2.data(),d_rbuf2);

  d_idrecv = DAT::t_cellint_1d((cellint *)d_rbuf2.data(),nglocal);

  MemKK::realloc_kokkos(d_map2,"fft/grid:map2",nglocal);

  // insure grid cell IDs are hashed, so can use them to build map2

  if (!grid->hashfilled) grid->rehash();

  gridKK->update_hash();
  auto hash_kk = gridKK->hash_kk;

  Kokkos::parallel_for(nglocal, SPARTA_CLASS_LAMBDA(int i) {
    const cellint gid = d_idrecv[i];
    const auto index = hash_kk.find(gid);
    d_map2[i] = hash_kk.value_at(index);
  });

  // clean up

  memoryKK->destroy_kokkos(k_proclist1,proclist1);
  memoryKK->destroy_kokkos(k_proclist2,proclist2);
  memoryKK->destroy_kokkos(k_proclist3,proclist3);

  copymode = 0;
}

/* ----------------------------------------------------------------------
   print out FFT precision and library
------------------------------------------------------------------------- */

void ComputeFFTGridKokkos::print_FFT_info()
{
  if (comm->me == 0) {
    char str[64];
    sprintf(str,"Using " SPARTA_FFT_PREC " precision " SPARTA_FFT_KOKKOS_LIB " for FFTs\n");
    if (screen) fprintf(screen,"%s",str);
    if (logfile) fprintf(logfile,"%s",str);
  }
}

