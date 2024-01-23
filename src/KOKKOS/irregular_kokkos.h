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

#ifndef SPARTA_IRREGULAR_KOKKOS_H
#define SPARTA_IRREGULAR_KOKKOS_H

#include "irregular.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

struct TagIrregularPackBuffer{};

struct TagIrregularUnpackBufferSelf{};

class IrregularKokkos : public Irregular {
 public:

  IrregularKokkos(class SPARTA *);
  ~IrregularKokkos();
  int create_data_uniform(int, int *, int sort = 0);
  int augment_data_uniform(int, int *);
  void exchange_uniform(DAT::t_char_1d, int, char *, DAT::t_char_1d);

  KOKKOS_INLINE_FUNCTION
  void operator()(TagIrregularPackBuffer, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagIrregularUnpackBufferSelf, const int&) const;

  inline
  void pack_buffer_serial(const int, const int) const;

 private:
  int offset_send;

  DAT::tdual_int_1d k_index_send;
  DAT::t_int_1d d_index_send;
  DAT::tdual_int_1d k_index_self;
  DAT::t_int_1d d_index_self;

  DAT::t_char_1d d_sendbuf;
  DAT::t_char_1d d_recvbuf;
  DAT::t_char_1d d_buf;
  HAT::t_char_1d h_recvbuf;
  HAT::t_char_1d h_buf;
  int nbytes;
};

}

#endif

/* ERROR/WARNING messages:

*/
