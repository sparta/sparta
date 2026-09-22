/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
------------------------------------------------------------------------- */

#ifndef SPARTA_NVTX_H
#define SPARTA_NVTX_H

#include "KOKKOS/nvtx3/nvToolsExt.h"

namespace SPARTA_NS {

class NvtxRegion {
 public:
  explicit NvtxRegion(const char *name) { nvtxRangePushA(name); }
  ~NvtxRegion() { nvtxRangePop(); }

 private:
  NvtxRegion(const NvtxRegion &);
  NvtxRegion &operator=(const NvtxRegion &);
};

}    // namespace SPARTA_NS

#define SPARTA_NVTX_JOIN_IMPL(a,b) a##b
#define SPARTA_NVTX_JOIN(a,b) SPARTA_NVTX_JOIN_IMPL(a,b)
#define SPARTA_NVTX_RANGE(name) \
  SPARTA_NS::NvtxRegion SPARTA_NVTX_JOIN(sparta_nvtx_region_,__LINE__)(name)

// Stringifying the destination gives each transfer a useful, stable name in
// Nsight Systems without requiring a hand-maintained label at every call site.
#define SPARTA_NVTX_DEEP_COPY(dst,...) do { \
  SPARTA_NS::NvtxRegion sparta_nvtx_deep_copy_region( \
    "Kokkos::deep_copy(" #dst ")"); \
  Kokkos::deep_copy(dst,__VA_ARGS__); \
} while (0)

#define SPARTA_NVTX_DEEP_COPY_ASYNC(exec,dst,...) do { \
  SPARTA_NS::NvtxRegion sparta_nvtx_deep_copy_region( \
    "Kokkos::deep_copy(" #dst ", async)"); \
  Kokkos::deep_copy(exec,dst,__VA_ARGS__); \
} while (0)

#endif
