#include "kokkos_type.h"

namespace SPARTA_NS {

template <typename Device>
struct ExclScan {
  using value_type = long;
  using view_type = Kokkos::View<int*, Device>;
  using total_type = Kokkos::View<int, Device>;
  KOKKOS_INLINE_FUNCTION void init(value_type& update) const { update = 0; }
  KOKKOS_INLINE_FUNCTION void join(
      value_type& update, const value_type& input) const {
    update = update + input;
  }
  view_type in_;
  view_type out_;
  total_type total_;
  ExclScan(view_type &in, view_type &out, total_type &total):
    in_(in), out_(out), total_(total) {}
  KOKKOS_INLINE_FUNCTION void operator()(int i, value_type& update, bool final_pass) const {
    update += in_[i];
    if (final_pass) {
      out_[0] = 0;
      out_[i + 1] = static_cast<int>(update);
      if (i + 1 == int(in_.extent(0))) total_() = static_cast<int>(update);
    }
  }
  using execution_space = Device;
};

template <typename Device>
Kokkos::View<int*, Device> offset_scan(Kokkos::View<int*, Device> in, int& total) {

  Kokkos::View<int*, Device> out;

  if (in.size() == 0) {
    total = 0;
    out = Kokkos::View<int*, Device>(in.label() + "_scan",1);
  } else {
    out = Kokkos::View<int*, Device>(Kokkos::view_alloc(in.label() + "_scan",Kokkos::WithoutInitializing), in.size() + 1);
    Kokkos::View<int, Device> total_dev(Kokkos::view_alloc("scan_total",Kokkos::WithoutInitializing));
    typename Kokkos::View<int, Device>::HostMirror total_host(Kokkos::view_alloc("scan_total_mirror",Kokkos::WithoutInitializing));
    Kokkos::parallel_scan(in.size(), ExclScan<Device>(in, out, total_dev));
    Kokkos::deep_copy(total_host, total_dev);
    total = total_host();
  }
  return out;
}

template Kokkos::View<int*, SPAHostType> offset_scan(
    Kokkos::View<int*, SPAHostType> in, int& total);
#ifdef SPARTA_KOKKOS_GPU
template Kokkos::View<int*, DeviceType> offset_scan(
    Kokkos::View<int*, DeviceType> in, int& total);
#endif

}
