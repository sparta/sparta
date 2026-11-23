4. Packages
===========

This section gives an overview of the optional packages that extend
SPARTA functionality with instructions on how to build SPARTA with
each of them.  Packages are groups of files that enable a specific set
of features.  For example, the KOKKOS package provides styles that
can run on different hardware such as GPUs.  You can see the list of all
packages and "make" commands to manage them by typing "make package"
from within the src directory of the SPARTA distribution or
"cmake -DSPARTA\_LIST\_PKGS" from within a build directory.  `Section 2.3 <Section_start.html#start_3>`_ gives general info on how to install
and un-install packages as part of the SPARTA build process.

Packages may require some
additional code compiled located in the lib folder, or may require
an external library to be downloaded, compiled, installed, and SPARTA
configured to know about its location and additional compiler flags.

Following the next two tables is a sub-section for each package.  It
lists authors (if applicable) and summarizes the package contents.  It
has specific instructions on how to install the package, including (if
necessary) downloading or building any extra library it requires. It
also gives links to documentation, example scripts, and
pictures/movies (if available) that illustrate use of the package.

NOTE: To see the complete list of commands a package adds to SPARTA,
just look at the files in its src directory, e.g. "ls src/KOKKOS".
Files with names that start with fix, compute, etc correspond to
commands with the same style names.

In these two tables, the "Example" column is a sub-directory in the
examples directory of the distribution which has an input script that
uses the package.  E.g. "fft" refers to the examples/fft
directory; The "Library" column indicates whether an extra library is needed to build
and use the package:

* dash = no library
* sys = system library: you likely have it on your machine
* int = internal library: provided with SPARTA, but you may need to build it
* ext = external library: you will need to download and install it on your machine


----------


.. _pkg_1:

.. raw:: html

   <span id="pkg_1"></span>

**SPARTA packages**

+---------------------+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------+------------+
| Package             | Description                   | Doc page                                                  | Example                                             | Library    |
+---------------------+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------+------------+
| `FFT <#FFT>`_       | fast Fourier transforms       | :doc:`compute\_style compute/fft/grid <compute_fft_grid>` | fft                                                 | int or ext |
+---------------------+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------+------------+
| `KOKKOS <#KOKKOS>`_ | Kokkos-enabled styles         | `Section 5.3 <Section_accelerate.html#acc_3>`_            | `Benchmarks <https://sparta.github.io/bench.html>`_ |\-          |
+---------------------+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------+------------+
| `VTK <#VTK>`_       | native VTK-format dump output | :doc:`dump particle/vtk <dump_vtk>`                       |\-                                                   | ext        |
+---------------------+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------+------------+


----------


.. raw:: html

   <span id="FFT"></span>

FFT package
-----------

**Contents:**

Apply Fast Fourier Transforms (FFTs) to simulation data. The FFT
library is specified in the Makefile.machine or CMake using the
FFT\_INC, FFT\_PATH, and FFT\_LIB variables. Supported external FFT
libraries that can be specified include FFTW3 or MKL. If no FFT
library is specified, SPARTA will use the internal KISS FFT library
that is included with SPARTA.

Similarly an external FFT library can be specified for the KOKKOS
package.  Options are CUFFT, HIPFFT, FFTW3, MKL, or MKL\_GPU. If no FFT
library is specified in CMake, SPARTA will use the internal
Kokkos version of the KISS FFT library that is included with SPARTA.

See the see discussion in `Section 2.2 <Section_start.html#start_2>`_ (step 6).

**Install or un-install with make:**


.. parsed-literal::

   make yes-fft
   make machine

   make no-fft
   make machine

**Install or un-install with CMake:**


.. parsed-literal::

   cd build
   cmake -C /path/to/sparta/cmake/presets/machine.cmake -DPKG_FFT=ON /path/to/sparta/cmake
   make

   cmake -C /path/to/sparta/cmake/presets/machine.cmake -DPKG_FFT=OFF /path/to/sparta/cmake
   make

**Supporting info:**

* :doc:`compute fft/grid <compute_fft_grid>`
* examples/fft


----------


.. raw:: html

   <span id="KOKKOS"></span>

KOKKOS package
--------------

**Contents:**

Styles adapted to compile using the Kokkos library which can convert
them to OpenMP, CUDA, HIP, or SYCL code so that they run efficiently on
multicore CPUs, NVIDIA GPUs, AMD GPUs, or Intel
GPUs.  All the styles have a "kk" as a suffix in their style name.
`Section 5.3 <Section_accelerate.html#acc_3>`_ gives details of what
hardware and software is required on your system, and how to build and
use this package.  Its styles can be invoked at run time via the "-sf
kk" or "-suffix kk" `command-line switches <Section_start.html#start_7>`_.

You must have a C++20 compatible compiler to use this package.

IMPORTANT NOTE: The KOKKOS package must be built using CMake. GNU
Makefile builds are not supported.

**Authors:** The KOKKOS package was created primarily by Stan Moore (Sandia),
with contributions from other folks as well.
It uses the open-source `Kokkos library <https://github.com/kokkos>`_
which was developed by Carter Edwards, Christian Trott, and others at
Sandia, and which is included in the SPARTA distribution in
lib/kokkos.

**Install or un-install:**

The KOKKOS package is built with CMake using a preset file that
specifies the Kokkos backend (OpenMP, CUDA, HIP, or SYCL) and the
target hardware architecture. Preset files are in cmake/presets/.

NOTE: You cannot build one executable that runs on multiple hardware
targets. Build SPARTA separately for each target.

**For multicore CPUs using OpenMP:**


.. parsed-literal::

   cmake -C /path/to/sparta/cmake/presets/kokkos_omp.cmake /path/to/sparta/cmake
   make

To select a specific CPU architecture (e.g. Haswell):


.. parsed-literal::

   cmake -C /path/to/sparta/cmake/presets/kokkos_omp.cmake   -DKokkos_ARCH_HSW=ON /path/to/sparta/cmake
   make

**For MPI-only (no threading):**


.. parsed-literal::

   cmake -C /path/to/sparta/cmake/presets/kokkos_mpi_only.cmake /path/to/sparta/cmake
   make

**For NVIDIA GPUs using CUDA:**


.. parsed-literal::

   cmake -C /path/to/sparta/cmake/presets/kokkos_cuda.cmake /path/to/sparta/cmake
   make

The kokkos\_cuda.cmake preset defaults to Hopper (H100) GPUs. To target
a different architecture, override the arch flag:

For A100 GPUs (Ampere):


.. parsed-literal::

   cmake -C /path/to/sparta/cmake/presets/kokkos_cuda.cmake   -DKokkos_ARCH_HOPPER90=OFF -DKokkos_ARCH_AMPERE80=ON   /path/to/sparta/cmake
   make

For V100 GPUs (Volta):


.. parsed-literal::

   cmake -C /path/to/sparta/cmake/presets/kokkos_cuda.cmake   -DKokkos_ARCH_HOPPER90=OFF -DKokkos_ARCH_VOLTA70=ON   /path/to/sparta/cmake
   make

**For AMD GPUs using HIP:**


.. parsed-literal::

   cmake -C /path/to/sparta/cmake/presets/kokkos_hip.cmake /path/to/sparta/cmake
   make

The kokkos\_hip.cmake preset defaults to MI250X GPUs. To target
MI300X/MI300A GPUs:


.. parsed-literal::

   cmake -C /path/to/sparta/cmake/presets/kokkos_hip.cmake   -DKokkos_ARCH_VEGA90A=OFF -DKokkos_ARCH_AMD_GFX942=ON   /path/to/sparta/cmake
   make

**For Intel GPUs using SYCL:**


.. parsed-literal::

   cmake -C /path/to/sparta/cmake/presets/kokkos_sycl.cmake /path/to/sparta/cmake
   make

The GPU architecture is detected automatically at runtime; no arch
override is needed for Intel GPUs.

To uninstall (disable) the KOKKOS package:


.. parsed-literal::

   cmake -C /path/to/sparta/cmake/presets/machine.cmake -DPKG_KOKKOS=OFF /path/to/sparta/cmake
   make

**Supporting info:**

* src/KOKKOS: filenames -> commands
* src/KOKKOS/README
* lib/kokkos/README
* the `Accelerating SPARTA <Section_accelerate.html#acc_3>`_ section
* `Section 5.3 <Section_accelerate.html#acc_3>`_
* `Section 2.6 -k on ... <Section_start.html#start_7>`_
* `Section 2.6 -sf kk <Section_start.html#start_7>`_
* `Section 2.6 -pk kokkos <Section_start.html#start_7>`_
* :doc:`package kokkos <package>`
* `Benchmarks page <https://sparta.github.io/bench.html>`_ of web site


----------


.. raw:: html

   <span id="VTK"></span>

VTK package
-----------

**Contents:**

Dump styles that write native VTK files (legacy .vtk, XML .vtp/.vtu,
and their parallel .pvtp/.pvtu variants) for direct visualization of
SPARTA data in ParaView or VisIt. The package adds three dump styles:
"particle/vtk", "grid/vtk", and "surf/vtk", which write particles as
points, grid cells as voxels/pixels, and surface elements as
triangles/lines, respectively, with requested per-particle, per-grid,
or per-surf attributes attached as data arrays.

This package requires the external `VTK library <https://vtk.org>`_ from
Kitware (version 9 or later recommended). It is distinct from the
"tools/paraview" Python converters, which post-process existing SPARTA
dump files rather than writing VTK directly.

IMPORTANT NOTE: The VTK package must be built using CMake. GNU Makefile
builds are not supported.

**Install or un-install with CMake:**

Set VTK\_ROOT in your environment (or pass -DVTK\_ROOT=/path) so CMake
can locate your VTK install, then enable both the package and its TPL:


.. parsed-literal::

   cd build
   cmake -C /path/to/sparta/cmake/presets/machine.cmake   -DPKG_VTK=ON -DBUILD_VTK=ON /path/to/sparta/cmake
   make

To uninstall (disable) the VTK package:


.. parsed-literal::

   cmake -C /path/to/sparta/cmake/presets/machine.cmake   -DPKG_VTK=OFF -DBUILD_VTK=OFF /path/to/sparta/cmake
   make

**Supporting info:**

* src/VTK: filenames -> commands
* :doc:`dump particle/vtk, dump grid/vtk, dump surf/vtk <dump_vtk>`


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
