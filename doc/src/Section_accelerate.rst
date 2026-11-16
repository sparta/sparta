5. Accelerating SPARTA performance
==================================

This section describes various methods for improving SPARTA
performance for different classes of problems running on different
kinds of machines.

Currently the only option is to use the KOKKOS accelerator
packages provided with SPARTA that
contains code optimized for certain kinds of hardware, including
multi-core CPUs and GPUs.

* 5.1 `Measuring performance <#acc_1>`_
* 5.2 `Accelerator packages with optimized styles <#acc_2>`_
* 5.3 `KOKKOS package <#acc_3>`_

The `Benchmark page <https://sparta.github.io/bench.html>`_ of the SPARTA
web site gives performance results for the various accelerator
packages discussed in Section 5.2, for several of the standard SPARTA
benchmark problems, as a function of problem size and number of
compute nodes, on different hardware platforms.


----------


.. _acc_1:

.. raw:: html

   <span id="acc_1"></span>

5.1 Measuring performance
---------------------------------------------------------------------------------

Before trying to make your simulation run faster, you should
understand how it currently performs and where the bottlenecks are.

The best way to do this is run the your system (actual number of
particles) for a modest number of timesteps (say 100 steps) on several
different processor counts, including a single processor if possible.
Do this for an equilibrium version of your system, so that the
100-step timings are representative of a much longer run.  There is
typically no need to run for 1000s of timesteps to get accurate
timings; you can simply extrapolate from short runs.

For the set of runs, look at the timing data printed to the screen and
log file at the end of each SPARTA run.  `This section <Section_start.html#start_8>`_ of the manual has an overview.

Running on one (or a few processors) should give a good estimate of
the serial performance and what portions of the timestep are taking
the most time.  Running the same problem on a few different processor
counts should give an estimate of parallel scalability.  I.e. if the
simulation runs 16x faster on 16 processors, its 100% parallel
efficient; if it runs 8x faster on 16 processors, it's 50% efficient.

The most important data to look at in the timing info is the timing
breakdown and relative percentages.  For example, trying different
options for speeding up the FFTs will have little impact
if they only consume 10% of the run time.  If the collide time is
dominating, you may want to look at the KOKKOS package, as discussed
below.  Comparing how the percentages change as
you increase the processor count gives you a sense of how different
operations within the timestep are scaling.

Another important detail in the timing info are the histograms of
particles counts and neighbor counts.  If these vary widely across
processors, you have a load-imbalance issue.  This often results in
inaccurate relative timing data, because processors have to wait when
communication occurs for other processors to catch up.  Thus the
reported times for "Communication" or "Other" may be higher than they
really are, due to load-imbalance.  If this is an issue, you can
uncomment the MPI\_Barrier() lines in src/timer.cpp, and recompile
SPARTA, to obtain synchronized timings.


----------


.. _acc_2:

.. raw:: html

   <span id="acc_2"></span>

5.2 Packages with optimized styles
------------------------------------------------------------------------------------------

Accelerated versions of various :doc:`collide\_style <collide>`,
:doc:`fixes <fix>`, :doc:`computes <compute>`, and other commands have
been added to SPARTA via the KOKKOS package, which may run faster than
the standard non-accelerated versions.

All of these commands are in the KOKKOS package provided with SPARTA.
An overview of packages is give in :doc:`Section packages <Section_packages>`.

SPARTA currently has acceleration support for different kinds of hardware
via the KOKKOS package: many-core CPUs, NVIDIA GPUs, AMD GPUs, and
Intel GPUs.

Whether you will see speedup for your hardware may depend on the size
problem you are running and what commands (accelerated and
non-accelerated) are invoked by your input script.  While these doc
pages include performance guidelines, there is no substitute for
trying out the KOKKOS package.

Any accelerated style has the same name as the corresponding standard
style, except that a suffix is appended.  Otherwise, the syntax for
the command that uses the style is identical, their functionality is
the same, and the numerical results it produces should also be the
same, except for precision and round-off effects, and differences in
random numbers.

For example, the KOKKOS package provides an accelerated variant of the
Temperature Compute :doc:`compute temp <compute_temp>`, namely :doc:`compute temp/kk <compute_temp>`

To see what accelerate styles are currently available, see `Section 3.5 <Section_commands.html#cmd_5>`_ of the manual.  The doc pages for
individual commands (e.g. :doc:`compute temp <compute_temp>`) also list
any accelerated variants available for that style.

To use an accelerator package in SPARTA, and one or more of the styles
it provides, follow these general steps:

using CMake from a build directory:

+---------------------------------+---------------------------------------------------------------------------------------+
| install the accelerator package | cmake -DPKG\_FFT=ON -DPKG\_KOKKOS=ON, etc                                             |
+---------------------------------+---------------------------------------------------------------------------------------+
| add compile/link flags          | cmake -C /path/to/sparta/cmake/presets/kokkos\_cuda.cmake -DKokkos\_ARCH\_PASCAL60=ON |
+---------------------------------+---------------------------------------------------------------------------------------+
| re-build SPARTA                 | make                                                                                  |
+---------------------------------+---------------------------------------------------------------------------------------+

Then do the following:

+-------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------+
| prepare and test a regular SPARTA simulation                                                                                              | lmp\_kokkos\_cuda -in in.script; mpirun -np 32 lmp\_kokkos\_cuda -in in.script |
+-------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------+
| enable specific accelerator support via '-k on' `command-line switch <Section_start.html#start_7>`_,                                      | -k on g 1                                                                      |
+-------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------+
| set any needed options for the package via "-pk" `command-line switch <Section_start.html#start_7>`_ or :doc:`package <package>` command, | only if defaults need to be changed, -pk kokkos react/retry yes                |
+-------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------+
| use accelerated styles in your input via "-sf" `command-line switch <Section_start.html#start_7>`_ or :doc:`suffix <suffix>` command      | lmp\_kokkos\_cuda -in in.script -sf kk                                         |
+-------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------+

Note that the first 3 steps can be done as a single command with
suitable make command invocations. This is discussed in :doc:`Section 4 <Section_packages>` of the manual, and its use is illustrated in
the individual accelerator sections.  Typically these steps only need
to be done once, to create an executable that uses one or more
accelerator packages.

The last 4 steps can all be done from the command-line when SPARTA is
launched, without changing your input script, as illustrated in the
individual accelerator sections.  Or you can add
:doc:`package <package>` and :doc:`suffix <suffix>` commands to your input
script.

The `Benchmark page <https://sparta.github.io/bench.html>`_ of the SPARTA
web site gives performance results for the various accelerator
packages for several of the standard SPARTA benchmark problems, as a
function of problem size and number of compute nodes, on different
hardware platforms.

Here is a brief summary of what the KOKKOS package provides.

* Styles with a "kk" suffix are part of the KOKKOS package, and can be
  run using OpenMP on multicore CPUs, on an NVIDIA GPU, on an AMD GPU,
  or on an Intel GPU.  The speed-up depends on a variety of
  factors, as discussed on the KOKKOS accelerator page.


The KOKKOS accelerator package doc page explains:

* what hardware and software the accelerated package requires
* how to build SPARTA with the accelerated package
* how to run with the accelerated package either via command-line switches or modifying the input script
* speed-ups to expect
* guidelines for best performance
* restrictions


----------


.. _acc_3:

.. raw:: html

   <span id="acc_3"></span>

5.3 KOKKOS package
--------------------------------------------------------------------------

Kokkos is a templated C++ library that provides abstractions to allow
a single implementation of an application kernel (e.g. a collision
style) to run efficiently on different kinds of hardware, such as
GPUs or many-core CPUs. Kokkos maps the C++ kernel
onto different backend languages such as CUDA, OpenMP, or Pthreads.
The Kokkos library also provides data abstractions to adjust (at
compile time) the memory layout of data structures like 2d and 3d
arrays to optimize performance on different hardware. For more
information on Kokkos, see
`Github <https://github.com/kokkos/kokkos>`_. Kokkos is part of
`Trilinos <https://github.com/trilinos/Trilinos>`_. The Kokkos
library was written primarily by Carter Edwards, Christian Trott, and
Dan Sunderland (all Sandia).

The SPARTA KOKKOS package contains versions of collide, fix, and
compute styles that use data structures and macros provided by the
Kokkos library, which is included with SPARTA in /lib/kokkos. The
KOKKOS package was developed primarily by Stan Moore (Sandia) with
contributions of various styles by others, including Dan Ibanez
(Sandia), Tim Fuller (Sandia), and Sam Mish (Sandia). For more
information on developing using Kokkos abstractions see the Kokkos
programmers' guide at /lib/kokkos/doc/Kokkos\_PG.pdf.

The KOKKOS package supports multiple execution backends (per MPI task):
Serial (MPI-only for CPUs), OpenMP (threading for many-core CPUs),
CUDA (NVIDIA GPUs), HIP (AMD GPUs), and SYCL (Intel GPUs
and other SYCL-capable devices). You choose the backend at build time
to produce an executable compatible with specific hardware.

NOTE: The KOKKOS package must be built using CMake. GNU Makefile builds
are not supported for the KOKKOS package.

NOTE: Kokkos support within SPARTA must be built with a C++20
compatible compiler. For a list of compilers that have been tested with
the Kokkos library, see the Kokkos `README <https://github.com/kokkos/kokkos/blob/master/README.md>`_.

**Building SPARTA with the KOKKOS package with CMake:**

To build with the KOKKOS package, start with one of the provided preset
files in /cmake/presets/. Preset files encode all the settings needed
for a particular hardware target. The available KOKKOS presets are:

* kokkos\_mpi\_only.cmake - Serial backend (MPI-only, no threading)
* kokkos\_omp.cmake      - OpenMP backend (multi-core CPUs)
* kokkos\_cuda.cmake     - CUDA backend (NVIDIA Hopper GPUs, generic OpenMPI/MPICH)
* kokkos\_hip.cmake      - HIP backend for AMD MI250X GPUs
* elcapitan\_kokkos.cmake - HIP backend for AMD MI300A APU (Cray MPICH)
* kokkos\_sycl.cmake     - SYCL backend for Intel Ponte Vecchio GPUs

You may need to override -D Kokkos\_ARCH\_\ *TYPE*\ =ON to match your
specific hardware. For example:

* for Sandy Bridge CPUs, set -D Kokkos\_ARCH\_SNB=ON
* for Broadwell CPUs, set -D Kokkos\_ARCH\_BDW=ON
* for V100 GPUs and Power9 CPUs, set -D Kokkos\_ARCH\_VOLTA70=ON -D Kokkos\_ARCH\_POWER9=ON
* for A100 GPUs, set -D Kokkos\_ARCH\_AMPERE80=ON
* for H100 GPUs, set -D Kokkos\_ARCH\_HOPPER90=ON

See the **Advanced Kokkos Options** section below for a complete listing
of all Kokkos architecture options.

NOTE: If you are migrating from GNU Makefile builds, see the table of
preset files above and pick the one that matches your hardware. The
Kokkos KOKKOS\_DEVICES and KOKKOS\_ARCH Makefile variables map to CMake
options as follows: KOKKOS\_DEVICES=OpenMP becomes
-DKokkos\_ENABLE\_OPENMP=ON, KOKKOS\_DEVICES=Cuda becomes
-DKokkos\_ENABLE\_CUDA=ON, KOKKOS\_ARCH=Volta70 becomes
-DKokkos\_ARCH\_VOLTA70=ON, etc.

**Compile for CPU-only (MPI only, no threading):**

Use a C++20 compatible compiler. Then do the following:


.. parsed-literal::

   mkdir build
   cd build
   cmake -C /path/to/sparta/cmake/presets/kokkos_mpi_only.cmake /path/to/sparta/cmake
   make -j 4

The resulting executable will be named spa\_kokkos\_mpi\_only.

**Compile for CPU-only (MPI plus OpenMP threading):**

NOTE: To build with Kokkos support for OpenMP threading, your compiler
must support the OpenMP interface. You should have one or more
multi-core CPUs so that multiple threads can be launched by each MPI
task running on a CPU.


.. parsed-literal::

   mkdir build
   cd build
   cmake -C /path/to/sparta/cmake/presets/kokkos_omp.cmake /path/to/sparta/cmake
   make -j 4

The resulting executable will be named spa\_kokkos\_omp.

To select a specific CPU architecture (e.g. Haswell), add:


.. parsed-literal::

   cmake -C /path/to/sparta/cmake/presets/kokkos_omp.cmake -DKokkos_ARCH_HSW=ON /path/to/sparta/cmake

**Compile for NVIDIA GPUs using CUDA (with OpenMPI or MPICH):**

NOTE: To build with Kokkos support for NVIDIA GPUs, NVIDIA CUDA
software version 12.2 or later must be installed on your system.

The kokkos\_cuda.cmake preset defaults to NVIDIA Hopper (H100) GPUs.
To target a different GPU, override the architecture flag. Common NVIDIA
GPU models and their Kokkos arch flags are:

* V100 (Volta)           ->  -DKokkos\_ARCH\_VOLTA70=ON
* A100 (Ampere)          ->  -DKokkos\_ARCH\_AMPERE80=ON
* H100/H200 (Hopper)     ->  -DKokkos\_ARCH\_HOPPER90=ON (default in kokkos\_cuda.cmake)
* L40S/RTX 4090 (Ada)    ->  -DKokkos\_ARCH\_ADA89=ON
* GB200 (Blackwell)      ->  -DKokkos\_ARCH\_BLACKWELL100=ON

Build for the default Hopper (H100) GPU:


.. parsed-literal::

   mkdir build
   cd build
   cmake -C /path/to/sparta/cmake/presets/kokkos_cuda.cmake /path/to/sparta/cmake
   make -j 4

Build for A100 GPUs (override the default Hopper arch):


.. parsed-literal::

   mkdir build
   cd build
   cmake -C /path/to/sparta/cmake/presets/kokkos_cuda.cmake   -DKokkos_ARCH_HOPPER90=OFF -DKokkos_ARCH_AMPERE80=ON   /path/to/sparta/cmake
   make -j 4

Build for V100 GPUs:


.. parsed-literal::

   mkdir build
   cd build
   cmake -C /path/to/sparta/cmake/presets/kokkos_cuda.cmake   -DKokkos_ARCH_HOPPER90=OFF -DKokkos_ARCH_VOLTA70=ON   /path/to/sparta/cmake
   make -j 4

Build with CUDA unified memory, so that a simulation can use memory on
the host to supplement the memory on the GPU, and host code can read
and write GPU data directly. See the **Advanced Kokkos Options** section
below for what this option requires and costs:


.. parsed-literal::

   mkdir build
   cd build
   cmake -C /path/to/sparta/cmake/presets/kokkos_cuda.cmake   -DKokkos_ENABLE_IMPL_CUDA_UNIFIED_MEMORY=ON   /path/to/sparta/cmake
   make -j 4

The resulting executable will be named spa\_kokkos\_cuda.

**Compile for AMD GPUs using HIP:**

NOTE: To build with Kokkos support for AMD GPUs, ROCm software version
6.2 or later must be installed on your system.

The kokkos\_hip.cmake preset targets AMD MI250X GPUs (GFX90A). Common AMD
GPU models and their Kokkos arch flags are:

* MI250X/MI250/MI210      ->  -DKokkos\_ARCH\_VEGA90A=ON (default in kokkos\_hip.cmake)
* MI300X/MI300A           ->  -DKokkos\_ARCH\_AMD\_GFX942=ON
* RX 7900 XTX (RDNA3)    ->  -DKokkos\_ARCH\_AMD\_GFX1100=ON

Build for the default MI250X GPU:


.. parsed-literal::

   mkdir build
   cd build
   cmake -C /path/to/sparta/cmake/presets/kokkos_hip.cmake /path/to/sparta/cmake
   make -j 4

Build for MI300X GPUs (override the default MI250X arch):


.. parsed-literal::

   mkdir build
   cd build
   cmake -C /path/to/sparta/cmake/presets/kokkos_hip.cmake   -DKokkos_ARCH_VEGA90A=OFF -DKokkos_ARCH_AMD_GFX942=ON   /path/to/sparta/cmake
   make -j 4

For AMD MI300A APU systems using Cray MPICH (e.g. El Capitan):


.. parsed-literal::

   mkdir build
   cd build
   cmake -C /path/to/sparta/cmake/presets/elcapitan_kokkos.cmake /path/to/sparta/cmake
   make -j 4

The resulting executable will be named spa\_kokkos\_hip or
spa\_elcapitan\_kokkos respectively.

**Compile for Intel GPUs using SYCL:**

NOTE: To build with Kokkos support for Intel GPUs via SYCL, Intel's
oneAPI toolkit (version 2024.0 or later) must be installed on your system.
Use the icpx compiler from the Intel oneAPI toolkit. The GPU architecture
is detected automatically at runtime; no arch flag override is needed.

The kokkos\_sycl.cmake preset targets Intel Ponte Vecchio (PVC) GPUs and
uses Intel MKL for FFTs:


.. parsed-literal::

   mkdir build
   cd build
   cmake -C /path/to/sparta/cmake/presets/kokkos_sycl.cmake /path/to/sparta/cmake
   make -j 4

The resulting executable will be named spa\_kokkos\_sycl.

**Running SPARTA with the KOKKOS package:**

All Kokkos operations occur within the context of an individual MPI
task running on a single node of the machine. The total number of MPI
tasks used by SPARTA (one or multiple per compute node) is set in the
usual manner via the mpirun or mpiexec commands, and is independent of
Kokkos. The mpirun or mpiexec command sets the total number of MPI
tasks used by SPARTA (one or multiple per compute node) and the number
of MPI tasks used per node. E.g. the mpirun command in OpenMPI does
this via its -np and -npernode switches. Ditto for MPICH via -np and
-ppn.

**Running on a multi-core CPU:**

Here is a quick overview of how to use the KOKKOS package for CPU
acceleration, assuming one or more 16-core nodes.


.. parsed-literal::

   mpirun -np 16 spa_kokkos_mpi_only -k on -sf kk -in in.collide        # 1 node, 16 MPI tasks/node, no multi-threading
   mpirun -np 2 -ppn 1 spa_kokkos_omp -k on t 16 -sf kk -in in.collide  # 2 nodes, 1 MPI task/node, 16 threads/task
   mpirun -np 2 spa_kokkos_omp -k on t 8 -sf kk -in in.collide          # 1 node,  2 MPI tasks/node, 8 threads/task
   mpirun -np 32 -ppn 4 spa_kokkos_omp -k on t 4 -sf kk -in in.collide  # 8 nodes, 4 MPI tasks/node, 4 threads/task

To run using the KOKKOS package, use the "-k on", "-sf kk" and "-pk
kokkos" `command-line switches <Section_start.html#start_7>`_ in your
mpirun command.  You must use the "-k on" `command-line switch <Section_start.html#start_7>`_ to enable the KOKKOS package. It
takes additional arguments for hardware settings appropriate to your
system. Those arguments are `documented here <Section_start.html#start_7>`_. For OpenMP use:


.. parsed-literal::

   -k on t Nt

The "t Nt" option specifies how many OpenMP threads per MPI task to
use with a node. The default is Nt = 1, which is MPI-only mode.  Note
that the product of MPI tasks \* OpenMP threads/task should not exceed
the physical number of cores (on a node), otherwise performance will
suffer. If hyperthreading is enabled, then the product of MPI tasks \*
OpenMP threads/task should not exceed the physical number of cores \*
hardware threads.  The "-k on" switch also issues a "package kokkos"
command (with no additional arguments) which sets various KOKKOS
options to default values, as discussed on the :doc:`package <package>`
command doc page.

The "-sf kk" `command-line switch <Section_start.html#start_7>`_ will
automatically append the "/kk" suffix to styles that support it.  In
this manner no modification to the input script is
needed. Alternatively, one can run with the KOKKOS package by editing
the input script as described below.

NOTE: When using a single OpenMP thread, the Kokkos Serial backend will give better performance than the OpenMP 
backend because some of the overhead to make 
the code thread-safe is removed.

NOTE: The default for the :doc:`package kokkos <package>` command is to
use "threaded" communication. However, when running on CPUs, it will
typically be faster to use "classic" non-threaded communication.  Use
the "-pk kokkos" `command-line switch <Section_start.html#start_7>`_ to
change the default :doc:`package kokkos <package>` options. See its doc
page for details and default settings. Experimenting with its options
can provide a speed-up for specific calculations. For example:


.. parsed-literal::

   mpirun -np 16 spa_kokkos_mpi_only -k on -sf kk -pk kokkos comm classic -in in.collide       # non-threaded comm

For OpenMP, the KOKKOS package uses data duplication (i.e. 
thread-private arrays) by default to avoid thread-level write conflicts 
in some compute styles. Data duplication is typically fastest for small 
numbers of threads (i.e. 8 or less) but does increase memory footprint 
and is not scalable to large numbers of threads. An alternative to data 
duplication is to use thread-level atomics, which don't require 
duplication. When using the Kokkos Serial backend or the OpenMP backend 
with a single thread, no duplication or atomics are used. For CUDA, the 
KOKKOS package always uses atomics in these computes when necessary. The 
use of atomics instead of duplication can be forced by compiling with the 
"-DSPARTA\_KOKKOS\_USE\_ATOMICS" compile switch.

**Core and Thread Affinity:**

When using multi-threading, it is important for performance to bind
both MPI tasks to physical cores, and threads to physical cores, so
they do not migrate during a simulation.

If you are not certain MPI tasks are being bound (check the defaults
for your MPI installation), binding can be forced with these flags:


.. parsed-literal::

   OpenMPI 1.8: mpirun -np 2 -bind-to socket -map-by socket ./spa_openmpi ...
   Mvapich2 2.0: mpiexec -np 2 -bind-to socket -map-by socket ./spa_mvapich ...

For binding threads with KOKKOS OpenMP, use thread affinity
environment variables to force binding. With OpenMP 3.1 (gcc 4.7 or
later, intel 12 or later) setting the environment variable
OMP\_PROC\_BIND=true should be sufficient. In general, for best
performance with OpenMP 4.0 or better set OMP\_PROC\_BIND=spread and
OMP\_PLACES=threads.  For binding threads with the KOKKOS pthreads
option, compile SPARTA the KOKKOS HWLOC=yes option as described below.

**Running on GPUs:**

Use the "-k" `command-line switch <Section_start.html#start_7>`_ to
specify the number of GPUs per node, and the number of threads per MPI
task. Typically the -np setting of the mpirun command should set the
number of MPI tasks/node to be equal to the # of physical GPUs on the
node.  You can assign multiple MPI tasks to the same GPU with the
KOKKOS package, but this is usually only faster if significant
portions of the input script have not been ported to use Kokkos. Using
CUDA MPS is recommended in this scenario. As above for multi-core CPUs
(and no GPU), if N is the number of physical cores/node, then the
number of MPI tasks/node should not exceed N.


.. parsed-literal::

   -k on g Ng

Here are examples of how to use the KOKKOS package for GPUs, assuming
one or more nodes, each with two GPUs.

**NVIDIA GPUs (CUDA):**


.. parsed-literal::

   mpirun -np 2 spa_kokkos_cuda -k on g 2 -sf kk -in in.collide          # 1 node,   2 MPI tasks/node, 2 GPUs/node
   mpirun -np 32 -ppn 2 spa_kokkos_cuda -k on g 2 -sf kk -in in.collide  # 16 nodes, 2 MPI tasks/node, 2 GPUs/node (32 GPUs total)

NOTE: Use the "-pk kokkos" `command-line switch <Section_start.html#start_7>`_ to change the default :doc:`package kokkos <package>` options. See its doc page for details and default
settings. For example:


.. parsed-literal::

   mpirun -np 2 spa_kokkos_cuda -k on g 2 -sf kk -pk kokkos gpu/aware off -in in.collide      # set gpu/aware MPI support off

**AMD GPUs (HIP):**


.. parsed-literal::

   mpirun -np 2 spa_kokkos_hip -k on g 2 -sf kk -in in.collide          # 1 node,   2 MPI tasks/node, 2 GPUs/node
   mpirun -np 16 -ppn 8 spa_kokkos_hip -k on g 8 -sf kk -in in.collide  # 2 nodes,  8 MPI tasks/node, 8 GPUs/node

For Cray MPICH systems (e.g. El Capitan), set the GTL (GPU Transfer
Library) environment variable before running:


.. parsed-literal::

   export MPICH_GPU_SUPPORT_ENABLED=1
   mpirun -np 8 spa_elcapitan_kokkos -k on g 8 -sf kk -in in.collide

**Intel GPUs (SYCL):**


.. parsed-literal::

   mpirun -np 4 spa_kokkos_sycl -k on g 4 -sf kk -in in.collide          # 1 node, 4 MPI tasks/node, 4 GPUs/node

NOTE: Using OpenMP threading and CUDA/HIP/SYCL together is currently not
possible with the SPARTA KOKKOS package.

NOTE: For good performance of the KOKKOS package on GPUs, you must
have Kepler generation GPUs (or later). The Kokkos library exploits
texture cache options not supported by Tesla generation GPUs (or
older).

NOTE: When using a GPU, you will achieve the best performance if your
input script does not use fix or compute styles which are not yet
Kokkos-enabled. This allows data to stay on the GPU for multiple
timesteps, without being copied back to the host CPU. Invoking a
non-Kokkos fix or compute, or performing I/O for :doc:`stat <stats>` or
:doc:`dump <dump>` output will cause data to be copied back to the CPU
incurring a performance penalty.

NOTE: Most non-Kokkos styles degrade this way, costing performance but
still producing correct results.  A few, however, are rejected outright
and will stop the run, rather than falling back to the host.  Those
restrictions are documented on the page for the style that imposes them.

NOTE: A run may define any number of instances of the styles that the
device kernels capture -- :doc:`surf_collide <surf_collide>`,
:doc:`surf_react <surf_react>`, and the per-event tally computes such as
:doc:`compute boundary <compute_boundary>`, :doc:`compute surf <compute_surf>`,
:doc:`compute isurf/grid <compute_isurf_grid>`,
:doc:`compute react/boundary <compute_react_boundary>`,
:doc:`compute react/surf <compute_react_surf>` and
:doc:`compute react/isurf/grid <compute_react_isurf_grid>`.  Each is held in a
runtime-sized device buffer, so there is no compile-time cap.  The one
remaining exception is that at most two :doc:`compute surf <compute_surf>`
instances may tally for a single :doc:`fix emit/surf <fix_emit_surf>`;
exceeding that stops the run with an explanatory message.

NOTE: Building with -DSPARTA_KOKKOS_FIXED_LISTS instead holds those
styles in fixed-size arrays inside the kernel functor, which restores
the former caps: at most two instances of each
:doc:`surf_react <surf_react>` style, at most four surf react models in
total, at most two active instances of each per-event tally compute
listed above, and at most four gas-phase tally computes.  The two
layouts produce identical results; which is faster is a
hardware-dependent question, so both are kept so they can be compared
by rebuilding.

**Run with the KOKKOS package by editing an input script:**

Alternatively the effect of the "-sf" or "-pk" switches can be
duplicated by adding the :doc:`package kokkos <package>` or :doc:`suffix kk <suffix>` commands to your input script.

The discussion above for building SPARTA with the KOKKOS package, the
mpirun/mpiexec command, and setting appropriate thread are the same.

You must still use the "-k on" `command-line switch <Section_start.html#start_7>`_ to enable the KOKKOS package, and
specify its additional arguments for hardware options appropriate to
your system, as documented above.

You can use the :doc:`suffix kk <suffix>` command, or you can explicitly add a
"kk" suffix to individual styles in your input script, e.g.


.. parsed-literal::

   collide vss/kk air ar.vss

You only need to use the :doc:`package kokkos <package>` command if you
wish to change any of its option defaults, as set by the "-k on"
`command-line switch <Section_start.html#start_7>`_.

**Speed-ups to expect:**

The performance of KOKKOS running in different modes is a function of
your hardware, which KOKKOS-enable styles are used, and the problem
size.

Generally speaking, when running on CPUs only, with a single thread per MPI task, the
performance difference of a KOKKOS style and (un-accelerated) styles
(MPI-only mode) is typically small (less than 20%).

See the `Benchmark page <https://sparta.github.io/bench.html>`_ of the
SPARTA web site for performance of the KOKKOS package on different
hardware.

**Advanced Kokkos options:**

There are other allowed options when building with the KOKKOS package.
A few options are listed here; for a full list of all options, please
refer to the Kokkos documentation.  As above, these options can be set
as variables on the command line or in a CMake presets file. For
default CMake values, see cmake -LH \| grep -i kokkos.

The CMake option Kokkos\_ENABLE\_\ *OPTION* sets the 
parallelization method used for Kokkos code (within SPARTA). 
For example, the CMake option Kokkos\_ENABLE\_SERIAL=ON 
means that no threading will be used.  The CMake option Kokkos\_ENABLE\_OPENMP=ON
means that OpenMP threading will be
used. The CMake option Kokkos\_ENABLE\_CUDA=ON
means an NVIDIA GPU running CUDA will be used.

As described above, the CMake option Kokkos\_ARCH\_\ *TYPE*\ =ON enables compiler switches needed when compiling for a specific hardware:

+------------------+-----------------+---------------------------------------------------------+
| **Arch-ID**      | **HOST or GPU** | **Description**                                         |
+------------------+-----------------+---------------------------------------------------------+
| NATIVE           | HOST            | Local machine                                           |
+------------------+-----------------+---------------------------------------------------------+
| AMDAVX           | HOST            | AMD chip                                                |
+------------------+-----------------+---------------------------------------------------------+
| ARMV80           | HOST            | ARMv8.0 Compatible CPU                                  |
+------------------+-----------------+---------------------------------------------------------+
| ARMV81           | HOST            | ARMv8.1 Compatible CPU                                  |
+------------------+-----------------+---------------------------------------------------------+
| ARMV84           | HOST            | ARMv8.4 Compatible CPU                                  |
+------------------+-----------------+---------------------------------------------------------+
| ARMV84\_SVE      | HOST            | Generic ARMv8.4 with SVE support (-march=armv8.4-a+sve) |
+------------------+-----------------+---------------------------------------------------------+
| ARMV8\_THUNDERX  | HOST            | ARMv8 Cavium ThunderX CPU                               |
+------------------+-----------------+---------------------------------------------------------+
| ARMV8\_THUNDERX2 | HOST            | ARMv8 Cavium ThunderX2 CPU                              |
+------------------+-----------------+---------------------------------------------------------+
| A64FX            | HOST            | ARMv8.2 with SVE Support                                |
+------------------+-----------------+---------------------------------------------------------+
| ARMV9\_GRACE     | HOST            | ARMv9 NVIDIA Grace CPU                                  |
+------------------+-----------------+---------------------------------------------------------+
| SNB              | HOST            | Intel Sandy/Ivy Bridge CPUs                             |
+------------------+-----------------+---------------------------------------------------------+
| HSW              | HOST            | Intel Haswell CPUs                                      |
+------------------+-----------------+---------------------------------------------------------+
| BDW              | HOST            | Intel Broadwell Xeon E-class CPUs                       |
+------------------+-----------------+---------------------------------------------------------+
| ICL              | HOST            | Intel Ice Lake Client CPUs (AVX512)                     |
+------------------+-----------------+---------------------------------------------------------+
| ICX              | HOST            | Intel Ice Lake Xeon Server CPUs (AVX512)                |
+------------------+-----------------+---------------------------------------------------------+
| SKL              | HOST            | Intel Skylake Client CPUs                               |
+------------------+-----------------+---------------------------------------------------------+
| SKX              | HOST            | Intel Skylake Xeon Server CPUs (AVX512)                 |
+------------------+-----------------+---------------------------------------------------------+
| KNC              | HOST            | Intel Knights Corner Xeon Phi                           |
+------------------+-----------------+---------------------------------------------------------+
| KNL              | HOST            | Intel Knights Landing Xeon Phi                          |
+------------------+-----------------+---------------------------------------------------------+
| SPR              | HOST            | Intel Sapphire Rapids Xeon Server CPUs (AVX512)         |
+------------------+-----------------+---------------------------------------------------------+
| POWER8           | HOST            | IBM POWER8 CPUs                                         |
+------------------+-----------------+---------------------------------------------------------+
| POWER9           | HOST            | IBM POWER9 CPUs                                         |
+------------------+-----------------+---------------------------------------------------------+
| ZEN              | HOST            | AMD Zen architecture                                    |
+------------------+-----------------+---------------------------------------------------------+
| ZEN2             | HOST            | AMD Zen2 architecture                                   |
+------------------+-----------------+---------------------------------------------------------+
| ZEN3             | HOST            | AMD Zen3 architecture                                   |
+------------------+-----------------+---------------------------------------------------------+
| ZEN4             | HOST            | AMD Zen4 architecture                                   |
+------------------+-----------------+---------------------------------------------------------+
| ZEN5             | HOST            | AMD Zen5 architecture                                   |
+------------------+-----------------+---------------------------------------------------------+
| RISCV\_SG2042    | HOST            | SG2042 (RISC-V) CPUs                                    |
+------------------+-----------------+---------------------------------------------------------+
| RISCV\_RVA22V    | HOST            | RVA22V (RISC-V) CPUs                                    |
+------------------+-----------------+---------------------------------------------------------+
| RISCV\_U74MC     | HOST            | U74MC (RISC-V) CPUs                                     |
+------------------+-----------------+---------------------------------------------------------+
| MAXWELL50        | GPU             | NVIDIA Maxwell generation CC 5.0                        |
+------------------+-----------------+---------------------------------------------------------+
| MAXWELL52        | GPU             | NVIDIA Maxwell generation CC 5.2                        |
+------------------+-----------------+---------------------------------------------------------+
| MAXWELL53        | GPU             | NVIDIA Maxwell generation CC 5.3                        |
+------------------+-----------------+---------------------------------------------------------+
| PASCAL60         | GPU             | NVIDIA Pascal generation CC 6.0                         |
+------------------+-----------------+---------------------------------------------------------+
| PASCAL61         | GPU             | NVIDIA Pascal generation CC 6.1                         |
+------------------+-----------------+---------------------------------------------------------+
| VOLTA70          | GPU             | NVIDIA Volta generation CC 7.0                          |
+------------------+-----------------+---------------------------------------------------------+
| VOLTA72          | GPU             | NVIDIA Volta generation CC 7.2                          |
+------------------+-----------------+---------------------------------------------------------+
| TURING75         | GPU             | NVIDIA Turing generation CC 7.5                         |
+------------------+-----------------+---------------------------------------------------------+
| AMPERE80         | GPU             | NVIDIA Ampere generation CC 8.0                         |
+------------------+-----------------+---------------------------------------------------------+
| AMPERE86         | GPU             | NVIDIA Ampere generation CC 8.6                         |
+------------------+-----------------+---------------------------------------------------------+
| AMPERE87         | GPU             | NVIDIA Ampere generation CC 8.7                         |
+------------------+-----------------+---------------------------------------------------------+
| ADA89            | GPU             | NVIDIA Ada generation CC 8.9                            |
+------------------+-----------------+---------------------------------------------------------+
| HOPPER90         | GPU             | NVIDIA Hopper generation CC 9.0                         |
+------------------+-----------------+---------------------------------------------------------+
| BLACKWELL100     | GPU             | NVIDIA Blackwell generation CC 10.0                     |
+------------------+-----------------+---------------------------------------------------------+
| BLACKWELL103     | GPU             | NVIDIA Blackwell generation CC 10.3                     |
+------------------+-----------------+---------------------------------------------------------+
| BLACKWELL120     | GPU             | NVIDIA Blackwell generation CC 12.0                     |
+------------------+-----------------+---------------------------------------------------------+
| BLACKWELL121     | GPU             | NVIDIA Blackwell generation CC 12.1                     |
+------------------+-----------------+---------------------------------------------------------+
| AMD\_GFX906      | GPU             | AMD GPU MI50/60                                         |
+------------------+-----------------+---------------------------------------------------------+
| AMD\_GFX908      | GPU             | AMD GPU MI100                                           |
+------------------+-----------------+---------------------------------------------------------+
| AMD\_GFX90A      | GPU             | AMD GPU MI200                                           |
+------------------+-----------------+---------------------------------------------------------+
| AMD\_GFX940      | GPU             | AMD GPU MI300                                           |
+------------------+-----------------+---------------------------------------------------------+
| AMD\_GFX942      | GPU             | AMD GPU MI300                                           |
+------------------+-----------------+---------------------------------------------------------+
| AMD\_GFX942\_APU | GPU             | AMD APU MI300A                                          |
+------------------+-----------------+---------------------------------------------------------+
| AMD\_GFX950      | GPU             | AMD GPU MI350                                           |
+------------------+-----------------+---------------------------------------------------------+
| AMD\_GFX1030     | GPU             | AMD GPU V620/W6800                                      |
+------------------+-----------------+---------------------------------------------------------+
| AMD\_GFX1100     | GPU             | AMD GPU RX7900XTX                                       |
+------------------+-----------------+---------------------------------------------------------+
| AMD\_GFX1101     | GPU             | AMD GPU RX7800XT/RX7700XT                               |
+------------------+-----------------+---------------------------------------------------------+
| AMD\_GFX1103     | GPU             | AMD APU Phoenix                                         |
+------------------+-----------------+---------------------------------------------------------+
| AMD\_GFX1151     | GPU             | AMD APU Strix Halo                                      |
+------------------+-----------------+---------------------------------------------------------+
| AMD\_GFX1152     | GPU             | AMD GPU Radeon 860M                                     |
+------------------+-----------------+---------------------------------------------------------+
| AMD\_GFX1201     | GPU             | AMD GPU RX9070XT                                        |
+------------------+-----------------+---------------------------------------------------------+
| INTEL\_GEN       | GPU             | SPIR64-based devices, e.g. Intel GPUs, using JIT        |
+------------------+-----------------+---------------------------------------------------------+
| INTEL\_DG1       | GPU             | Intel Iris XeMAX GPU                                    |
+------------------+-----------------+---------------------------------------------------------+
| INTEL\_GEN9      | GPU             | Intel GPU Gen9                                          |
+------------------+-----------------+---------------------------------------------------------+
| INTEL\_GEN11     | GPU             | Intel GPU Gen11                                         |
+------------------+-----------------+---------------------------------------------------------+
| INTEL\_GEN12LP   | GPU             | Intel GPU Gen12LP                                       |
+------------------+-----------------+---------------------------------------------------------+
| INTEL\_XEHP      | GPU             | Intel GPU Xe-HP                                         |
+------------------+-----------------+---------------------------------------------------------+
| INTEL\_PVC       | GPU             | Intel GPU Ponte Vecchio                                 |
+------------------+-----------------+---------------------------------------------------------+
| INTEL\_DG2       | GPU             | Intel GPU DG2                                           |
+------------------+-----------------+---------------------------------------------------------+
|                  |                 |                                                         |
+------------------+-----------------+---------------------------------------------------------+

This list was last updated for version 5.2.2 of the Kokkos library.

The CMake option Kokkos\_ENABLE\_CUDA\_\ *OPTION* enables additional
options for CUDA. For example, the CMake option
Kokkos\_ENABLE\_IMPL\_CUDA\_UNIFIED\_MEMORY=ON makes Kokkos allocate all GPU
memory as CUDA managed memory, which the host can read and write
directly. This allows a simulation to use memory on the host to
supplement the memory on the GPU (with some performance penalty), so
that larger problems can be run than would otherwise fit on the GPU,
and it allows host code to access GPU data directly, which is useful
when developing or debugging a Kokkos-enabled style within SPARTA. It
requires CUDA 12.2 or later and a GPU with support for concurrent
managed access, which is any NVIDIA GPU since the Pascal generation
running under Linux. Kokkos classifies this as an internal option that
may change in a future release. It replaces the option
Kokkos\_ENABLE\_CUDA\_UVM=ON, which Kokkos no longer supports;
configuring with that option now stops with an error.

The CMake option Kokkos\_ENABLE\_DEBUG=ON is useful
when developing a Kokkos-enabled style within SPARTA. This option enables printing of run-time debugging
information that can be useful and also enables runtime bounds
checking on Kokkos data structures, but may slow down performance.

**Restrictions:**

Currently, there are no precision options with the KOKKOS package. All
compilation and computation is performed in double precision.


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
