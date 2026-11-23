14. Security considerations
===========================

This section first describes what a SPARTA input script may do on the
machine where SPARTA is run.  It then describes which software from
outside of the SPARTA distribution is used when building SPARTA and the
manual, whether anything is downloaded, and how that is checked.  It
also covers SPARTA provided by other projects.  It is meant for users
that have to follow security rules of their institution or computing
facility, or that want to know what exactly they are compiling and
running.

| 14.1 `What a SPARTA input script can do <#sec_1>`_
| 14.2 `What is included in the SPARTA distribution <#sec_2>`_
| 14.3 `What is downloaded when building SPARTA <#sec_3>`_
| 14.4 `SPARTA provided by other projects <#sec_4>`_
| 14.5 `Checking a download <#sec_5>`_
|

General information about the security of SPARTA, including how to
report problems, is in the file `SECURITY.md
<https://github.com/sparta/sparta/blob/master/SECURITY.md>`_ of the
SPARTA repository on GitHub.


----------


.. _sec_1:

.. raw:: html

   <span id="sec_1"></span>

14.1 What a SPARTA input script can do
--------------------------------------

SPARTA is controlled by the commands in its input script.  An input
script should therefore be regarded as a program and not as data:
running an input script that you have received from somebody else is the
same as running a program or a script from that person.  SPARTA has no
means to limit what the commands in an input script may do, so an input
can do everything that you are allowed to do on the machine where SPARTA
is run.  In particular:

* **Files:** An input can read, create, overwrite, and delete any files
  and folders that you have access to, for example with the :doc:`shell
  <shell>` command or with the :doc:`dump <dump>`, :doc:`restart
  <restart>`, and :doc:`write_surf <write_surf>` commands.
* **Other programs:** The :doc:`shell <shell>` command hands its
  arguments to the command-line interpreter of the operating system and
  can thus run any program.  SPARTA also runs other programs on its own:
  ``gzip`` and ``gunzip`` (through ``popen()``) for reading and writing
  compressed files (when compiled with ``-DSPARTA_GZIP``, which is the
  default of the CMake build), and ``ffmpeg`` for creating movies with
  the :doc:`dump movie <dump_image>` command (when compiled with
  ``-DSPARTA_FFMPEG``).
* **Program code:** The :doc:`python <python>` command and
  python-style :doc:`variables <variable>`, which are part of the PYTHON
  package, run Python code.

This has the following consequences:

* Only run input scripts and Python code from sources that you trust,
  and look at input scripts from other people before running them.
* Never run SPARTA as superuser or administrator.  The SPARTA executable
  prints a :doc:`warning <Section_errors>` ("SPARTA is run with superuser or
  administrator privileges") when it is started that way, and CMake
  prints a warning when it is used to configure SPARTA as superuser.
  The warning is printed by the executable only and not when SPARTA is
  used as a library.  The :doc:`shell <shell>` command is always
  available.
* Compression and movie support can be left out by not defining
  ``SPARTA_GZIP`` or ``SPARTA_FFMPEG`` in the compiler flags, and the
  PYTHON package is only included if you ask for it (``-DPKG_PYTHON=ON``
  for CMake).
* If you need to restrict what SPARTA can do beyond that, this has to be
  done with the means of the operating system, for example by running
  SPARTA from a separate account with limited permissions or inside of a
  container.

SPARTA also does not check all of its input data in every detail.  For
performance reasons, some large data files (for example surface files
read with the :doc:`read_surf <read_surf>` command, or restart files)
are only checked partially.  A malformed or malicious file can therefore
crash SPARTA or lead to undefined behavior.


----------


.. _sec_2:

.. raw:: html

   <span id="sec_2"></span>

14.2 What is included in the SPARTA distribution
------------------------------------------------

Changes to SPARTA are proposed as pull requests on GitHub, and the
automated tests defined in ``.github/workflows`` are run on them.

The SPARTA distribution also contains copies of source code from other
projects:

* The Kokkos library is bundled in the ``lib/kokkos`` folder and is used
  by the KOKKOS package.  The copy is updated by the SPARTA developers
  from releases of Kokkos.  Kokkos itself bundles two more libraries in
  ``lib/kokkos/tpls``: mdspan and desul.  The files
  ``mdspan-hash.txt`` and ``desul-hash.txt`` there record the git commit
  of each copy.
* A stripped down copy of the KISS FFT library (version 1.2.9) is part
  of the source code of the FFT package (``src/FFT/kissfft.h``) and of
  the KOKKOS package (``src/KOKKOS/kissfft_kokkos.h``).
* Code from other projects is also part of the core source code in
  ``src``:

  * the 3d Marching Cubes code (``marching_cubes.cpp``), which is
    adapted from the Marching Cubes 33 implementation by Thomas Lewiner,
    and its lookup tables (``lookup_table.h`` and
    ``decrement_lookup_table.h``)
  * the public domain lookup3 hash function by Bob Jenkins
    (``hashlittle.cpp``)
  * a regular expression matcher adapted from the public domain
    tiny-regex-c library (``utils.cpp``)
  * the Jacobi eigenvalue solver from the ``jacobi_pd`` library by
    Andrew Jewett (``math_eigen.cpp`` and ``math_eigen_impl.h``)

Since these copies are part of the SPARTA source code, they are part of
what is downloaded when you obtain SPARTA.  Nothing else is needed to
compile SPARTA in the default configuration besides a C++ compiler and
MPI (or the MPI stubs that are included in the distribution).


----------


.. _sec_3:

.. raw:: html

   <span id="sec_3"></span>

14.3 What is downloaded when building SPARTA
--------------------------------------------

**Building the executable.**  The CMake build files (``cmake``) and the
makefiles (``src/Makefile`` and ``src/MAKE``) of SPARTA do **not**
download any source code or binary files.  The only use of CMake's
``ExternalProject`` is in the "multibuild" mode, where it runs CMake on
the SPARTA source tree itself.  Optional software is searched for on your
machine with ``find_package()`` and must be installed beforehand.  The
following external software can be used.  All of it is developed
independently of SPARTA, and what version is used depends on what is
installed on your machine:

* MPI (required for parallel runs, ``BUILD_MPI``)
* FFTW3 or Intel MKL (optional alternatives to the included KISS FFT
  for the FFT package)
* JPEG and PNG libraries (``BUILD_JPEG``, ``BUILD_PNG``, for the
  :doc:`dump image <dump_image>` command)
* VTK (``PKG_VTK``, ``BUILD_VTK``, for the dump vtk style)
* Python (``PKG_PYTHON``, for the PYTHON package; also used to run the
  tests)
* CUDA, HIP, SYCL, OpenMP, or hwloc, when the KOKKOS package is
  configured to use them

Because nothing is downloaded, there are no checksums to verify for the
build of the executable.  When you install these libraries with your
operating system's package manager or with a module system on a computing
facility, their integrity is covered by the mechanisms of those tools.

**Building the manual.**  The manual is built with ``make -C doc html``.
This downloads software into the ``doc`` folder:

* The Sphinx documentation tools are installed with ``pip`` from the
  Python package index into a virtual environment (``doc/docenv``),
  using the version ranges in ``doc/utils/requirements.txt``.  These
  downloads are **not** checked against fixed checksums.
* The MathJax package is downloaded from the npm registry over HTTPS.
  Its version is fixed in ``doc/Makefile`` together with its SHA-256
  checksum.  The checksum is compared before the archive is unpacked and
  the manual build stops with an error if it does not match.

None of this is needed to compile or run SPARTA.

**Continuous integration.**  The workflows in ``.github/workflows`` run
GitHub Actions and install packages from the operating system.  They
only have read access to the repository, and the version updates of the
actions are proposed by Dependabot.


----------


.. _sec_4:

.. raw:: html

   <span id="sec_4"></span>

14.4 SPARTA provided by other projects
--------------------------------------

The SPARTA developers only distribute SPARTA through the `GitHub
repository <https://github.com/sparta/sparta>`_ and the `SPARTA web site
<https://sparta.github.io>`_.  The SPARTA documentation does not list or
endorse packages of SPARTA from other projects (for example Linux
distributions or package managers such as Spack).  If you use such a
package or any other pre-compiled executable, the SPARTA developers have
no control over how it was built, which options and patches were used,
or how it was checked, and you have to decide whether you trust that
source.  If you need to know exactly what you are running, compile SPARTA
yourself from the source code and review which packages and libraries
are included.


----------


.. _sec_5:

.. raw:: html

   <span id="sec_5"></span>

14.5 Checking a download
------------------------

SPARTA releases are dated versions of the source code with a tag in the
git repository.  Currently there are no signed tags, no lists of
checksums, and no other cryptographic attestations for SPARTA releases.
If you want to make sure that the source code that you have is the same
as what is on GitHub, the most reliable way is to use git over HTTPS:

.. code-block:: bash

   git clone https://github.com/sparta/sparta.git
   cd sparta
   git checkout 27Aug2026
   git status

In a git clone ``git status`` shows if any file differs from the tagged
version.  For a tarball you can unpack it next to such a clone and compare
the two folders, e.g. with ``diff -r``.  If the SPARTA developers start to
sign tags or publish checksums for releases, this section will describe
how to check them.

.. _sws: https://sparta.github.io
