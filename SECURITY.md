# Security Policy

SPARTA is designed as a user-level application to conduct computer
simulations for research using the Direct Simulation Monte Carlo (DSMC)
method.  As such SPARTA depends to some degrees on users providing
correctly formatted input and SPARTA needs to read and write files based
on uncontrolled user input.  As a parallel application for use in
high-performance computing environments, performance critical steps are
also done without checking data.

SPARTA also uses a number of external libraries (for example MPI, FFTW,
JPEG, PNG, VTK, and Python) and the Kokkos library that is bundled in the
`lib/kokkos` folder.  The external libraries are not validated and tested
by the SPARTA developers, so it is easy to import bad behavior from
calling functions in one of those libraries.  The section [Security
considerations](doc/src/Section_security.rst) of the SPARTA manual (in
the `doc/src` folder) lists the external software that may be used when
compiling SPARTA, explains what is downloaded during the build and how
that is checked, and what applies to SPARTA provided by other projects.

Thus it is quite easy to crash SPARTA through malicious input and do all
kinds of file system manipulations.  A SPARTA input can also run other
programs through the `shell` command, run Python code (PYTHON package,
python-style variables), and run external programs for reading and
writing compressed files (`gzip`, `gunzip`) and for creating movies
(`ffmpeg` with the `dump movie` command).  And because of that SPARTA
should **NEVER** be compiled or **run** as superuser, either from a
"root" or "administrator" account directly or indirectly via "sudo" or
"su".  The SPARTA executable prints a warning when it is started that
way, and so does CMake when configuring SPARTA.

Therefore what could be seen as a security vulnerability is usually
either a user mistake or a bug in the code.  Bugs can be reported in the
SPARTA project [issue tracker on
GitHub](https://github.com/sparta/sparta/issues).

To mitigate issues with using homoglyphs or bidirectional reordering in
unicode, which have been demonstrated as a vector to obfuscate and hide
malicious changes to the source code, all SPARTA pull requests are
checked for unicode characters and only all-ASCII text files are
accepted.  The check can also be run locally with `make check-unicode`
in the `src` folder.  The bundled Kokkos library in `lib/kokkos` is
maintained by the Kokkos project and is excluded from this check.

# Reporting a Vulnerability

If you have found a problem that could be used to harm other SPARTA
users or the SPARTA project, please do **not** report it in a public
issue.  Examples for such problems are ways to get malicious changes
into the SPARTA source code or into the software that is bundled with or
downloaded when compiling SPARTA, and also passwords or access keys that
have become public by accident.  Please report those problems privately
with the "Report a vulnerability" button on the [Security
page](https://github.com/sparta/sparta/security) of the SPARTA
repository on GitHub.

Please include which version of SPARTA is affected and how the problem
can be reproduced.

# Development and Release Model

SPARTA follows a continuous release development model.  Whenever a bug
is fixed or a feature is added, the change is merged to the `master`
branch of the repository at <https://github.com/sparta/sparta>.  A SPARTA
"version" is the date on which it was released (for example
`7Mar2024`), and each dated version contains all features and bug fixes
up to and including that date.  The version date is printed every time
SPARTA is run and is also stored in the file `src/version.h`.  Dated
versions are marked with a tag in the git repository and are announced
on the SPARTA web site.

Changes are proposed as pull requests, which use the pull request
template of the repository, and the continuous integration tests defined
in `.github/workflows/main.yml` are run on pull requests to the `master`
branch.

# Supported Versions

Corrections are applied to the development version (the `master` branch)
and are thus part of the next dated release.  Older releases are not
updated, so please upgrade to a current version.

# Integrity of Downloaded Source Code

The most reliable way to obtain SPARTA is a clone of the git repository
from GitHub over https.  The git object ids of a clone make it possible
to compare what you have with what is on GitHub, e.g. with `git
rev-parse <tag>` and `git status`.  Currently SPARTA releases do **not**
come with signed tags, lists of checksums, or other cryptographic
attestations.  If this changes, this section will be updated.
