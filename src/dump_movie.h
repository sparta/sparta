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

#ifdef DUMP_CLASS

DumpStyle(movie,DumpMovie)

#else

#ifndef SPARTA_DUMP_MOVIE_H
#define SPARTA_DUMP_MOVIE_H

#include "dump_image.h"

namespace SPARTA_NS {

class DumpMovie : public DumpImage {
 public:
  DumpMovie(class SPARTA *, int, char**);

  virtual void openfile();
  virtual void init_style();
  virtual int modify_param(int, char **);

 protected:
  double framerate;             // frame rate of animation
  int bitrate;                  // bitrate of video file in kbps
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid dump movie filename

The file produced by dump movie cannot be binary or compressed
and must be a single file for a single processor.

E: Support for writing movies not included

SPARTA was not built with the -DSPARTA_FFMPEG switch in the Makefile

E: Failed to open FFmpeg pipeline to file %s

The specified file cannot be opened.  Check that the path and name are
correct and writable and that the FFmpeg executable can be found and run.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/
