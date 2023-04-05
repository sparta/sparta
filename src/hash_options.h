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

// structs for specialized maps/hashes

#ifdef SPARTA_MAP

  struct OnePoint2d {
    double pt[2];

    bool operator <(const OnePoint2d& other) const {
      if (pt[0] < other.pt[0]) return 1;
      else if (pt[0] > other.pt[0]) return 0;
      if (pt[1] < other.pt[1]) return 1;
      else if (pt[1] > other.pt[1]) return 0;
      return 0;
    }
  };

  struct OnePoint3d {
    double pt[3];

    bool operator <(const OnePoint3d& other) const {
      if (pt[0] < other.pt[0]) return 1;
      else if (pt[0] > other.pt[0]) return 0;
      if (pt[1] < other.pt[1]) return 1;
      else if (pt[1] > other.pt[1]) return 0;
      if (pt[2] < other.pt[2]) return 1;
      else if (pt[2] > other.pt[2]) return 0;
      return 0;
    }
  };

  struct TwoPoint3d {
    double pts[6];

    bool operator <(const TwoPoint3d& other) const {
      for (int i = 0; i < 6; i++) {
        if (pts[i] < other.pts[i]) return 1;
        else if (pts[i] > other.pts[i]) return 0;
      }
      return 0;
    }
  };

#else

  struct OnePoint2d {
    double pt[2];

    bool operator ==(const OnePoint2d &other) const {
      if (pt[0] != other.pt[0]) return 0;
      if (pt[1] != other.pt[1]) return 0;
      return 1;
    }
  };

  struct OnePoint2dHash {
    uint32_t operator ()(const OnePoint2d& one) const {
      return hashlittle(one.pt,2*sizeof(double),0);
    }
  };

  struct OnePoint3d {
    double pt[3];

    bool operator ==(const OnePoint3d &other) const {
      if (pt[0] != other.pt[0]) return 0;
      if (pt[1] != other.pt[1]) return 0;
      if (pt[2] != other.pt[2]) return 0;
      return 1;
    }
  };

  struct OnePoint3dHash {
    uint32_t operator ()(const OnePoint3d& one) const {
      return hashlittle(one.pt,3*sizeof(double),0);
    }
  };

  struct TwoPoint3d {
    double pts[6];

    bool operator ==(const TwoPoint3d &other) const {
      for (int i = 0; i < 6; i++)
        if (pts[i] != other.pts[i]) return 0;
      return 1;
    }
  };

  struct TwoPoint3dHash {
    uint32_t operator ()(const TwoPoint3d& two) const {
      return hashlittle(two.pts,6*sizeof(double),0);
    }
  };

#endif
