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

/* ----------------------------------------------------------------------
   File adapted from LAMMPS (https://www.lammps.org), October 2024
   Ported to SPARTA by: Stan Moore (SNL)
   Original Author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifndef SPARTA_UTILS_H
#define SPARTA_UTILS_H

#include "spatype.h"

#include <mpi.h>

#include <vector>
#include <string>

namespace SPARTA_NS {

// forward declarations
class Error;
class SPARTA;

namespace utils {

  /*! Match text against a simplified regex pattern
   *
   *  \param text the text to be matched against the pattern
   *  \param pattern the search pattern, which may contain regexp markers
   *  \return true if the pattern matches, false if not */

  bool strmatch(const std::string &text, const std::string &pattern);

  /*! Find sub-string that matches a simplified regex pattern
   *
   *  \param text the text to be matched against the pattern
   *  \param pattern the search pattern, which may contain regexp markers
   *  \return the string that matches the pattern or an empty one */

  std::string strfind(const std::string &text, const std::string &pattern);

  /*! Print error message about missing arguments for command
   *
   * This function simplifies the repetitive reporting missing arguments to a command.
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param cmd      name of the failing command
   *  \param error    pointer to Error class instance (for abort) or nullptr */

  void missing_cmd_args(const std::string &file, int line, const std::string &cmd, Error *error);

  /*! \overload
   *
   *  \param sparta    pointer to SPARTA class instance
   *  \param mesg   string with message to be printed */

  void logmesg(SPARTA *sparta, const std::string &mesg);

  int logical(const char *file, int line, const std::string &str, bool do_abort, SPARTA *sparta);

  /*! \overload
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to logical
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param sparta      pointer to top-level SPARTA class instance
   *  \return         1 if string resolves to "true", otherwise 0 */

  int logical(const char *file, int line, const char *str, bool do_abort, SPARTA *sparta);

  /*! Convert a string to a floating point number while checking
   *  if it is a valid floating point or integer number
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to number
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param sparta      pointer to top-level SPARTA class instance
   *  \return         double precision floating point number */

  double numeric(const char *file, int line, const std::string &str, bool do_abort, SPARTA *sparta);

  /*! \overload
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to number
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param sparta      pointer to top-level SPARTA class instance
   *  \return         double precision floating point number */

  double numeric(const char *file, int line, const char *str, bool do_abort, SPARTA *sparta);

  /*! Convert a string to an integer number while checking
   *  if it is a valid integer number (regular int)
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to number
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param sparta      pointer to top-level SPARTA class instance
   *  \return         integer number (regular int)  */

  int inumeric(const char *file, int line, const std::string &str, bool do_abort, SPARTA *sparta);

  /*! \overload
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to number
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param sparta      pointer to top-level SPARTA class instance
   *  \return         double precision floating point number */

  int inumeric(const char *file, int line, const char *str, bool do_abort, SPARTA *sparta);

  /*! Convert a string to an integer number while checking
   *  if it is a valid integer number (bigint)
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to number
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param sparta      pointer to top-level SPARTA class instance
   *  \return         integer number (bigint) */

  bigint bnumeric(const char *file, int line, const std::string &str, bool do_abort, SPARTA *sparta);

  /*! \overload
   *
   *  \param file     name of source file for error message
   *  \param line     line number in source file for error message
   *  \param str      string to be converted to number
   *  \param do_abort determines whether to call Error::one() or Error::all()
   *  \param sparta      pointer to top-level SPARTA class instance
   *  \return         double precision floating point number */

  bigint bnumeric(const char *file, int line, const char *str, bool do_abort, SPARTA *sparta);

  /*! Make C-style copy of string in new storage
   *
   * This allocates a storage buffer and copies the C-style or
   * C++ style string into it.  The buffer is allocated with "new"
   * and thus needs to be deallocated with "delete[]".
   *
   * \param text  string that should be copied
   * \return new buffer with copy of string */

  char *strdup(const std::string &text);

  /*! Check if a string will likely have UTF-8 encoded characters
   *
   * UTF-8 uses the 7-bit standard ASCII table for the first 127 characters and
   * all other characters are encoded as multiple bytes.  For the multi-byte
   * characters the first byte has either the highest two, three, or four bits
   * set followed by a zero bit and followed by one, two, or three more bytes,
   * respectively, where the highest bit is set and the second highest bit set
   * to 0.  The remaining bits combined are the character code, which is thus
   * limited to 21-bits.
   *
   * For the sake of efficiency this test only checks if a character in the string
   * has the highest bit set and thus is very likely an UTF-8 character.  It will
   * not be able to tell this this is a valid UTF-8 character or whether it is a
   * 2-byte, 3-byte, or 4-byte character.
   *
\verbatim embed:rst

*See also*
   :cpp:func:`utils::utf8_subst`

\endverbatim
   * \param line  string that should be checked
   * \return true if string contains UTF-8 encoded characters (bool) */

  inline bool has_utf8(const std::string &line)
  {
    for (auto c : line)
      if (c & 0x80U) return true;
    return false;
  }

  /*! Replace known UTF-8 characters with ASCII equivalents
   *
\verbatim embed:rst

*See also*
   :cpp:func:`utils::has_utf8`

\endverbatim
   * \param line  string that should be converted
   * \return new string with ascii replacements (string) */

  std::string utf8_subst(const std::string &line);

  /*! Check if string can be converted to valid integer
   *
   * \param str string that should be checked
   * \return true, if string contains valid a integer, false otherwise */

  bool is_integer(const std::string &str);

  /*! Check if string can be converted to valid floating-point number
   *
   * \param str string that should be checked
   * \return true, if string contains valid number, false otherwise */

  bool is_double(const std::string &str);

}    // namespace utils
}    // namespace SPARTA_NS

#endif
