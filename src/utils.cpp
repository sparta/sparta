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
   File copied and adapted from LAMMPS (lammps.org), October 2024
   Author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "utils.h"

#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "universe.h"
#include "update.h"
#include "variable.h"

#include <cctype>
#include <cerrno>
#include <cstring>
#include <ctime>
#include <stdexcept>

/*
 * Mini regex-module adapted from https://github.com/kokke/tiny-regex-c
 * which is in the public domain.
 *
 * Supports:
 * ---------
 *   '.'        Dot, matches any character
 *   '^'        Start anchor, matches beginning of string
 *   '$'        End anchor, matches end of string
 *   '*'        Asterisk, match zero or more (greedy)
 *   '+'        Plus, match one or more (greedy)
 *   '?'        Question, match zero or one (non-greedy)
 *   '[abc]'    Character class, match if one of {'a', 'b', 'c'}
 *   '[a-zA-Z]' Character ranges, the character set of the ranges { a-z | A-Z }
 *   '\s'       Whitespace, \t \f \r \n \v and spaces
 *   '\S'       Non-whitespace
 *   '\w'       Alphanumeric, [a-zA-Z0-9_]
 *   '\W'       Non-alphanumeric
 *   '\d'       Digits, [0-9]
 *   '\D'       Non-digits
 *   '\i'       Integer chars, [0-9], '+' and '-'
 *   '\I'       Non-integers
 *   '\f'       Floating point number chars, [0-9], '.', 'e', 'E', '+' and '-'
 *   '\F'       Non-floats
 *
 * *NOT* supported:
 *   '[^abc]'   Inverted class
 *   'a|b'      Branches
 *   '(abc)+'   Groups
 */

extern "C" {
/** Match text against a (simplified) regular expression
   * (regexp will be compiled automatically). */
static int re_match(const char *text, const char *pattern);

/** Match find substring that matches a (simplified) regular expression
   * (regexp will be compiled automatically). */
static int re_find(const char *text, const char *pattern, int *matchlen);
}

////////////////////////////////////////////////////////////////////////

using namespace SPARTA_NS;

/** More flexible and specific matching of a string against a pattern.
 *  This function is supposed to be a more safe, more specific and
 *  simple to use API to find pattern matches. The purpose is to replace
 *  uses of either strncmp() or strstr() in the code base to find
 *  sub-strings safely. With strncmp() finding prefixes, the number of
 *  characters to match must be counted, which can lead to errors,
 *  while using "^pattern" will do the same with less problems.
 *  Matching for suffixes using strstr() is not as specific as 'pattern$',
 *  and complex matches, e.g. "^rigid.*\/small.*", to match all small
 *  body optimized rigid fixes require only one test.
 *
 *  The use of std::string arguments allows for simple concatenation
 *  even with char * type variables.
 *  Example: utils::strmatch(text, std::string("^") + charptr)
 */
bool utils::strmatch(const std::string &text, const std::string &pattern)
{
  const int pos = re_match(text.c_str(), pattern.c_str());
  return (pos >= 0);
}

/** This function is a companion function to utils::strmatch(). Arguments
 *  and logic is the same, but instead of a boolean, it returns the
 *  sub-string that matches the regex pattern.  There can be only one match.
 *  This can be used as a more flexible alternative to strstr().
 */
std::string utils::strfind(const std::string &text, const std::string &pattern)
{
  int matchlen;
  const int pos = re_find(text.c_str(), pattern.c_str(), &matchlen);
  if ((pos >= 0) && (matchlen > 0))
    return text.substr(pos, matchlen);
  else
    return "";
}

void utils::missing_cmd_args(const std::string &file, int line, const std::string &cmd,
                             Error *error)
{
  char msg[128];
  sprintf(msg,"Illegal %s command: missing argument(s)",cmd.c_str());
  if (error) error->all(file.c_str(), line, msg);
}

/* specialization for the case of just a single string argument */

void utils::logmesg(SPARTA *sparta, const std::string &mesg)
{
  if (sparta->screen) fputs(mesg.c_str(), sparta->screen);
  if (sparta->logfile) fputs(mesg.c_str(), sparta->logfile);
}

/* ----------------------------------------------------------------------
   read a boolean value from a string
   transform to lower case before checking
   generate an error if is not a legitimate boolean
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

int utils::logical(const char *file, int line, const std::string &str, bool do_abort, SPARTA *sparta)
{
  if (str.empty()) {
    const char msg[] = "Expected boolean parameter instead of NULL or empty string "
                       "in input script or data file";
    if (do_abort)
      sparta->error->one(file, line, msg);
    else
      sparta->error->all(file, line, msg);
  }

  // convert to ascii
  std::string buf(str);
  if (has_utf8(buf)) buf = utf8_subst(buf);

  int rv = 0;
  if ((buf == "yes") || (buf == "on") || (buf == "true") || (buf == "1")) {
    rv = 1;
  } else if ((buf == "no") || (buf == "off") || (buf == "false") || (buf == "0")) {
    rv = 0;
  } else {
    std::string msg("Expected boolean parameter instead of '");
    msg += buf + "' in input script or data file";
    if (do_abort)
      sparta->error->one(file, line, msg.c_str());
    else
      sparta->error->all(file, line, msg.c_str());
  }
  return rv;
}

/* ----------------------------------------------------------------------
   wrapper for logical() that accepts a char pointer instead of a string
------------------------------------------------------------------------- */

int utils::logical(const char *file, int line, const char *str, bool do_abort, SPARTA *sparta)
{
  if (str)
    return logical(file, line, std::string(str), do_abort, sparta);
  else
    return logical(file, line, std::string(""), do_abort, sparta);
}

/* ----------------------------------------------------------------------
   read a floating point value from a string
   generate an error if not a legitimate floating point value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

double utils::numeric(const char *file, int line, const std::string &str, bool do_abort,
                      SPARTA *sparta)
{
  if (str.empty()) {
    const char msg[] = "Expected floating point parameter instead of"
                       " NULL or empty string in input script or data file";
    if (do_abort)
      sparta->error->one(file, line, msg);
    else
      sparta->error->all(file, line, msg);
  }

  std::string buf(str);
  if (has_utf8(buf)) buf = utf8_subst(buf);

  if (!is_double(buf)) {
    std::string msg("Expected floating point parameter instead of '");
    msg += buf + "' in input script or data file";
    if (do_abort)
      sparta->error->one(file, line, msg.c_str());
    else
      sparta->error->all(file, line, msg.c_str());
  }

  double rv = 0;
  char msg[128];
  sprintf(msg,"Floating point number %s in input script or data file is invalid", buf.c_str());
  try {
    std::size_t endpos;
    rv = std::stod(buf, &endpos);
    if (buf.size() != endpos) {
      if (do_abort)
        sparta->error->one(file, line, msg);
      else
        sparta->error->all(file, line, msg);
    }
  } catch (std::invalid_argument const &) {
    if (do_abort)
      sparta->error->one(file, line, msg);
    else
      sparta->error->all(file, line, msg);
  } catch (std::out_of_range const &) {
    sprintf(msg,"Floating point number %s in input script or data file is out of range", buf.c_str());
    if (do_abort)
      sparta->error->one(file, line, msg);
    else
      sparta->error->all(file, line, msg);
  }
  return rv;
}

/* ----------------------------------------------------------------------
   wrapper for numeric() that accepts a char pointer instead of a string
------------------------------------------------------------------------- */

double utils::numeric(const char *file, int line, const char *str, bool do_abort, SPARTA *sparta)
{
  if (str)
    return numeric(file, line, std::string(str), do_abort, sparta);
  else
    return numeric(file, line, std::string(""), do_abort, sparta);
}

/* ----------------------------------------------------------------------
   read an integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

int utils::inumeric(const char *file, int line, const std::string &str, bool do_abort, SPARTA *sparta)
{
  if (str.empty()) {
    const char msg[] = "Expected integer parameter instead of"
                       " NULL or empty string in input script or data file";
    if (do_abort)
      sparta->error->one(file, line, msg);
    else
      sparta->error->all(file, line, msg);
  }

  std::string buf(str);
  if (has_utf8(buf)) buf = utf8_subst(buf);

  if (!is_integer(buf)) {
    std::string msg("Expected integer parameter instead of '");
    msg += buf + "' in input script or data file";
    if (do_abort)
      sparta->error->one(file, line, msg.c_str());
    else
      sparta->error->all(file, line, msg.c_str());
  }

  int rv = 0;
  char msg[128];
  sprintf(msg,"Integer %s in input script or data file is invalid", buf.c_str());
  try {
    std::size_t endpos;
    rv = std::stoi(buf, &endpos);
    if (buf.size() != endpos) {
      if (do_abort)
        sparta->error->one(file, line, msg);
      else
        sparta->error->all(file, line, msg);
    }
  } catch (std::invalid_argument const &) {
    if (do_abort)
      sparta->error->one(file, line, msg);
    else
      sparta->error->all(file, line, msg);
  } catch (std::out_of_range const &) {
    char msg[128];
    sprintf(msg,"Integer %s in input script or data file is out of range", buf.c_str());
    if (do_abort)
      sparta->error->one(file, line, msg);
    else
      sparta->error->all(file, line, msg);
  }
  return rv;
}

/* ----------------------------------------------------------------------
   wrapper for inumeric() that accepts a char pointer instead of a string
------------------------------------------------------------------------- */

int utils::inumeric(const char *file, int line, const char *str, bool do_abort, SPARTA *sparta)
{
  if (str)
    return inumeric(file, line, std::string(str), do_abort, sparta);
  else
    return inumeric(file, line, std::string(""), do_abort, sparta);
}

/* ----------------------------------------------------------------------
   read a big integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

bigint utils::bnumeric(const char *file, int line, const std::string &str, bool do_abort,
                       SPARTA *sparta)
{
  if (str.empty()) {
    const char msg[] = "Expected integer parameter instead of"
                       " NULL or empty string in input script or data file";
    if (do_abort)
      sparta->error->one(file, line, msg);
    else
      sparta->error->all(file, line, msg);
  }

  std::string buf(str);
  if (has_utf8(buf)) buf = utf8_subst(buf);

  if (!is_integer(buf)) {
    std::string msg("Expected integer parameter instead of '");
    msg += buf + "' in input script or data file";
    if (do_abort)
      sparta->error->one(file, line, msg.c_str());
    else
      sparta->error->all(file, line, msg.c_str());
  }

  long long rv = 0;
  char msg[128];
  sprintf(msg,"Integer %s in input script or data file is invalid", buf.c_str());
  try {
    std::size_t endpos;
    rv = std::stoll(buf, &endpos);
    if (buf.size() != endpos) {
      if (do_abort)
        sparta->error->one(file, line, msg);
      else
        sparta->error->all(file, line, msg);
    }
    if ((rv < (-MAXBIGINT - 1) || (rv > MAXBIGINT))) throw std::out_of_range("bigint");
  } catch (std::invalid_argument const &) {
    if (do_abort)
      sparta->error->one(file, line, msg);
    else
      sparta->error->all(file, line, msg);
  } catch (std::out_of_range const &) {
    char msg[128];
    sprintf(msg,"Integer %s in input script or data file is out of range", buf.c_str());
    if (do_abort)
      sparta->error->one(file, line, msg);
    else
      sparta->error->all(file, line, msg);
  }
  return static_cast<bigint>(rv);
}

/* ----------------------------------------------------------------------
   wrapper for bnumeric() that accepts a char pointer instead of a string
------------------------------------------------------------------------- */

bigint utils::bnumeric(const char *file, int line, const char *str, bool do_abort, SPARTA *sparta)
{
  if (str)
    return bnumeric(file, line, std::string(str), do_abort, sparta);
  else
    return bnumeric(file, line, std::string(""), do_abort, sparta);
}

/* ----------------------------------------------------------------------
   Make copy of string in new storage. Works like the (non-portable)
   C-style strdup() but also accepts a C++ string as argument.
------------------------------------------------------------------------- */

char *utils::strdup(const std::string &text)
{
  auto tmp = new char[text.size() + 1];
  strcpy(tmp, text.c_str());    // NOLINT
  return tmp;
}

/* ----------------------------------------------------------------------
   Replace UTF-8 encoded chars with known ASCII equivalents
------------------------------------------------------------------------- */

std::string utils::utf8_subst(const std::string &line)
{
  const auto *const in = (const unsigned char *) line.c_str();
  const int len = line.size();
  std::string out;

  for (int i = 0; i < len; ++i) {

    // UTF-8 2-byte character
    if ((in[i] & 0xe0U) == 0xc0U) {
      if ((i + 1) < len) {
        // NON-BREAKING SPACE (U+00A0)
        if ((in[i] == 0xc2U) && (in[i + 1] == 0xa0U)) out += ' ', ++i;
        // MODIFIER LETTER PLUS SIGN (U+02D6)
        if ((in[i] == 0xcbU) && (in[i + 1] == 0x96U)) out += '+', ++i;
        // MODIFIER LETTER MINUS SIGN (U+02D7)
        if ((in[i] == 0xcbU) && (in[i + 1] == 0x97U)) out += '-', ++i;
      }
      // UTF-8 3-byte character
    } else if ((in[i] & 0xf0U) == 0xe0U) {
      if ((i + 2) < len) {
        // EN QUAD (U+2000)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x80U)) out += ' ', i += 2;
        // EM QUAD (U+2001)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x81U)) out += ' ', i += 2;
        // EN SPACE (U+2002)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x82U)) out += ' ', i += 2;
        // EM SPACE (U+2003)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x83U)) out += ' ', i += 2;
        // THREE-PER-EM SPACE (U+2004)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x84U)) out += ' ', i += 2;
        // FOUR-PER-EM SPACE (U+2005)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x85U)) out += ' ', i += 2;
        // SIX-PER-EM SPACE (U+2006)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x86U)) out += ' ', i += 2;
        // FIGURE SPACE (U+2007)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x87U)) out += ' ', i += 2;
        // PUNCTUATION SPACE (U+2008)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x88U)) out += ' ', i += 2;
        // THIN SPACE (U+2009)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x89U)) out += ' ', i += 2;
        // HAIR SPACE (U+200A)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x8aU)) out += ' ', i += 2;
        // ZERO WIDTH SPACE (U+200B)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x8bU)) out += ' ', i += 2;
        // LEFT SINGLE QUOTATION MARK (U+2018)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x98U)) out += '\'', i += 2;
        // RIGHT SINGLE QUOTATION MARK (U+2019)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x99U)) out += '\'', i += 2;
        // LEFT DOUBLE QUOTATION MARK (U+201C)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x9cU)) out += '"', i += 2;
        // RIGHT DOUBLE QUOTATION MARK (U+201D)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x9dU)) out += '"', i += 2;
        // NARROW NO-BREAK SPACE (U+202F)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0xafU)) out += ' ', i += 2;
        // WORD JOINER (U+2060)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x81U) && (in[i + 2] == 0xa0U)) out += ' ', i += 2;
        // INVISIBLE SEPARATOR (U+2063)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x81U) && (in[i + 2] == 0xa3U)) out += ' ', i += 2;
        // INVISIBLE PLUS (U+2064)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x81U) && (in[i + 2] == 0xa4U)) out += '+', i += 2;
        // MINUS SIGN (U+2212)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x88U) && (in[i + 2] == 0x92U)) out += '-', i += 2;
        // ZERO WIDTH NO-BREAK SPACE (U+FEFF)
        if ((in[i] == 0xefU) && (in[i + 1] == 0xbbU) && (in[i + 2] == 0xbfU)) out += ' ', i += 2;
      }
      // UTF-8 4-byte character
    } else if ((in[i] & 0xf8U) == 0xf0U) {
      if ((i + 3) < len) { ; }
    } else
      out += in[i];
  }
  return out;
}

/* ----------------------------------------------------------------------
   Return whether string is a valid integer number
------------------------------------------------------------------------- */

bool utils::is_integer(const std::string &str)
{
  if (str.empty()) return false;

  return strmatch(str, "^[+-]?\\d+$");
}

/* ----------------------------------------------------------------------
   Return whether string is a valid floating-point number
------------------------------------------------------------------------- */

bool utils::is_double(const std::string &str)
{
  if (str.empty()) return false;

  return strmatch(str, "^[+-]?\\d+\\.?\\d*$") ||
      strmatch(str, "^[+-]?\\d+\\.?\\d*[eE][+-]?\\d+$") || strmatch(str, "^[+-]?\\d*\\.?\\d+$") ||
      strmatch(str, "^[+-]?\\d*\\.?\\d+[eE][+-]?\\d+$");
}

/* ------------------------------------------------------------------ */

extern "C" {

/* Typedef'd pointer to get abstract datatype. */
typedef struct regex_t *re_t;
typedef struct regex_context_t *re_ctx_t;

/* Compile regex string pattern to a regex_t-array. */
static re_t re_compile(re_ctx_t context, const char *pattern);

/* Find matches of the compiled pattern inside text. */
static int re_matchp(const char *text, re_t pattern, int *matchlen);

/* Definitions: */

#define MAX_REGEXP_OBJECTS 256 /* Max number of regex symbols in expression. */
#define MAX_CHAR_CLASS_LEN 256 /* Max length of character-class buffer in.   */

enum {
  RX_UNUSED,
  RX_DOT,
  RX_BEGIN,
  RX_END,
  RX_QUESTIONMARK,
  RX_STAR,
  RX_PLUS,
  RX_CHAR,
  RX_CHAR_CLASS,
  RX_INV_CHAR_CLASS,
  RX_DIGIT,
  RX_NOT_DIGIT,
  RX_INTEGER,
  RX_NOT_INTEGER,
  RX_FLOAT,
  RX_NOT_FLOAT,
  RX_ALPHA,
  RX_NOT_ALPHA,
  RX_WHITESPACE,
  RX_NOT_WHITESPACE /*, BRANCH */
};

typedef struct regex_t {
  unsigned char type; /* CHAR, STAR, etc.                      */
  union {
    unsigned char ch;   /*      the character itself             */
    unsigned char *ccl; /*  OR  a pointer to characters in class */
  } u;
} regex_t;

typedef struct regex_context_t {
  /* MAX_REGEXP_OBJECTS is the max number of symbols in the expression.
       MAX_CHAR_CLASS_LEN determines the size of buffer for chars in all char-classes in the expression. */
  regex_t re_compiled[MAX_REGEXP_OBJECTS];
  unsigned char ccl_buf[MAX_CHAR_CLASS_LEN];
} regex_context_t;

int re_match(const char *text, const char *pattern)
{
  regex_context_t context;
  int dummy;
  return re_matchp(text, re_compile(&context, pattern), &dummy);
}

int re_find(const char *text, const char *pattern, int *matchlen)
{
  regex_context_t context;
  return re_matchp(text, re_compile(&context, pattern), matchlen);
}

/* Private function declarations: */
static int matchpattern(regex_t *pattern, const char *text, int *matchlen);
static int matchcharclass(char c, const char *str);
static int matchstar(regex_t p, regex_t *pattern, const char *text, int *matchlen);
static int matchplus(regex_t p, regex_t *pattern, const char *text, int *matchlen);
static int matchone(regex_t p, char c);
static int matchdigit(char c);
static int matchint(char c);
static int matchfloat(char c);
static int matchalpha(char c);
static int matchwhitespace(char c);
static int matchmetachar(char c, const char *str);
static int matchrange(char c, const char *str);
static int matchdot(char c);
static int ismetachar(char c);

/* Semi-public functions: */
int re_matchp(const char *text, re_t pattern, int *matchlen)
{
  *matchlen = 0;
  if (pattern != nullptr) {
    if (pattern[0].type == RX_BEGIN) {
      return ((matchpattern(&pattern[1], text, matchlen)) ? 0 : -1);
    } else {
      int idx = -1;

      do {
        idx += 1;

        if (matchpattern(pattern, text, matchlen)) {
          if (text[0] == '\0') return -1;

          return idx;
        }
      } while (*text++ != '\0');
    }
  }
  return -1;
}

re_t re_compile(re_ctx_t context, const char *pattern)
{
  regex_t *const re_compiled = context->re_compiled;
  unsigned char *const ccl_buf = context->ccl_buf;
  int ccl_bufidx = 1;

  char c;    /* current char in pattern   */
  int i = 0; /* index into pattern        */
  int j = 0; /* index into re_compiled    */

  while (pattern[i] != '\0' && (j + 1 < MAX_REGEXP_OBJECTS)) {
    c = pattern[i];

    switch (c) {
        /* Meta-characters: */
      case '^': {
        re_compiled[j].type = RX_BEGIN;
      } break;
      case '$': {
        re_compiled[j].type = RX_END;
      } break;
      case '.': {
        re_compiled[j].type = RX_DOT;
      } break;
      case '*': {
        re_compiled[j].type = RX_STAR;
      } break;
      case '+': {
        re_compiled[j].type = RX_PLUS;
      } break;
      case '?': {
        re_compiled[j].type = RX_QUESTIONMARK;
      } break;

        /* Escaped character-classes (\s \w ...): */
      case '\\': {
        if (pattern[i + 1] != '\0') {
          /* Skip the escape-char '\\' */
          i += 1;
          /* ... and check the next */
          switch (pattern[i]) {
              /* Meta-character: */
            case 'd': {
              re_compiled[j].type = RX_DIGIT;
            } break;
            case 'D': {
              re_compiled[j].type = RX_NOT_DIGIT;
            } break;
            case 'i': {
              re_compiled[j].type = RX_INTEGER;
            } break;
            case 'I': {
              re_compiled[j].type = RX_NOT_INTEGER;
            } break;
            case 'f': {
              re_compiled[j].type = RX_FLOAT;
            } break;
            case 'F': {
              re_compiled[j].type = RX_NOT_FLOAT;
            } break;
            case 'w': {
              re_compiled[j].type = RX_ALPHA;
            } break;
            case 'W': {
              re_compiled[j].type = RX_NOT_ALPHA;
            } break;
            case 's': {
              re_compiled[j].type = RX_WHITESPACE;
            } break;
            case 'S': {
              re_compiled[j].type = RX_NOT_WHITESPACE;
            } break;

              /* Escaped character, e.g. '.' or '$' */
            default: {
              re_compiled[j].type = RX_CHAR;
              re_compiled[j].u.ch = pattern[i];
            } break;
          }
        }
        /* '\\' as last char in pattern -> invalid regular expression. */
      } break;

        /* Character class: */
      case '[': {
        /* Remember where the char-buffer starts. */
        int buf_begin = ccl_bufidx;

        /* Look-ahead to determine if negated */
        if (pattern[i + 1] == '^') {
          re_compiled[j].type = RX_INV_CHAR_CLASS;
          i += 1;                  /* Increment i to avoid including '^' in the char-buffer */
          if (pattern[i + 1] == 0) /* incomplete pattern, missing non-zero char after '^' */
          {
            return nullptr;
          }
        } else {
          re_compiled[j].type = RX_CHAR_CLASS;
        }

        /* Copy characters inside [..] to buffer */
        while ((pattern[++i] != ']') && (pattern[i] != '\0')) {
          /* Missing ] */
          if (pattern[i] == '\\') {
            if (ccl_bufidx >= MAX_CHAR_CLASS_LEN - 1) { return nullptr; }
            if (pattern[i + 1] == 0) /* incomplete pattern, missing non-zero char after '\\' */
            {
              return nullptr;
            }
            ccl_buf[ccl_bufidx++] = pattern[i++];
          } else if (ccl_bufidx >= MAX_CHAR_CLASS_LEN) {
            return nullptr;
          }
          ccl_buf[ccl_bufidx++] = pattern[i];
        }
        if (ccl_bufidx >= MAX_CHAR_CLASS_LEN) {
          /* Catches cases such as [00000000000000000000000000000000000000][ */
          return nullptr;
        }
        /* Null-terminate string end */
        ccl_buf[ccl_bufidx++] = 0;
        re_compiled[j].u.ccl = &ccl_buf[buf_begin];
      } break;

        /* Other characters: */
      default: {
        re_compiled[j].type = RX_CHAR;
        re_compiled[j].u.ch = c;
      } break;
    }
    /* no buffer-out-of-bounds access on invalid patterns -
     * see https://github.com/kokke/tiny-regex-c/commit/1a279e04014b70b0695fba559a7c05d55e6ee90b */
    if (pattern[i] == 0) { return nullptr; }

    i += 1;
    j += 1;
  }
  /* 'RX_UNUSED' is a sentinel used to indicate end-of-pattern */
  re_compiled[j].type = RX_UNUSED;

  return (re_t) re_compiled;
}

/* Private functions: */
static int matchdigit(char c)
{
  return isdigit(c);
}

static int matchint(char c)
{
  return (matchdigit(c) || (c == '-') || (c == '+'));
}

static int matchfloat(char c)
{
  return (matchint(c) || (c == '.') || (c == 'e') || (c == 'E'));
}

static int matchalpha(char c)
{
  return isalpha(c);
}

static int matchwhitespace(char c)
{
  return isspace(c);
}

static int matchalphanum(char c)
{
  return ((c == '_') || matchalpha(c) || matchdigit(c));
}

static int matchrange(char c, const char *str)
{
  return ((c != '-') && (str[0] != '\0') && (str[0] != '-') && (str[1] == '-') &&
          (str[1] != '\0') && (str[2] != '\0') && ((c >= str[0]) && (c <= str[2])));
}

static int matchdot(char c)
{
#if defined(RE_DOT_MATCHES_NEWLINE) && (RE_DOT_MATCHES_NEWLINE == 1)
  (void) c;
  return 1;
#else
  return c != '\n' && c != '\r';
#endif
}

static int ismetachar(char c)
{
  return ((c == 's') || (c == 'S') || (c == 'w') || (c == 'W') || (c == 'd') || (c == 'D'));
}

static int matchmetachar(char c, const char *str)
{
  switch (str[0]) {
    case 'd':
      return matchdigit(c);
    case 'D':
      return !matchdigit(c);
    case 'i':
      return matchint(c);
    case 'I':
      return !matchint(c);
    case 'f':
      return matchfloat(c);
    case 'F':
      return !matchfloat(c);
    case 'w':
      return matchalphanum(c);
    case 'W':
      return !matchalphanum(c);
    case 's':
      return matchwhitespace(c);
    case 'S':
      return !matchwhitespace(c);
    default:
      return (c == str[0]);
  }
}

static int matchcharclass(char c, const char *str)
{
  do {
    if (matchrange(c, str)) {
      return 1;
    } else if (str[0] == '\\') {
      /* Escape-char: increment str-ptr and match on next char */
      str += 1;
      if (matchmetachar(c, str)) {
        return 1;
      } else if ((c == str[0]) && !ismetachar(c)) {
        return 1;
      }
    } else if (c == str[0]) {
      if (c == '-') {
        return ((str[-1] == '\0') || (str[1] == '\0'));
      } else {
        return 1;
      }
    }
  } while (*str++ != '\0');

  return 0;
}

static int matchone(regex_t p, char c)
{
  switch (p.type) {
    case RX_DOT:
      return matchdot(c);
    case RX_CHAR_CLASS:
      return matchcharclass(c, (const char *) p.u.ccl);
    case RX_INV_CHAR_CLASS:
      return !matchcharclass(c, (const char *) p.u.ccl);
    case RX_DIGIT:
      return matchdigit(c);
    case RX_NOT_DIGIT:
      return !matchdigit(c);
    case RX_INTEGER:
      return matchint(c);
    case RX_NOT_INTEGER:
      return !matchint(c);
    case RX_FLOAT:
      return matchfloat(c);
    case RX_NOT_FLOAT:
      return !matchfloat(c);
    case RX_ALPHA:
      return matchalphanum(c);
    case RX_NOT_ALPHA:
      return !matchalphanum(c);
    case RX_WHITESPACE:
      return matchwhitespace(c);
    case RX_NOT_WHITESPACE:
      return !matchwhitespace(c);
    default:
      return (p.u.ch == c);
  }
}

static int matchstar(regex_t p, regex_t *pattern, const char *text, int *matchlen)
{
  int prelen = *matchlen;
  const char *prepos = text;
  while ((text[0] != '\0') && matchone(p, *text)) {
    text++;
    (*matchlen)++;
  }
  while (text >= prepos) {
    if (matchpattern(pattern, text--, matchlen)) return 1;
    (*matchlen)--;
  }

  *matchlen = prelen;
  return 0;
}

static int matchplus(regex_t p, regex_t *pattern, const char *text, int *matchlen)
{
  const char *prepos = text;
  while ((text[0] != '\0') && matchone(p, *text)) {
    text++;
    (*matchlen)++;
  }
  while (text > prepos) {
    if (matchpattern(pattern, text--, matchlen)) return 1;
    (*matchlen)--;
  }
  return 0;
}

static int matchquestion(regex_t p, regex_t *pattern, const char *text, int *matchlen)
{
  if (p.type == RX_UNUSED) return 1;
  if (matchpattern(pattern, text, matchlen)) return 1;
  if (*text && matchone(p, *text++)) {
    if (matchpattern(pattern, text, matchlen)) {
      (*matchlen)++;
      return 1;
    }
  }
  return 0;
}

/* Iterative matching */
static int matchpattern(regex_t *pattern, const char *text, int *matchlen)
{
  int pre = *matchlen;
  do {
    if ((pattern[0].type == RX_UNUSED) || (pattern[1].type == RX_QUESTIONMARK)) {
      return matchquestion(pattern[0], &pattern[2], text, matchlen);
    } else if (pattern[1].type == RX_STAR) {
      return matchstar(pattern[0], &pattern[2], text, matchlen);
    } else if (pattern[1].type == RX_PLUS) {
      return matchplus(pattern[0], &pattern[2], text, matchlen);
    } else if ((pattern[0].type == RX_END) && pattern[1].type == RX_UNUSED) {
      return (text[0] == '\0');
    }
    (*matchlen)++;
  } while ((text[0] != '\0') && matchone(*pattern++, *text++));

  *matchlen = pre;
  return 0;
}
}
