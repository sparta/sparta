"""Pygments lexer for SPARTA input scripts.

Adapted from LAMMPSLexer.py in the LAMMPS documentation, which is also
GPL-2.0 software from Sandia.  Used by ``.. code-block:: SPARTA`` in the
manual.

As in LAMMPS, a command is recognised only as the first word of a line,
and the commands that define something are followed through their
arguments, so their IDs are coloured as identifiers and their style
names as names -- "fix 1 ave/grid ...", "collide vss ...".  The argument
layouts are SPARTA's, not LAMMPS's: a SPARTA fix has no group-ID before
its style, a dump has its style before its group, and so on.

The command names are not listed here.  conf.py reads them from the
source tree -- every CommandStyle() in src/*.h plus the commands
src/input.cpp handles itself -- and passes them to make_lexer(), so a new
command is highlighted as soon as it exists, and a removed one stops
being highlighted, without anyone remembering to edit this file.

Sphinx does not track this file: after changing it, run "make clean" (or
delete doc/doctrees) so every page is highlighted again, not just the ones
whose source changed.
"""

from pygments.lexer import RegexLexer, bygroups, default, words
from pygments.token import (Comment, Keyword, Name, Number, Operator,
                            String, Text, Whitespace)

# How each defining command's leading arguments are laid out, from the
# Syntax line of its page:
#
#   fix ID style args                     compute ID style args
#   dump ID style select-ID N file args   region ID style args
#   surf_collide ID style args            surf_react ID style args
#   variable name style args              group ID which style args
#   collide style args                    react style args
#   balance_grid style args               mixture ID species ...
#   unfix/uncompute/undump/dump_modify ID ...
#   label name   jump file label   next variables ...
LAYOUTS = {
    'fix': 'id_style', 'compute': 'id_style', 'dump': 'id_style',
    'region': 'id_style', 'surf_collide': 'id_style',
    'surf_react': 'id_style', 'variable': 'id_style',
    'group': 'group',
    'collide': 'style', 'react': 'style', 'balance_grid': 'style',
    'mixture': 'id', 'unfix': 'id', 'uncompute': 'id', 'undump': 'id',
    'dump_modify': 'id', 'next': 'id',
    'label': 'label', 'jump': 'jump',
}

# special argument words several commands accept
CONSTANTS = ('NULL', 'EDGE', 'INF', 'SELF')

# one identifier or style word
WORD = r'[\w\-\.\[\]/]+'


def _command(name, state):
    """A command at the start of a line, then the state for its arguments."""
    return (r'^([ \t]*)(' + name + r')([ \t]+)',
            bygroups(Whitespace, Keyword, Whitespace), state)


def make_lexer(commands):
    """Return a lexer class that highlights the given command names."""
    missing = set(LAYOUTS) - set(commands)
    if missing:
        raise ValueError('LAYOUTS names commands the source does not '
                         'define: ' + ', '.join(sorted(missing)))
    plain = sorted(set(commands) - set(LAYOUTS), key=len, reverse=True)

    class SPARTALexer(RegexLexer):
        name = 'SPARTA'
        aliases = ['sparta']

        tokens = {
            'root': [
                (r'#.*?$', Comment.Single),
            ] + [_command(cmd, LAYOUTS[cmd]) for cmd in
                 sorted(LAYOUTS, key=len, reverse=True)] + [
                (r'^([ \t]*)(' + '|'.join(plain) + r')\b',
                 bygroups(Whitespace, Keyword)),
                (words(('if', 'then', 'elif', 'else'), prefix=r'\b',
                       suffix=r'\b'), Keyword),
                (words(CONSTANTS, prefix=r'\b', suffix=r'\b'),
                 Keyword.Constant),
                # a continuation "&" at the end of a line
                (r'&[ \t]*$', String.Char),
                # python command bodies may span lines
                (r'"""', String, 'triple_string'),
                (r'"', String, 'string'),
                (r"'", String, 'single_string'),
                # references: v_name, c_ID, f_ID, and s_, g_, p_ custom
                # attributes; a trailing [i] or [*] is part of it, an
                # operator after it is not
                (r'\b[vcfsgp]_[\w.-]+(\[[\w*]*\])?', Name.Variable),
                (r'\$\(', Name.Variable, 'expression'),
                (r'\$\{', Name.Variable, 'variable'),
                (r'\$\w', Name.Variable),
                (r'[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?\b', Number),
                (r'[\w.\[\]\-/]+', Text),
                # a newline on its own, so the next line's indentation is
                # left for the line-start command rules
                (r'\n', Whitespace),
                (r'[ \t]+', Whitespace),
                (r'[=<>!&|^~+\-*/%]+', Operator),
                (r'.', Text),
            ],
            # ID, then style:  fix ID style, region ID style, ...
            'id_style': [
                (WORD, Name.Variable.Identifier),
                (r'[ \t]+', Whitespace, ('#pop', 'style')),
                default('#pop'),
            ],
            # group ID which style
            'group': [
                (WORD, Name.Variable.Identifier),
                (r'[ \t]+', Whitespace, ('#pop', 'style', 'skip_word')),
                default('#pop'),
            ],
            'skip_word': [
                (WORD, Text),
                (r'[ \t]+', Whitespace, '#pop'),
                default('#pop'),
            ],
            # one style name
            'style': [
                (r'[A-Za-z][\w/\-]*', Name.Class, '#pop'),
                default('#pop'),
            ],
            # one ID
            'id': [
                (WORD, Name.Variable.Identifier, '#pop'),
                default('#pop'),
            ],
            # label name
            'label': [
                (WORD, String.Char, '#pop'),
                default('#pop'),
            ],
            # jump file label
            'jump': [
                (r'[^ \t\n]+', Text),
                (r'[ \t]+', Whitespace, ('#pop', 'label')),
                default('#pop'),
            ],
            'string': [
                (r'\$\{[^}]*\}', String.Interpol),
                (r'[^"$]+', String),
                (r'\$', String),
                (r'"', String, '#pop'),
            ],
            'single_string': [
                (r"[^']+", String),
                (r"'", String, '#pop'),
            ],
            'triple_string': [
                (r'"""', String, '#pop'),
                (r'[^"]+', String),
                (r'"', String),
            ],
            'variable': [
                (r'[^}]+', Name.Variable),
                (r'\}', Name.Variable, '#pop'),
            ],
            'expression': [
                (r'[^()]+', Name.Variable),
                (r'\(', Name.Variable, '#push'),
                (r'\)', Name.Variable, '#pop'),
            ],
        }

    return SPARTALexer
