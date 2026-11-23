# SPARTA documentation build configuration file.

import glob
import os
import re

from sphinx.util import logging as sphinx_logging

logger = sphinx_logging.getLogger(__name__)

# Derived, not configured.  This file used to be a template with the two
# directories substituted in by the Makefile, which meant conf.py was
# generated, gitignored, had to be cleaned, and carried absolute paths that
# went stale if the tree was ever moved -- silently, in the case of the
# image copy at the bottom.  This file sits in doc/utils/sphinx-config, so
# both paths follow from where it is.
SPARTA_DOC_DIR = os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))))
SPARTA_SOURCE_DIR = os.path.join(os.path.dirname(SPARTA_DOC_DIR), 'src')

# -- Anchor names ----------------------------------------------------------
#
# docutils builds an HTML id by lowercasing a string and replacing every run
# of non-alphanumerics with a hyphen, so the txt2html anchor "start_7" would
# become "start-7" and "Bird94" would become "bird94".  Those anchor names
# are the manual's URLs -- they are linked to from other pages, from the
# SPARTA web site, and from outside -- so renaming them silently breaks
# inbound links.
#
# All 141 of them are kept by the sources themselves: alongside each
# ".. _Bird94:" label the conversion emitted a raw-HTML
# <span id="Bird94"></span>, which docutils passes through untouched.  Grep
# doc/src for "raw:: html" to see them.
#
# Two other mechanisms for the same job used to live here -- a monkeypatch
# of docutils' make_id and fully_normalize_name, and an anchor_compat
# extension driven by a manifest of the published names.  Removing each in
# turn produced a byte-identical build, so both were dead, and a dormant
# second mechanism is worse than none: nobody would notice it had stopped
# working.  doc/utils/parity-check.py is what actually holds the line, by
# comparing every anchor against the txt2html manual on every run.

# -- General ---------------------------------------------------------------

extensions = [
    'sphinx.ext.mathjax',
    'sphinxcontrib.jquery',
]

# Sphinx turns straight quotes into typographic ones by default.  The
# manual quotes input-script syntax constantly, and the txt2html pages used
# straight quotes throughout, so leave the characters as written.
smartquotes = False

# Blocks that do not name a language are not highlighted.  Sphinx's default
# is to colour them as Python, so the hundreds of syntax descriptions and
# file excerpts in "parsed-literal" blocks had words like "for", "with" and
# "or" coloured as Python keywords.  Input scripts are marked
# ".. code-block:: SPARTA" and real Python ".. code-block:: python" (or
# pycon for interactive sessions), and only those are coloured.
highlight_language = 'none'

source_suffix = '.rst'
master_doc = 'Manual'
exclude_patterns = ['_build']

project = 'SPARTA'
copyright = 'Sandia National Laboratories'
author = 'The SPARTA Developers'


def get_sparta_version():
    """Read the release date out of src/version.h.

    SPARTA's "version" is the release date, e.g. "24 Sep 2025".  It is
    taken from the source rather than maintained separately here, and shown
    in the page title and the sidebar; sphinx_rtd_theme 3 dropped the
    footer version line it used to appear in, so setting `version` alone no
    longer puts it anywhere a reader can see.
    """
    try:
        with open(os.path.join(SPARTA_SOURCE_DIR, 'version.h'), 'r') as f:
            line = f.readline()
        return line.split('"')[1].strip()
    except (IOError, IndexError):
        return 'unknown'


version = get_sparta_version()
release = version

# -- HTML ------------------------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'collapse_navigation': False,
    'navigation_depth': 3,
    # the logo stands in for the project name at the top of the sidebar
    'logo_only': True,
}
# what the sidebar shows under the title
html_context = {'display_version': True}

# The SPARTA logo, as the LAMMPS manual has its own: at the top of the
# sidebar on every page, small at the right of the breadcrumbs beside a
# link to the command list (_templates/breadcrumbs.html), and as the
# favicon.
html_logo = '_static/sparta-logo.png'
html_favicon = '_static/sparta.ico'
templates_path = [os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               '_templates')]
html_css_files = ['sparta.css']

html_title = f'SPARTA Documentation ({release})'
html_short_title = 'SPARTA'
html_static_path = ['_static']

# MathJax is served from the built tree rather than a CDN.  The manual is
# read offline -- from a release tarball, or on a machine with no outbound
# network -- and equations that used to be images must not turn into raw
# LaTeX because a CDN is unreachable.  The path is relative, so the built
# tree can be copied anywhere.  The Makefile fetches it into _static.
mathjax_path = 'mathjax/es5/tex-mml-chtml.js'
html_show_sourcelink = False

# Where the manual is published.  sphinx-sitemap writes sitemap.xml beside
# the pages, listing each one under this URL, as the LAMMPS manual does,
# and each page gets a <link rel="canonical"> to its published copy.  The
# site's robots.txt, at the root of sparta.github.io, points to the
# sitemap.  The search and index pages are not content.
extensions.append('sphinx_sitemap')
html_baseurl = 'https://sparta.github.io/doc/'
sitemap_url_scheme = '{link}'
sitemap_excludes = ['search.html', 'genindex.html']

# Equations were pre-rendered images under doc/Eqs in the txt2html manual.
# They are now typeset by MathJax from their LaTeX source, which is the one
# rendering difference this migration introduces; see
# doc/utils/parity-check.py.


# -- PDF -------------------------------------------------------------------
#
# The manual is read on the web; the PDF is the offline copy that ships in
# the release tarball, which is why doc/ builds it at all.
latex_documents = [
    (master_doc, 'SPARTA.tex', 'SPARTA Documentation', author, 'manual'),
]
latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '10pt',
    # The manual quotes input-script syntax constantly, and a long unbroken
    # command name or path runs off the right edge of the page; \sloppy lets
    # LaTeX loosen the line rather than overflow it.  (A \usepackage for
    # seqsplit used to sit here too.  Nothing ever called \seqsplit, so all
    # it did was require texlive-latex-extra, and fail the build outright
    # where only latex-base and latex-recommended are installed.)
    'preamble': r'''
\setcounter{tocdepth}{2}
\sloppy
''',
}
# The SPARTA logo on the title page, as the LAMMPS PDF has its own.
latex_logo = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          '_static', 'sparta-logo.png')
# Sphinx repeats the whole toctree as a "domain index" the manual has no use
# for, and its module index is empty because there are no Python modules.
latex_domain_indices = False


# -- Full-size images ------------------------------------------------------
#
# Several pages show a scaled-down picture that links to the full-size one,
# as "click for larger version".  The thumbnail is a separate file, so the
# full-size image is named only by the link and Sphinx, which copies only
# the images a page actually displays, never copies it.  Left alone every
# one of those links 404s.
#
# The whole JPG directory is copied into the build instead of just the
# linked files: the manual has always published these under JPG/, the web
# site links to them from outside the manual, and keeping the directory
# whole keeps those URLs working.
# Loud rather than silent if the directory is missing, so that a manual
# built without it says why instead of quietly publishing dead links.
def copy_full_size_images(app, exception):
    import shutil

    if exception is not None or app.builder.name not in ('html', 'dirhtml'):
        return
    src = os.path.join(SPARTA_DOC_DIR, 'src', 'JPG')
    if not os.path.isdir(src):
        logger.warning('full-size images not copied: no such directory %s. '
                       'The "click for larger version" links will 404. '
                       'Run "make clean-all" to regenerate conf.py if the '
                       'tree has moved.', src)
        return
    shutil.copytree(src, os.path.join(app.outdir, 'JPG'), dirs_exist_ok=True)


# -- Spell checking ----------------------------------------------------------
#
# "make -C doc spelling" runs the spelling builder from
# sphinxcontrib-spelling, as the LAMMPS manual does.  It needs the enchant
# library and an English dictionary installed on the system, which an HTML
# build does not, so the extension is loaded only when that builder is the
# one asked for.  Every misspelling is a warning, so with -W the target
# fails until each is either fixed or, if it is a real word the dictionary
# does not know, added to false_positives.txt beside this file.
def _requested_builder(argv):
    import sys

    argv = argv or sys.argv
    for i, arg in enumerate(argv):
        if arg in ('-b', '--builder') and i + 1 < len(argv):
            return argv[i + 1]
        if arg.startswith('-b') and len(arg) > 2:
            return arg[2:]
        if arg.startswith('--builder='):
            return arg.split('=', 1)[1]
    return None


if _requested_builder(None) == 'spelling':
    try:
        import enchant  # noqa: F401
        enchant.Dict('en_US')
    except Exception as err:
        raise RuntimeError(
            '"make spelling" needs the enchant library and an en_US '
            'dictionary installed on the system, e.g. the enchant-2 and '
            'hunspell-en-us packages: %s' % err) from err
    extensions.append('sphinxcontrib.spelling')
    spelling_lang = 'en_US'
    spelling_word_list_filename = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'false_positives.txt')
    spelling_warning = True
    spelling_show_suggestions = True
    spelling_ignore_pypi_package_names = False
    spelling_ignore_acronyms = True
    # Section_errors quotes SPARTA's error and warning messages, and their
    # explanations, from the ERROR/WARNING blocks in the src/ headers.  A
    # misspelling there has to be fixed in src/ as well, so the page keeps
    # matching the text a user is shown.


# -- Bibliography -------------------------------------------------------------
#
# Bibliography.rst lists every work the manual cites, once, as LAMMPS's
# manual does.  The pages keep their own references -- their #Bird94-style
# anchors are published URLs -- so the list is generated from them when the
# manual is built rather than kept by hand.  A reference on a page is:
#
#     .. _Bird94:
#
#     .. raw:: html
#
#        <span id="Bird94"></span>
#
#     **(Bird94)** G. A. Bird, Molecular Gas Dynamics ...
#
# Copies that differ only in punctuation are one work.  A key used for two
# different works gives two entries, each linking to the pages that mean
# it -- but should be split into two keys, as Bird09 was into Bird09a and
# Bird09b.
BIB_MARKER = '.. sparta-bibliography-entries'
BIB_LABEL = re.compile(r'^\.\. _([A-Z][A-Za-z-]*\d\d(?:\d\d)?[a-z]?):', re.M)
BIB_ENTRY = re.compile(
    r'^\.\. _(?P<key>[A-Z][A-Za-z-]*\d\d(?:\d\d)?[a-z]?):\n\n'
    r'\.\. raw:: html\n\n   <span id="(?P=key)"></span>\n\n+'
    r'\*\*\((?P=key)\)\*\* (?P<text>.+?)\n\n', re.M | re.S)


def bibliography_entries(src_dir):
    """[(key, text, [page, ...])], one per distinct work, sorted by key."""
    works = {}
    for path in sorted(glob.glob(os.path.join(src_dir, '*.rst'))):
        page = os.path.basename(path)[:-4]
        if page == 'Bibliography':
            continue
        with open(path) as f:
            text = f.read()
        labels = BIB_LABEL.findall(text)
        found = [(m.group('key'), ' '.join(m.group('text').split()))
                 for m in BIB_ENTRY.finditer(text)]
        if sorted(labels) != sorted(k for k, _t in found):
            raise RuntimeError(
                'Bibliography: %s.rst has citation labels %s but readable '
                'references only for %s; each label must be followed by its '
                'raw-html span and a "**(Key)** text" paragraph'
                % (page, sorted(labels), sorted(k for k, _t in found)))
        for key, body in found:
            same = (key, re.sub(r'\W+', '', body).lower())
            work = works.setdefault(same, [key, body, []])
            if len(body) > len(work[1]):
                work[1] = body
            work[2].append(page)
    return sorted(works.values(), key=lambda w: (w[0].lower(), w[1]))


def bibliography_rst(src_dir):
    entries = bibliography_entries(src_dir)
    per_key = {}
    for key, _body, _pages in entries:
        per_key[key] = per_key.get(key, 0) + 1
    out, used = [], {}
    for key, body, pages in entries:
        used[key] = used.get(key, 0) + 1
        anchor = key if per_key[key] == 1 else '%s-%d' % (key, used[key])
        # anonymous links (__): the same page name links to a different
        # anchor in each entry, which docutils < 0.22 rejects for named ones
        cites = ', '.join('`%s <%s.html#%s>`__' % (p.replace('_', '\\_'), p, key)
                          for p in pages)
        out += ['.. raw:: html', '', '   <span id="%s"></span>' % anchor, '',
                '**(%s)** %s' % (key, body), '',
                '   Cited on: %s' % cites, '']
    return '\n'.join(out)


def bibliography_always_outdated(app, env, added, changed, removed):
    # Generated from the other pages, so an incremental build has to
    # re-read it whenever any of them could have changed.
    return ['Bibliography'] if 'Bibliography' in env.found_docs else []


def fill_bibliography(app, docname, source):
    if docname != 'Bibliography':
        return
    if BIB_MARKER not in source[0]:
        raise RuntimeError('Bibliography.rst has lost its "%s" line, so '
                           'the generated entries have nowhere to go'
                           % BIB_MARKER)
    source[0] = source[0].replace(BIB_MARKER, bibliography_rst(app.srcdir))


# -- SPARTA input-script highlighting --------------------------------------
#
# ".. code-block:: SPARTA" uses the lexer in SPARTALexer.py, modelled on the
# one LAMMPS uses.  The command names come from the source tree rather than
# a list kept by hand: every CommandStyle() in src/*.h, plus the commands
# src/input.cpp handles itself.  (The /kk variants in src/KOKKOS are the same
# commands with an accelerator suffix, so src/*.h is enough.)
#
# Any failure here stops the build.  A lexer that fails to load would
# otherwise leave every block unhighlighted and the build would still pass.
def sparta_commands():
    import glob
    import re

    names = set()
    for header in glob.glob(os.path.join(SPARTA_SOURCE_DIR, '*.h')):
        with open(header) as f:
            names.update(re.findall(r'^CommandStyle\((\w+),', f.read(), re.M))
    with open(os.path.join(SPARTA_SOURCE_DIR, 'input.cpp')) as f:
        code = re.sub(r'//[^\n]*', '', f.read())
        names.update(re.findall(r'strcmp\(command,"(\w+)"\)', code))
    if len(names) < 50:
        raise RuntimeError('found only %d SPARTA commands under %s; the '
                           'lexer would miss most of them'
                           % (len(names), SPARTA_SOURCE_DIR))
    return names


# -- Anchors in the PDF ----------------------------------------------------
#
# The manual keeps its published anchor names -- Bird94, start_3, cmd_5 --
# as raw-HTML <span id="..."> beside each reST label (see "Anchor names" at
# the top of this file), because docutils would turn the label itself into
# "bird94" or "start-3".  The links to them say "#Bird94" and
# "Section_start.html#start_3".  LaTeX drops raw HTML, so in the PDF those
# anchors did not exist and every such link went nowhere (155 of them).
#
# For the LaTeX builder only, a target node with the same ids is put in
# front of each raw-HTML node that holds nothing but such spans, and the
# LaTeX writer turns it into the label the links already ask for.  The HTML build is not
# touched: it keeps the raw spans exactly as before.
RAW_ANCHOR = re.compile(r'<span id="([^"]+)"></span>')


def latex_raw_anchors(app, doctree):
    from docutils import nodes

    if app.builder.format != 'latex':
        return
    # ids that become labels; an external link's own target node carries an
    # id too ("stl2surf" for `stl2surf <...>`_) but is never a label
    have = {i for n in doctree.findall(nodes.Element)
            if not (isinstance(n, nodes.target) and 'refuri' in n)
            for i in n.get('ids', [])}
    for raw in list(doctree.findall(nodes.raw)):
        if 'html' not in raw.get('format', '').split():
            continue
        text = raw.astext()
        ids = RAW_ANCHOR.findall(text)
        if not ids or RAW_ANCHOR.sub('', text).strip():
            continue
        new = [i for i in ids if i not in have]
        if not new:
            continue
        have.update(new)
        # inserted before the raw node rather than replacing it: the raw
        # node may carry the ids of the reST label beside it, which the
        # label already gives LaTeX by reference, and the LaTeX writer
        # skips the raw node itself
        #
        # Right before a section heading the LaTeX writer skips targets,
        # expecting the section to carry their ids, so there the ids go on
        # the section.
        following = raw.next_node(descend=False, siblings=True, ascend=True)
        if isinstance(following, nodes.section):
            following['ids'].extend(new)
        else:
            raw.parent.insert(raw.parent.index(raw),
                              nodes.target('', '', ids=new))


def setup(app):
    import sys

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import SPARTALexer

    app.add_lexer('SPARTA', SPARTALexer.make_lexer(sparta_commands()))
    app.connect('source-read', fill_bibliography)
    app.connect('doctree-read', latex_raw_anchors)
    app.connect('env-get-outdated', bibliography_always_outdated)
    app.connect('build-finished', copy_full_size_images)
