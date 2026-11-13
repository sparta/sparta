# SPARTA documentation build configuration file.

import os

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
}
# what the sidebar shows under the title
html_context = {'display_version': True}
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


def setup(app):
    app.connect('build-finished', copy_full_size_images)
