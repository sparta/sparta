# Build the SPARTA manual with Sphinx.
#
# Sphinx and its dependencies are installed into a virtualenv under
# doc/docenv the first time a build runs, so the only prerequisite is a
# python3 with venv available.

SHELL          = /bin/bash
BUILDDIR       = $(shell pwd)
RSTDIR         = $(BUILDDIR)/src
VENV           = $(BUILDDIR)/docenv
SPHINXCONFIG   = $(BUILDDIR)/utils/sphinx-config
MATHJAX        = $(SPHINXCONFIG)/_static/mathjax
MATHJAXSTAMP   = $(MATHJAX)/.version-$(MATHJAXVER)
MATHJAXVER     = 3.2.2
MATHJAXURL     = https://registry.npmjs.org/mathjax/-/mathjax-$(MATHJAXVER).tgz
MATHJAXSHA     = 1b9c0a1c44df864e915690558e72adb9cc5203360daefd385084ced3b6c64c09
PYTHON         = python3
SPHINXEXTRA    =

.PHONY: help clean clean-all html pdf fasthtml fasthtml-init check

HAS_PDFLATEX   = $(shell command -v latexmk > /dev/null 2>&1 && \
                         command -v pdflatex > /dev/null 2>&1 && echo YES || echo NO)
HAS_PANDOC     = $(shell command -v pandoc > /dev/null 2>&1 && echo YES || echo NO)
# macOS has shasum but not GNU coreutils' sha256sum, and SPARTA builds on
# macOS (src/MAKE/Makefile.mac); take whichever is there.
SHA256SUM      = $(shell command -v sha256sum > /dev/null 2>&1 \
                         && echo "sha256sum" || echo "shasum -a 256")
FASTHTMLFILES  = $(patsubst $(RSTDIR)/%.rst,fasthtml/%.html,$(wildcard $(RSTDIR)/*.rst))

help:
	@echo "Usage: make <target>"
	@echo ""
	@echo "  html          build the HTML manual in doc/html"
	@echo "  pdf           build doc/html/Manual.pdf (needs a LaTeX installation)"
	@echo "  fasthtml      quick per-page preview in doc/fasthtml (needs pandoc)"
	@echo "  check         verify every internal link and anchor in doc/html"
	@echo "  clean         remove the build tree, keep the virtualenv"
	@echo "  clean-all     also remove the virtualenv and the MathJax copy"
	@echo ""
	@echo "  SPHINXEXTRA=... passes extra flags to sphinx-build, e.g."
	@echo "  make html SPHINXEXTRA=\"-j auto -W --keep-going\""

# The one check that outlives the migration: parity-check.py and
# equation-check.py both need a txt2html manual to compare against, and
# there will not be one once doc/*.txt is gone from the release.
check: html
	@$(PYTHON) $(BUILDDIR)/utils/link-check.py $(BUILDDIR)/html

clean:
	rm -rf $(BUILDDIR)/html $(BUILDDIR)/doctrees \
	       $(BUILDDIR)/latex $(BUILDDIR)/latexdoctrees \
	       $(BUILDDIR)/fasthtml

clean-all: clean
	rm -rf $(VENV) $(MATHJAX)

# the virtualenv is a directory target: it is built once, then only rebuilt
# if requirements.txt changes
#
# Chained with && rather than ; so a failure stops here.  With ; a failed
# "python3 -m venv" would go on to install Sphinx into whatever interpreter
# is on PATH, and a failed "pip install" would still reach the touch below,
# leaving a half-built virtualenv that make then considers up to date.
# A partial virtualenv is also removed on the way out: it would otherwise be
# left newer than requirements.txt, so the next make would take it for a
# finished one and build the manual against half-installed dependencies.
$(VENV): utils/requirements.txt
	@rm -rf $(VENV)
	@{ $(PYTHON) -m venv $(VENV) && \
	   . $(VENV)/bin/activate && \
	   pip install --upgrade pip wheel > /dev/null && \
	   pip install -r $(BUILDDIR)/utils/requirements.txt; \
	 } || { rm -rf $(VENV); exit 1; }
	@touch $(VENV)

# MathJax is fetched rather than committed, so the repository does not
# carry six megabytes of minified JavaScript and fonts.  Version 3 is used
# because it bundles its fonts in the one package; version 4 keeps them in
# a separate 13 MB package, and without it every equation renders with
# broken glyphs.  Pinned and checksummed so the build is reproducible.
# The target is a stamp named for the version, not the directory.  A
# directory target is satisfied by mere existence, so bumping MATHJAXVER
# would have left the old copy in place and never checked the new sum; a
# stamp with the version in its name does not exist yet after a bump, so
# the fetch runs.  It is written last, so an interrupted unpack is not
# mistaken for a finished one either.
$(MATHJAXSTAMP):
	@echo "Fetching MathJax $(MATHJAXVER)"
	@rm -rf $(MATHJAX)
	@mkdir -p $(SPHINXCONFIG)/_static
	@tmp=$$(mktemp -d -p $(SPHINXCONFIG)/_static); \
	 trap 'rm -rf $$tmp' EXIT; \
	 curl -sSL $(MATHJAXURL) -o $$tmp/mathjax.tgz && \
	 echo "$(MATHJAXSHA)  $$tmp/mathjax.tgz" | $(SHA256SUM) -c - > /dev/null && \
	 tar -xzf $$tmp/mathjax.tgz -C $$tmp && \
	 mv $$tmp/package $(MATHJAX) && \
	 echo "$(MATHJAXSHA)" > $(MATHJAXSTAMP)

# && rather than ; here too: with ; the recipe's exit status is the status
# of the last command in the list, so a failed sphinx-build reported success
# and a broken manual looked like a good one.
html: $(VENV) $(MATHJAXSTAMP)
	@. $(VENV)/bin/activate && \
	 sphinx-build $(SPHINXEXTRA) -b html \
	   -c $(SPHINXCONFIG) -d $(BUILDDIR)/doctrees $(RSTDIR) html
	@echo "Build finished. The HTML pages are in doc/html."

# The PDF is the copy that ships in the release tarball.  Sphinx writes the
# LaTeX and its own Makefile into doc/latex; that Makefile drives latexmk,
# which reruns LaTeX until the cross-references and the table of contents
# settle.
#
# It lands in doc/html beside the manual, because that is where the manual
# links to it: Manual.rst says `PDF file <Manual.pdf>`_, which used to
# resolve when the pages and the PDF shared doc/ and would 404 with the
# pages one directory down.
#
# Known limitation.  2266 same-page anchors resolve inside the PDF, but the
# 350 cross-page ones do not: the sources write those as literal URIs
# ("Section_start.html#start_7"), which Sphinx passes through as external
# links, so in the PDF they point at HTML files that are not there.  Making
# them :ref: would fix the PDF and change the anchor spelling in the HTML,
# and those anchor names are published URLs -- so it is a change for the
# step that is allowed to alter them, not for this one.  The HTML manual is
# unaffected; there every one of them resolves.
pdf: $(VENV) $(MATHJAXSTAMP)
	@if [ "$(HAS_PDFLATEX)" = "NO" ]; then \
	   echo "make pdf needs latexmk and pdflatex" 1>&2; exit 1; fi
	@. $(VENV)/bin/activate && \
	 sphinx-build $(SPHINXEXTRA) -b latex \
	   -c $(SPHINXCONFIG) -d $(BUILDDIR)/latexdoctrees $(RSTDIR) latex
	@$(MAKE) -C $(BUILDDIR)/latex LATEXMKOPTS="-interaction=nonstopmode -halt-on-error"
	@mkdir -p $(BUILDDIR)/html
	@mv $(BUILDDIR)/latex/SPARTA.pdf $(BUILDDIR)/html/Manual.pdf
	@echo "Build finished. The manual is in doc/html/Manual.pdf."

# A quick look at one page while editing it: pandoc renders each .rst on its
# own, with no cross-references and no theme, in about a second.  It is a
# preview, not the manual -- "make html" is the manual.
#
# It warns on a handful of equations that use \rm or \textnormal, which its
# MathML writer does not implement, and falls back to showing them as TeX.
# Sphinx and the PDF both typeset them correctly; there is nothing to fix in
# the sources.
fasthtml: $(FASTHTMLFILES)
	@echo "Preview built. The pages are in doc/fasthtml."

fasthtml-init:
	@if [ "$(HAS_PANDOC)" = "NO" ]; then \
	   echo "make fasthtml needs pandoc" 1>&2; exit 1; fi
	@mkdir -p fasthtml
	@cp -r $(RSTDIR)/JPG fasthtml/

# pandoc's MathML writer takes only em as a length unit, so a \hspace given
# in mm makes it give up on that equation.  Sphinx and MathJax both take it
# as written; substituting the thin space it stands for, only for the
# preview, keeps the equation readable without touching the source.
#
# fasthtml-init is an order-only prerequisite: as an ordinary one, make -j
# would start converting pages alongside the pandoc check rather than after
# it, and report "pandoc: command not found" for every page instead of
# saying once what is missing.
fasthtml/%.html: $(RSTDIR)/%.rst | fasthtml-init
	@sed -e 's/\\hspace{[0-9.]*mm}/\\,/g' $< > fasthtml/$*.tmp.rst
	@pandoc -s --mathml --metadata title="$*" -o $@ fasthtml/$*.tmp.rst
	@rm -f fasthtml/$*.tmp.rst
