#!/usr/bin/env python3
"""Check that every internal link in a built manual resolves.

Both halves of a link matter.  The file has to be there, and so does the
anchor: an anchor name is part of a URL, renaming one breaks an inbound
link just as surely as deleting the page, and it is the failure this
migration was most at risk of.

This runs against the build alone, with no reference to the txt2html
manual, so unlike parity-check.py and equation-check.py it keeps working
after the migration is over.  It is what CI runs.

Usage:
    link-check.py HTML_DIR
"""
import argparse
import html
import pathlib
import re
import sys

# <a name="..."> is the anchor form the conversion kept for the names the
# old manual published; id="..." on any element is Sphinx's own.  Scoped to
# inside a tag, because the manual's prose contains things like
# 'the default mixture has an ID = "all"', which is text, not an anchor.
NAME = re.compile(r'<a\b[^>]*?\bname\s*=\s*"([^"]+)"', re.I)
TAG = re.compile(r'<[a-zA-Z][^>]*>')
ID = re.compile(r'\bid\s*=\s*"([^"]+)"', re.I)
HREF = re.compile(r'<a\b[^>]*?href\s*=\s*"([^"]*)"', re.I)

# built by the pdf target, into this directory; see the note below
PDF = 'Manual.pdf'


def anchors(text):
    found = set(NAME.findall(text))
    for tag in TAG.finditer(text):
        found |= set(ID.findall(tag.group(0)))
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('html_dir')
    args = ap.parse_args()
    root = pathlib.Path(args.html_dir)
    if not root.is_dir():
        print(f'no such directory: {root}', file=sys.stderr)
        return 2

    have = {p.name: anchors(p.read_text(errors='replace'))
            for p in root.glob('*.html')}

    bad = []
    unbuilt = set()
    for page in sorted(root.glob('*.html')):
        text = page.read_text(errors='replace')
        for m in HREF.finditer(text):
            href = html.unescape(m.group(1))
            if href.startswith(('http', 'mailto', 'ftp', 'javascript:')):
                continue
            target, _, frag = href.partition('#')
            if target and not (root / target).exists():
                # "make pdf" is a separate target with a LaTeX toolchain
                # behind it, so the PDF is legitimately absent from an
                # html-only build.  Named exactly, and only at the path the
                # pdf target writes to -- an earlier version of this check
                # skipped Manual.pdf wherever it appeared, which is why
                # nobody noticed the pdf target was writing it one
                # directory above the page that links to it.
                if target == PDF:
                    unbuilt.add(target)
                    continue
                bad.append(f'{page.name} -> {target} (no such file)')
                continue
            dest = target or page.name
            if frag and dest in have and frag not in have[dest]:
                bad.append(f'{page.name} -> {href} (no such anchor)')

    for b in sorted(set(bad)):
        print(f'  {b}')
    for u in sorted(unbuilt):
        print(f'  note: {u} is not in this build; "make pdf" puts it here')
    n = len(set(bad))
    print(f'  {n} dead internal link(s) or anchor(s) in {len(have)} pages')
    return 1 if n else 0


if __name__ == '__main__':
    sys.exit(main())
