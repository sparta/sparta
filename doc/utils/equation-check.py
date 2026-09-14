#!/usr/bin/env python3
"""Check each MathJax equation against the .tex source of the image it replaced.

The txt2html manual embedded pre-rendered images from doc/Eqs, and all but
one had its one-equation LaTeX source in the same directory.  The Sphinx
manual typesets that source directly, so the two can be compared exactly
rather than by eye: this pairs the Nth Eqs/ image on a page with the Nth
equation on the same page of the build and diffs the LaTeX, ignoring
whitespace and the wrappers each side supplies.

parity-check.py cannot do this.  An image contributes no text, so where one
stood it has nothing to compare against and has to allow whatever replaced
it; this is what closes that hole.

Usage:
    equation-check.py OLD_DOC_DIR NEW_HTML_DIR

OLD_DOC_DIR is a doc/ tree from before the migration -- the .html pages and
the Eqs directory beside them.  Exits non-zero on any mismatch that is not
listed in ALLOWED.

Two images are not equations and are not checked here: GS_list and PS_list
depicted tables of surface reactions and became reST tables.  Their .tex
sources are LaTeX tabulars, so they are compared by their reaction list
instead, which this also does.
"""
import argparse
import html as H
import pathlib
import re
import sys

# Differences from the .tex that are intended.  Each entry names the exact
# substitution, so the equation is still compared -- with the substitution
# applied to the .tex first.  Listing only a reason would switch the check
# off for the one equation it has been told to expect a difference in.
ALLOWED = {
    'impulsive_softsphere': [
        (r'\textlangle', r'\langle', ''),
        (r'\textrangle', r'\rangle',
         r'\textlangle and \textrangle need a LaTeX package MathJax does not '
         r'load, so it printed the macro names in red instead of the angle '
         r'brackets.  \langle and \rangle are the same characters.'),
    ],
}

# The one image with no .tex beside it; its LaTeX was written during the
# migration and checked against the published image by eye.
NO_SOURCE = {'lambda_old'}

# Depicted a table, not an equation; checked by reaction list below.
TABLES = {'GS_list', 'PS_list'}

IMG = re.compile(r'<img\b[^>]*?src\s*=\s*"(?:\./)?Eqs/([^"]+)"', re.I)
TAG = re.compile(r'<[^>]+>')
MATH_DIV = re.compile(r'<div class="math[^"]*">(.*?)</div>', re.S)


def tex_body(path):
    """The maths out of a standalone one-equation .tex document."""
    t = re.sub(r'%.*', '', path.read_text(errors='replace'))
    m = re.search(r'\\begin\{document\}(.*?)\\end\{document\}', t, re.S)
    if m:
        t = m.group(1)
    t = re.sub(r'\$\$(.*?)\$\$', r'\1', t, flags=re.S)
    return re.sub(r'^\s*\$(.*)\$\s*$', r'\1', t.strip(), flags=re.S)


def norm(s):
    """Compare the maths, not the spacing or the delimiters."""
    s = s.replace('\\[', ' ').replace('\\]', ' ')
    s = s.replace('\\(', ' ').replace('\\)', ' ')
    s = re.sub(r'\\(begin|end)\{(split|equation)\}', ' ', s)
    return re.sub(r'\s+', '', s)


def norm_reaction(r):
    r = re.sub(r'\s+', '', r)
    return r.replace('_{2}', '_2').replace('_{3}', '_3')


def tex_reactions(text):
    """The reactions in a LaTeX tabular, as a multiset."""
    body = text.split(r'\begin{tabular}')[-1].split(r'\end{tabular}')[0]
    return sorted(norm_reaction(r) for r in re.findall(r'\$(.*?)\$', body))


def html_tables(text):
    """Each <table> on a page, as its own multiset of reactions.

    Per table rather than per page: both reaction tables live on the same
    page, so scanning the whole page cannot tell whether a row landed in
    the right one.  A multiset rather than a set, so a duplicated or
    dropped repeat is caught.
    """
    out = []
    for m in re.finditer(r'<table\b.*?</table>', text, re.S | re.I):
        body = H.unescape(TAG.sub(' ', m.group(0)))
        out.append(sorted(norm_reaction(r)
                          for r in re.findall(r'\\\((.*?)\\\)', body)))
    return [t for t in out if t]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('old')
    ap.add_argument('new')
    args = ap.parse_args()
    old_dir, new_dir = pathlib.Path(args.old), pathlib.Path(args.new)

    checked = matched = 0
    problems = []
    skipped = []

    for page in sorted(p.name for p in old_dir.glob('*.html')):
        old = (old_dir / page).read_text(errors='replace')
        imgs = [i for i in IMG.findall(old)
                if pathlib.Path(i).stem not in TABLES]
        if not imgs:
            continue
        new = (new_dir / page).read_text(errors='replace')
        eqs = [H.unescape(TAG.sub('', m)) for m in MATH_DIV.findall(new)]

        for i, img in enumerate(imgs):
            stem = pathlib.Path(img).stem
            tex = old_dir / 'Eqs' / (stem + '.tex')
            if stem in NO_SOURCE:
                skipped.append((page, stem))
                continue
            if not tex.exists():
                problems.append((page, stem, 'no .tex in the old tree',
                                 'nothing to check this against; add it to '
                                 'NO_SOURCE with a reason if that is right'))
                continue
            checked += 1
            raw = tex_body(tex)
            for was, now, _why in ALLOWED.get(stem, ()):
                raw = raw.replace(was, now)
            want, got = norm(raw), norm(eqs[i]) if i < len(eqs) else ''
            if want == got:
                matched += 1
                if stem in ALLOWED:
                    skipped.append((page, stem + ' (agreed)'))
            else:
                problems.append((page, stem, want, got))

    # the two tables, by the reactions they list
    for stem in sorted(TABLES):
        tex = old_dir / 'Eqs' / (stem + '.tex')
        if not tex.exists():
            continue
        page = next((p.name for p in old_dir.glob('*.html')
                     if stem in (old_dir / p.name).read_text(errors='replace')), None)
        if not page:
            continue
        want = tex_reactions(tex.read_text(errors='replace'))
        tables = html_tables((new_dir / page).read_text(errors='replace'))
        checked += 1
        # exactly one table on the page must match this one exactly
        if want in tables:
            matched += 1
        else:
            near = min(tables, key=lambda x: len(set(x) ^ set(want)),
                       default=[])
            problems.append((page, stem,
                             f'{len(want)} reactions: ' + ', '.join(want),
                             f'closest table has {len(near)}: '
                             + ', '.join(near)))

    print(f'  {checked} equations and tables compared with their .tex source')
    print(f'  {matched} match')
    for page, stem in skipped:
        key = stem.split(' ')[0]
        why = next((w for _a, _b, w in ALLOWED.get(key, ()) if w), None)
        if key in NO_SOURCE:
            why = ('no .tex was published beside this image; its LaTeX was '
                   'written during the migration and checked against the '
                   'published image by eye')
        print(f'    ~ {page}: {stem} -- {why or "no reason recorded"}')
    for page, stem, want, got in problems:
        print(f'\n  MISMATCH {page}  Eqs/{stem}')
        print(f'    tex: {want[:200]}')
        print(f'    new: {got[:200]}')
    print('\n  RESULT: ' + ('EQUATIONS DIFFER' if problems else 'EQUATIONS OK'))
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
