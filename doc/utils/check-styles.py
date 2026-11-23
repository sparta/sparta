#!/usr/bin/env python3
"""Check the manual's command and style listings against the source.

Modelled on the check-styles script LAMMPS runs on its manual.  SPARTA
registers each style with a macro in src/*.h or src/<PACKAGE>/*.h:

    ComputeStyle(grid,ComputeGrid)
    ComputeStyle(grid/kk,ComputeGridKokkos)

and Section_commands.rst lists them: one table of input-script commands
("Individual commands") and one table per style category ("Fix styles",
"Compute styles", ...).  Each entry is a :doc: link whose text is the
name, followed by "(k)" when a KOKKOS variant exists.

For every category this checks, in both directions:

  - every style the source registers is listed in that category's table,
    and every listed style is registered -- by category, so that
    "compute grid" being listed says nothing about "dump grid";
  - "(k)" is shown exactly when a /kk variant is registered;
  - the page each entry links to exists.

The commands table is checked against the commands the source defines:
every CommandStyle() plus the commands src/input.cpp dispatches itself.
Its "dump image"-style entries name dump styles and are checked as such;
"dump vtk" stands for the family grid/vtk, particle/vtk and surf/vtk.

Usage:  check-styles.py SRC_DIR DOC_SRC_DIR
Exits 1 and lists every discrepancy if there is any.
"""

import glob
import os
import re
import sys

# Section heading in Section_commands.rst -> the macro prefix in src/.
SECTIONS = {
    'Fix styles': 'Fix',
    'Compute styles': 'Compute',
    'Collide styles': 'Collide',
    'Surface collide styles': 'SurfCollide',
    'Surface reaction styles': 'SurfReact',
    'Reaction styles': 'React',
    'Region styles': 'Region',
    'Dump styles': 'Dump',
}
COMMANDS_SECTION = '3.5 Individual commands'

STYLE_RE = re.compile(r'^(\w+)Style\(([^,)]+),', re.M)
ENTRY_RE = re.compile(r':doc:`((?:[^`<\\]|\\.)+?)\s*<([^>]+)>`')
HEADING_RE = re.compile(r'^(\S[^\n]*)\n([-=~^"]+)\n', re.M)


def registered(src_dir):
    """{category: {name: has_kokkos_variant}} from the Style() macros."""
    styles = {}
    for path in glob.glob(os.path.join(src_dir, '*.h')) + \
            glob.glob(os.path.join(src_dir, '*', '*.h')):
        if os.path.basename(os.path.dirname(path)) in ('MAKE', 'STUBS'):
            continue
        with open(path, errors='replace') as f:
            for category, name in STYLE_RE.findall(f.read()):
                name = name.strip()
                base = name[:-3] if name.endswith('/kk') else name
                cat = styles.setdefault(category, {})
                cat[base] = cat.get(base, False) or name.endswith('/kk')
    with open(os.path.join(src_dir, 'input.cpp')) as f:
        code = re.sub(r'//[^\n]*', '', f.read())
        for name in re.findall(r'strcmp\(command,"(\w+)"\)', code):
            styles.setdefault('Command', {}).setdefault(name, False)
    return styles


def sections(rst):
    """{heading: body} for every section of a reST file."""
    heads = [(m.group(1).strip(), m.start(), m.end())
             for m in HEADING_RE.finditer(rst)
             if len(m.group(2)) >= len(m.group(1).strip())]
    return {h: rst[end:(heads[i + 1][1] if i + 1 < len(heads) else len(rst))]
            for i, (h, _start, end) in enumerate(heads)}


def entries(body):
    """[(name, kokkos_marked, page)] from the :doc: links in a table."""
    table = '\n'.join(l for l in body.splitlines() if l.startswith('|'))
    out = []
    for text, page in ENTRY_RE.findall(table):
        text = text.replace('\\', '').strip()
        marked = text.endswith('(k)')
        out.append((text[:-3].strip() if marked else text, marked, page))
    return out


def check(src_dir, doc_dir):
    styles = registered(src_dir)
    with open(os.path.join(doc_dir, 'Section_commands.rst')) as f:
        secs = sections(f.read())
    problems = []

    def compare(label, category, listed):
        have = styles.get(category, {})
        seen = {}
        for name, marked, page in listed:
            if name in seen:
                problems.append('%s: %s is listed twice' % (label, name))
            seen[name] = marked
            if not os.path.exists(os.path.join(doc_dir, page + '.rst')):
                problems.append('%s: %s links to %s.rst, which does not exist'
                                % (label, name, page))
        for name in sorted(set(have) - set(seen)):
            problems.append('%s: %s is registered in src/ but not listed'
                            % (label, name))
        for name in sorted(set(seen) - set(have)):
            problems.append('%s: %s is listed but src/ does not define it'
                            % (label, name))
        for name in sorted(set(seen) & set(have)):
            if seen[name] and not have[name]:
                problems.append('%s: %s is marked (k) but has no /kk variant'
                                % (label, name))
            if have[name] and not seen[name]:
                problems.append('%s: %s has a /kk variant but is not marked '
                                '(k)' % (label, name))

    for heading, category in SECTIONS.items():
        if heading not in secs:
            if styles.get(category):
                problems.append('%s: no such section in Section_commands.rst, '
                                'so none of the %d %s styles are listed'
                                % (heading, len(styles[category]), category))
            continue
        compare(heading, category, entries(secs[heading]))

    if COMMANDS_SECTION not in secs:
        problems.append('%s: no such section' % COMMANDS_SECTION)
    else:
        commands, dumps = [], []
        for name, marked, page in entries(secs[COMMANDS_SECTION]):
            if name.startswith('dump ') and ' ' in name:
                style = name.split(None, 1)[1]
                dumps.append((style, marked, page))
            else:
                commands.append((name, marked, page))
        compare(COMMANDS_SECTION, 'Command', commands)
        have = styles.get('Dump', {})
        for style, _marked, _page in dumps:
            # "dump vtk" names the page for grid/vtk, particle/vtk, surf/vtk
            family = [n for n in have if n == style or n.endswith('/' + style)]
            if not family:
                problems.append('%s: "dump %s" is listed but src/ defines no '
                                'dump style %s' % (COMMANDS_SECTION, style,
                                                   style))
    return styles, problems


def main(argv):
    if len(argv) != 3:
        print(__doc__.strip().splitlines()[-3], file=sys.stderr)
        return 2
    styles, problems = check(argv[1], argv[2])
    total = sum(len(v) for v in styles.values())
    if problems:
        print('check-styles: %d problem(s) in the listings of %d commands '
              'and styles:' % (len(problems), total))
        for p in problems:
            print('  ' + p)
        return 1
    print('check-styles: %d commands and styles, every one listed once, '
          'with the right (k) marker and a page that exists' % total)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
