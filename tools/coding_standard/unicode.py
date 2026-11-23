#!/usr/bin/env python3
# Utility for detecting and fixing non-ASCII characters in SPARTA
#
# Unicode homoglyphs and bidirectional control characters can be used to
# hide malicious changes in source code ("Trojan Source", CVE-2021-42574),
# so SPARTA files must be plain ASCII.  Common typographic characters
# (curly quotes, dashes, ...) are replaced with ASCII equivalents by -f.
# Invisible and bidirectional control characters are never fixed
# automatically: they must be inspected and removed by hand.
from __future__ import print_function
import sys

if sys.version_info.major < 3:
    sys.exit('This script must be run with Python 3.5 or later')

if sys.version_info.minor < 5:
    sys.exit('This script must be run with Python 3.5 or later')

import os
import glob
import yaml
import argparse
import shutil
import unicodedata

DEFAULT_CONFIG = """
recursive: true
include:
    - .
    - bench/**
    - cmake/**
    - data/**
    - doc/**
    - examples/**
    - python/**
    - src/**
    - tools/**
exclude:
    - lib/kokkos
patterns:
    - "*"
    - ".gitignore"
"""

# ASCII replacements applied by -f
REPLACEMENTS = {
    '\u2018': "'", '\u2019': "'", '\u201a': "'", '\u201b': "'",
    '\u201c': '"', '\u201d': '"', '\u201e': '"', '\u201f': '"',
    '\u2032': "'", '\u2033': '"',
    '\u2010': '-', '\u2011': '-', '\u2012': '-', '\u2013': '-',
    '\u2014': '-', '\u2015': '-', '\u2212': '-',
    '\u2026': '...', '\u00a0': ' ', '\u00b7': '*', '\u00d7': 'x',
    '\u03b8': 'theta', '\u03a3': 'Sum',
}

# never fixed automatically: invisible or reordering characters
def is_dangerous(c):
    o = ord(c)
    return (0x202a <= o <= 0x202e or 0x2066 <= o <= 0x2069 or
            0x200b <= o <= 0x200f or o in (0x061c, 0x2060, 0xfeff) or
            unicodedata.category(c) in ('Cf', 'Co', 'Cs') or
            (unicodedata.category(c) == 'Cc' and c not in '\t\n\r'))

def is_binary(path):
    with open(path, 'rb') as f:
        return b'\0' in f.read(8192)

def check_file(path):
    errors = []
    try:
        with open(path, 'r', encoding='UTF-8', newline='') as f:
            for lineno, line in enumerate(f, 1):
                for col, c in enumerate(line, 1):
                    if c in '\t\n\r' or ' ' <= c <= '~':
                        continue
                    errors.append((lineno, col, c))
    except UnicodeDecodeError:
        return {'errors': errors, 'encoding': 'unknown'}
    return {'errors': errors, 'encoding': 'UTF-8'}

def fix_file(path):
    newfile = path + ".modified"
    with open(path, 'r', encoding='UTF-8', newline='') as src:
        text = src.read()
    with open(newfile, 'w', encoding='UTF-8', newline='') as out:
        out.write(''.join(REPLACEMENTS.get(c, c) for c in text))
    shutil.copymode(path, newfile)
    shutil.move(newfile, path)

def check_folder(directory, config, fix=False, verbose=False):
    success = True
    files = set()

    for base_path in config['include']:
        for pattern in config['patterns']:
            path = os.path.join(directory, base_path, pattern)
            files.update(glob.glob(path, recursive=config['recursive']))
    for exclude in config['exclude']:
        files = [f for f in files if not os.path.normpath(f).startswith(os.path.join(directory, exclude))]

    for f in sorted(files):
        path = os.path.normpath(f)
        if os.path.islink(path) or not os.path.isfile(path) or is_binary(path):
            continue

        if verbose:
            print("Checking file:", path)

        result = check_file(path)

        if result['encoding'] == 'unknown':
            print("[Error] Not UTF-8 or ASCII text @ {}".format(path))
            success = False
            continue

        has_resolvable_errors = False
        has_manual_errors = False
        for lineno, col, c in result['errors']:
            name = unicodedata.name(c, 'UNNAMED')
            if is_dangerous(c):
                print("[Error] Invisible/bidi character U+{:04X} {} (fix manually) @ {}:{}:{}"
                      .format(ord(c), name, path, lineno, col))
                has_manual_errors = True
            else:
                print("[Error] Non-ASCII character U+{:04X} {} @ {}:{}:{}"
                      .format(ord(c), name, path, lineno, col))
                if c in REPLACEMENTS:
                    has_resolvable_errors = True
                else:
                    has_manual_errors = True

        if has_resolvable_errors:
            if fix:
                print("Applying automatic fixes to file:", path)
                fix_file(path)
            else:
                success = False
        if has_manual_errors:
            success = False

    return success

def main():
    parser = argparse.ArgumentParser(description='Utility for detecting and fixing non-ASCII characters in SPARTA')
    parser.add_argument('-c', '--config', metavar='CONFIG_FILE', help='location of a optional configuration file')
    parser.add_argument('-f', '--fix', action='store_true', help='replace common typographic characters with ASCII')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')
    parser.add_argument('DIRECTORY', help='directory that should be checked')
    args = parser.parse_args()
    spartadir = os.path.abspath(os.path.expanduser(args.DIRECTORY))

    if args.config:
        with open(args.config, 'r') as cfile:
            config = yaml.load(cfile, Loader=yaml.FullLoader)
    else:
        config = yaml.load(DEFAULT_CONFIG, Loader=yaml.FullLoader)

    if not check_folder(spartadir, config, args.fix, args.verbose):
        sys.exit(1)

if __name__ == "__main__":
    main()
