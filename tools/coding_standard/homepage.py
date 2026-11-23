#!/usr/bin/env python3
# Utility for detecting incorrect SPARTA homepage URLs
#
# Written by Richard Berger (Temple University)
from __future__ import print_function
import sys

if sys.version_info.major < 3:
    sys.exit('This script must be run with Python 3.5 or later')

if sys.version_info.minor < 5:
    sys.exit('This script must be run with Python 3.5 or later')

import os
import glob
import re
import yaml
import argparse
import shutil

DEFAULT_CONFIG = """
recursive: true
include:
    - .
    - bench
    - cmake/**
    - doc
    - python/**
    - src/**
    - tools/**
patterns:
    - "*.c"
    - "*.cmake"
    - "*.cpp"
    - "*.h"
    - "*.md"
    - "*.py"
    - "*.sh"
    - ".gitignore"
    - "CMakeLists.txt"
    - "README"
    - "requirements.txt"
"""

def check_homepage(f):
    pattern = re.compile(r'.*(sparta.sandia.gov|www.cs.sandia.gov/~sjplimp/sparta).*')
    lineno = 1
    errors = set()

    for line in f:
        if pattern.match(line):
            errors.add(lineno)
        lineno += 1

    return errors

def check_email(f):
    pattern = re.compile(r'.*sjplimp@sandia.gov.*')
    lineno = 1
    errors = set()

    for line in f:
        if pattern.match(line):
            errors.add(lineno)
        lineno += 1

    return errors

def check_file(path):
    if path.find('homepage.py') >= 0: return { 'homepage_errors' : '', 'email_errors' : '' }
    encoding = 'UTF-8'
    homepage_errors = set()
    email_errors = set()
    try:
        with open(path, 'r') as f:
            homepage_errors = check_homepage(f)
    except UnicodeDecodeError:
        encoding = 'ISO-8859-1'
        try:
            with open(path, 'r', encoding=encoding) as f:
                homepage_errors = check_homepage(f)
        except Exception:
            encoding = 'unknown'
    try:
        with open(path, 'r') as f:
            email_errors = check_email(f)
    except UnicodeDecodeError:
        encoding = 'ISO-8859-1'
        try:
            with open(path, 'r', encoding=encoding) as f:
                email_errors = check_email(f)
        except Exception:
            encoding = 'unknown'

    return {
        'homepage_errors': homepage_errors,
        'email_errors': email_errors,
        'encoding': encoding
    }

def fix_file(path, check_result):
    if path.find('homepage.py') >= 0: return
    newfile = path + ".modified"
    with open(newfile, 'w', encoding='UTF-8') as out:
        with open(path, 'r', encoding=check_result['encoding']) as src:
            for line in src:
                newline = line.replace("https://sparta.sandia.gov","https://sparta.github.io")
                newline = newline.replace("http://sparta.sandia.gov","https://sparta.github.io")
                newline = newline.replace("sparta.sandia.gov","sparta.github.io")
                newline = newline.replace("www.cs.sandia.gov/~sjplimp/sparta.html","sparta.github.io")
                newline = newline.replace("sjplimp@sandia.gov","sjplimp@gmail.com")
                print(newline, end='', file=out)
    shutil.copymode(path, newfile)
    shutil.move(newfile, path)

def check_folder(directory, config, fix=False, verbose=False):
    success = True
    files = []

    for base_path in config['include']:
        for pattern in config['patterns']:
            path = os.path.join(directory, base_path, pattern)
            files += glob.glob(path, recursive=config['recursive'])

    for f in files:
        path = os.path.normpath(f)
        if os.path.islink(path):
            continue

        if verbose:
            print("Checking file:", path)

        result = check_file(path)
        has_resolvable_errors = False

        for lineno in result['homepage_errors']:
            print("[Error] Incorrect SPARTA homepage @ {}:{}".format(path, lineno))
            has_resolvable_errors = True

        for lineno in result['email_errors']:
            print("[Error] Incorrect Developer email @ {}:{}".format(path, lineno))
            has_resolvable_errors = True

        if has_resolvable_errors:
            if fix:
                print("Applying automatic fixes to file:", path)
                fix_file(path, result)
            else:
                success = False

    return success

def main():
    parser = argparse.ArgumentParser(description='Utility for detecting and correct the SPARTA homepage URL and developer email in SPARTA source files')
    parser.add_argument('-c', '--config', metavar='CONFIG_FILE', help='location of a optional configuration file')
    parser.add_argument('-f', '--fix', action='store_true', help='automatically fix URLs')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')
    parser.add_argument('DIRECTORY', help='directory (or file) that should be checked')
    args = parser.parse_args()
    spartadir = os.path.abspath(os.path.expanduser(args.DIRECTORY))

    if args.config:
        with open(args.config, 'r') as cfile:
            config = yaml.load(cfile, Loader=yaml.FullLoader)
    else:
        config = yaml.load(DEFAULT_CONFIG, Loader=yaml.FullLoader)

    if os.path.isdir(spartadir):
        if not check_folder(spartadir, config, args.fix, args.verbose):
           sys.exit(1)
    else:
        success = True
        path = os.path.normpath(spartadir)

        if args.verbose:
            print("Checking file:", path)

        result = check_file(path)

        has_resolvable_errors = False

        for lineno in result['homepage_errors']:
            print("[Error] Incorrect SPARTA homepage @ {}:{}".format(path, lineno))
            has_resolvable_errors = True

        if has_resolvable_errors:
            if args.fix:
                print("Applying automatic fixes to file:", path)
                fix_file(path, result)
            else:
                success = False

        if not success:
            sys.exit(1)

if __name__ == "__main__":
    main()
