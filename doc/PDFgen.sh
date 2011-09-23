#!/bin/csh
# generate a PDF version of Manual

txt2html -b *.txt

/home/sjplimp/tools/htmldoc/bin/htmldoc --title --toctitle "Table of Contents" --tocfooter ..i --toclevels 4 --header ... --footer ..1 --size letter --linkstyle plain --linkcolor blue -f Manual.pdf Manual.html Section_intro.html Section_start.html Section_commands.html Section_howto.html Section_example.html Section_perf.html Section_tools.html Section_modify.html Section_python.html Section_errors.html Section_history.html [a-z]*.html

txt2html *.txt
