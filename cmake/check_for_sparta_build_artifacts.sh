#!/bin/bash
cd $1

installed_packages=$(make package-status | grep "Installed YES" 2>&1 > /dev/null; echo $?)
installed_style_files=$([ $(ls -lat style_*.h 2> /dev/null | wc -l) -ge 1 ] && true || false; echo $?)

if [ $installed_packages -eq 0 ]; then
    #if [ $2 == "y" ]; then
	#make no-all
    #else
	echo "CMake may not build properly! Please run 'make no-all' from $1."
    #fi
fi

if [ $installed_style_files -eq 0 ]; then
    #if [ $2 == "y" ]; then
	#rm -f style_*.h
    #else
	echo "CMake may not build properly! Please run 'make purge' from $1."
    #fi
fi
