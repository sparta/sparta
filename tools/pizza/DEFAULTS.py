# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@gmail.com, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# --------------
# --------------

# 3 variables used by Pizza.py
# to use TOOLS or SCRIPTS, edit and uncomment the line
# to use EXCLUDE, add to the existing list

# TOOLS = list of extra directories that contain Pizza.py tools
#   Pizza.py will load all *.py files as tools from TOOLS dirs, then pizza/src
#   this ordering means your tool can override a Pizza.py tool of the same name
# SCRIPTS = list of extra directories that contain Pizza.py scripts
#   Pizza.py will look in working dir, SCRIPTS dirs, then pizza/scripts
#   this ordering means your script can override a Pizza.py script of same name
# EXCLUDE = Python files to NOT load as tools when Pizza.py starts
#   typically done for auxiliary Python files that are not tools
#   any non-tool Python files from your TOOLS dirs should be added to list

#PIZZA_TOOLS = ["~/mystuff/new_pizza_tools"]
#PIZZA_SCRIPTS = ["~/mystuff/new_pizza_scripts"]
PIZZA_EXCLUDE = ["pizza", "DEFAULTS", "vizinfo"]

# --------------
# --------------

# Pathname for programs executed by various Pizza.py tools

# if you don't use a tool, it's settings can be ignored
# the default values are program names with no path
# to use a default value, the executable must therefore be in your path
# to change a default, uncomment and edit the PIZZA variable line

# --------------

# ImageMagick programs to manipulate image files
# DISPLAY = program to view GIF, PNG, SVG files
# tools that use it: rasmol, raster, svg
# CONVERT = program to convert one image format to another
# MONTAGE = program to stitch 2 images together
# tools that use it: image

#PIZZA_DISPLAY = "/usr/bin/display"
#PIZZA_CONVERT = "/usr/bin/convert"
#PIZZA_MONTAGE = "/usr/bin/montage"

# --------------

# GNUPLOT = the GnuPlot plotting package
# GNUTERM = terminal setting used by GnuPlot
# tools that use it: gnu

#PIZZA_GNUPLOT = "gnuplot"
#PIZZA_GNUTERM = "x11"
#PIZZA_GNUTERM = "aqua"                   # for Macs with Aquaterm installed

# --------------

# GUNZIP = program to uncompress gzipped files
# tools that use it: data dump log

#PIZZA_GUNZIP = "gunzip"

# --------------

# LABEL3D = program to put a label on a Raster3D image
# RENDER = the Raster3D visualization rendering engine
# tools that use it: raster

#PIZZA_LABEL3D = "label3d"
#PIZZA_RENDER = "render"

# --------------

# MATLAB = the MatLab numerical analysis and plotting package
# tools that use it: matlab

#PIZZA_MATLAB = "matlab -nosplash -nodesktop -nojvm"

# --------------

# RASMOL = the RasMol visualization package
# tools that use it: rasmol

#PIZZA_RASMOL = "rasmol"

# --------------

# VMD = the VMD visualization package
# tools that use it: vmd

#PIZZA_VMDNAME = "vmd"                 # good settings for a Linux box
#PIZZA_VMDDIR = "/usr/local/lib/vmd"
#PIZZA_VMDDEV = "win"
#PIZZA_VMDARCH = "LINUX"

#PIZZA_VMDNAME = "vmd"                 # good settings for a Mac
#PIZZA_VMDDIR = "/Applications/VMD\ 1.8.7.app/Contents/vmd"
#PIZZA_VMDDEV = "win"
#PIZZA_VMDARCH = "MACOSXX86"
