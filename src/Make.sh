# Make.sh = update Makefile.lib, Makefile.shlib, Makefile.list 
#           or style_*.h files
# Syntax: sh Make.sh style
#         sh Make.sh Makefile.lib
#         sh Make.sh Makefile.shlib
#         sh Make.sh Makefile.list

# function to create one style_*.h file
# must whack *.d files that depend on style_*.h file,
# else Make will not recreate them

style () {
  list=`grep -l $1 $2*.h`
  if (test -e style_$3.tmp) then
    rm -f style_$3.tmp
  fi
  for file in $list; do
    qfile="\"$file\""
    echo "#include $qfile" >> style_$3.tmp
  done
  if (test ! -e style_$3.tmp) then
    rm -f style_$3.h
    touch style_$3.h
  elif (test ! -e style_$3.h) then
    mv style_$3.tmp style_$3.h
    rm -f Obj_*/$4.d
    rm -f Obj_*/sparta.d
  elif (test "`diff --brief style_$3.h style_$3.tmp`" != "") then
    mv style_$3.tmp style_$3.h
    rm -f Obj_*/$4.d
    rm -f Obj_*/sparta.d
  else
    rm -f style_$3.tmp
  fi
}

# create individual style files
# called by "make machine"

if (test $1 = "style") then

  style COLLIDE_CLASS   collide_    collide    input
  style COMMAND_CLASS   ""          command    input
  style COMPUTE_CLASS   compute_    compute    modify
  style DUMP_CLASS      dump_       dump       output
  style FIX_CLASS       fix_        fix        modify
  style REACT_CLASS     react_      react      input
  style REGION_CLASS    region_     region     domain
  style SURF_COLLIDE_CLASS surf_collide_ surf_collide surf
  style SURF_REACT_CLASS   surf_react_   surf_react   surf

# edit Makefile.lib
# called by "make makelib"
# use current list of *.cpp and *.h files in src dir w/out main.cpp

elif (test $1 = "Makefile.lib") then

  list=`ls -1 *.cpp | sed s/^main\.cpp// | tr "[:cntrl:]" " "`
  sed -i -e "s/SRC =	.*/SRC =	$list/" Makefile.lib
  list=`ls -1 *.h | tr "[:cntrl:]" " "`
  sed -i -e "s/INC =	.*/INC =	$list/" Makefile.lib

# edit Makefile.lib, for creating non-shared lib
# called by "make makelib"
# use current list of *.cpp and *.h files in src dir w/out main.cpp

elif (test $1 = "Makefile.shlib") then

  list=`ls -1 *.cpp | sed s/^main\.cpp// | tr "[:cntrl:]" " "`
  sed -i -e "s/SRC =	.*/SRC =	$list/" Makefile.shlib
  list=`ls -1 *.h | tr "[:cntrl:]" " "`
  sed -i -e "s/INC =	.*/INC =	$list/" Makefile.shlib

# edit Makefile.list
# called by "make makelist"
# use current list of *.cpp and *.h files in src dir

elif (test $1 = "Makefile.list") then

  list=`ls -1 *.cpp | tr "[:cntrl:]" " "`
  sed -i -e "s/SRC =	.*/SRC =	$list/" Makefile.list
  list=`ls -1 *.h | tr "[:cntrl:]" " "`
  sed -i -e "s/INC =	.*/INC =	$list/" Makefile.list

fi
