# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# enforce using portable C locale
LC_ALL=C
export LC_ALL

# arg1 = file, arg2 = file it depends on

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

for file in *.cpp *.h; do
  test -f ${file} && action $file
done

# edit 2 Makefile.package files to include/exclude package info
SED="${SED:-sed}"
if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    $SED -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/lammps/src |' ../Makefile.package
    $SED -i -e '/^PKG_PATH *=/ {/-L..\/..\/lib\/lammps\/src/! s|$| -L../../lib/lammps/src -Wl,-rpath,$(abspath ../../lib/lammps/src)|}' ../Makefile.package
    $SED -i -e 's|^PKG_LIB =[ \t]*|& -llammps -ldl |' ../Makefile.package
  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
  $SED -i -e 's|[ \t]*-I\.\./\.\./lib/lammps/src||g' ../Makefile.package
  $SED -i -e 's|[ \t]*-L\.\./\.\./lib/lammps/src||g' ../Makefile.package
  $SED -i -e 's|[ \t]*-Wl,-rpath,\$(abspath \.\./\.\./lib/lammps/src)||g' ../Makefile.package
  $SED -i -e 's|[ \t]*-llammps||g' ../Makefile.package
  $SED -i -e 's|[ \t]*-ldl||g' ../Makefile.package
  fi

fi
