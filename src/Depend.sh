# Depend.sh = Install/unInstall files due to package dependencies
# this script is invoked after any package is installed/uninstalled

# enforce using portable C locale

LC_ALL=C
export LC_ALL

# all parent/child package dependencies should be listed below
# parent package = has files that files in another package derive from
# child package = has files that derive from files in another package

# update child packages that depend on the parent,
#   but only if the child package is already installed
# this is necessary to insure the child package installs
#   only child files whose parent package files are now installed
# decisions on (un)installing individual child files are made by
#   the Install.sh script in the child package

# depend function: arg = child-package
# checks if child-package is installed, if not just return
# otherwise invoke update of child package via its Install.sh

depend () {
  cd $1
  installed=0
  for file in *.cpp *.h; do
    if (test -e ../$file) then
      installed=1
    fi
  done

  cd ..
  if (test $installed = 0) then
    return
  fi

  echo "  updating package $1"
  if (test -e $1/Install.sh) then
    cd $1; /bin/sh Install.sh 2; cd ..
  else
    cd $1; /bin/sh ../Install.sh 2; cd ..
  fi
}

# add one if statement per parent package
# add one depend() call per child package that depends on that parent

#if (test $1 = "FFT") then
#  depend KOKKOS
#fi

