# Install/unInstall package files in SPARTA
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

# force rebuild of files with SPARTA_KOKKOS switch

KOKKOS_INSTALLED=0
if (test -e ../Makefile.package) then
  KOKKOS_INSTALLED=`grep DSPARTA_KOKKOS ../Makefile.package | wc -l`
fi

if (test $mode = 1) then
  if (test $KOKKOS_INSTALLED = 0) then
    touch ../accelerator_kokkos.h
  fi
elif (test $mode = 0) then
  if (test $KOKKOS_INSTALLED = 1) then
    touch ../accelerator_kokkos.h
  fi
fi

# list of files with optional dependcies

action collide_vss_kokkos.cpp
action collide_vss_kokkos.h
action compute_grid_kokkos.cpp
action compute_grid_kokkos.h
action compute_temp_kokkos.cpp
action compute_temp_kokkos.h
action compute_thermal_grid_kokkos.cpp
action compute_thermal_grid_kokkos.h
action comm_kokkos.cpp
action comm_kokkos.h
action domain_kokkos.cpp
action domain_kokkos.h
action fix_adapt_kokkos.cpp
action fix_adapt_kokkos.h
action fix_ave_grid_kokkos.cpp
action fix_ave_grid_kokkos.h
action fix_balance_kokkos.cpp
action fix_balance_kokkos.h
action fix_move_surf_kokkos.cpp
action fix_move_surf_kokkos.h
action geometry_kokkos.h
action grid_id_kokkos.cpp
action grid_kokkos.cpp
action grid_kokkos.h
action irregular_kokkos.cpp
action irregular_kokkos.h
action kokkos.cpp
action kokkos.h
action kokkos_base.h
action kokkos_copy.h
action kokkos_type.h
action math_extra_kokkos.h
action memory_kokkos.h
action modify_kokkos.cpp
action modify_kokkos.h
action particle_kokkos.cpp
action particle_kokkos.h
action rand_pool_wrap.cpp
action rand_pool_wrap.h
action surf_collide_diffuse_kokkos.cpp
action surf_collide_diffuse_kokkos.h
action surf_collide_specular_kokkos.cpp
action surf_collide_specular_kokkos.h
action surf_collide_vanish_kokkos.cpp
action surf_collide_vanish_kokkos.h
action surf_kokkos.cpp
action surf_kokkos.h
action update_kokkos.cpp
action update_kokkos.h
action kokkos_scan.cpp
action create_particles_kokkos.cpp
action create_particles_kokkos.h
action fix_emit_face_kokkos.cpp
action fix_emit_face_kokkos.h
action fix_grid_check_kokkos.cpp
action fix_grid_check_kokkos.h
action read_surf_kokkos.cpp
action read_surf_kokkos.h

# edit 2 Makefile.package files to include/exclude package info
# allow user to specify sed.  Useful on Mac OSX to specify SED=gsed
SED="${SED:-sed}"
if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    $SED -i -e 's/[^ \t]*kokkos[^ \t]* //g' ../Makefile.package
    $SED -i -e 's/[^ \t]*KOKKOS[^ \t]* //g' ../Makefile.package
    $SED -i -e 's|^PKG_INC =[ \t]*|&-DSPARTA_KOKKOS |' ../Makefile.package
#    $SED -i -e 's|^PKG_PATH =[ \t]*|&-L..\/..\/lib\/kokkos\/core\/src |' ../Makefile.package
    $SED -i -e 's|^PKG_CPP_DEPENDS =[ \t]*|&$(KOKKOS_CPP_DEPENDS) |' ../Makefile.package
    $SED -i -e 's|^PKG_LIB =[ \t]*|&$(KOKKOS_LIBS) |' ../Makefile.package
    $SED -i -e 's|^PKG_LINK_DEPENDS =[ \t]*|&$(KOKKOS_LINK_DEPENDS) |' ../Makefile.package
    $SED -i -e 's|^PKG_SYSINC =[ \t]*|&$(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) |' ../Makefile.package
    $SED -i -e 's|^PKG_SYSLIB =[ \t]*|&$(KOKKOS_LDFLAGS) |' ../Makefile.package
#    $SED -i -e 's|^PKG_SYSPATH =[ \t]*|&$(kokkos_SYSPATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    $SED -i -e '/CXX\ =\ \$(CC)/d' ../Makefile.package.settings
    $SED -i -e '/^include.*kokkos.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    $SED -i -e '4 i \CXX = $(CC)' ../Makefile.package.settings
    $SED -i -e '5 i \include ..\/..\/lib\/kokkos\/Makefile.kokkos' ../Makefile.package.settings
  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    $SED -i -e 's/[^ \t]*kokkos[^ \t]* //g' ../Makefile.package
    $SED -i -e 's/[^ \t]*KOKKOS[^ \t]* //g' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    $SED -i -e '/CXX\ =\ \$(CC)/d' ../Makefile.package.settings
    $SED -i -e '/^include.*kokkos.*$/d' ../Makefile.package.settings
  fi

fi
