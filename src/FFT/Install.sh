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

# all package files with no dependencies

for file in *.cpp *.h; do
  test -f ${file} && action $file
done

# no Makefile.package edits are needed for this package
# fftdata.h falls back to the bundled KISS FFT when no FFT library macro is
#   defined, so the package builds with no external library
# builds against MKL or FFTW instead set FFT_INC / FFT_PATH / FFT_LIB and the
#   -DFFT_MKL or -DFFT_FFTW3 switch in src/MAKE/Makefile.<machine>, which
#   every machine Makefile already provides
