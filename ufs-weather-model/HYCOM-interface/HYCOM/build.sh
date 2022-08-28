#!/bin/ksh -l

set -eu

# Formatting Variables
HR="----------------------------------------"
HR+="----------------------------------------"

# Script Help
help="NAME\n"
help+="\t$0 --- HYCOM build script\n\n"
help+="USAGE\n"
help+="\t$0 [ options ]\n\n"
help+="DESCRIPTION\n"
help+="\tThis script builds or cleans the HYCOM object files and "
help+="executable.\n\n"
help+="OPTIONS\n"
help+="\t-b,\tBuild object files and executable (default).\n"
help+="\t-a,\tBuild architecture.\n"
help+="\t-t,\tBuild type.\n"
help+="\t-e,\tBuild extras.\n"
help+="\t-n,\tBuild with NUOPC cap.\n"
help+="\t-i,\tNUOPC installation directory.\n"
help+="\t-c,\tClean object files and executable.\n"
help+="\t-v,\tIncrease output.\n"
help+="\t-h,\tHelp information.\n"

# Script Usage
usage="ba:t:e:ni:cvh"

# Default Options
verbosity=0
build=true clean=false
build_arch="intelsse-impi-sm-relo"
build_type="mpi"
build_extras=""
nuopc=false
nuopc_dir="NUOPC"
nuopc_install=false
nuopc_install_dir="."

# Process Command Line Options
while getopts "$usage" optchar ; do
  case $optchar in
    b ) ;;
    a ) build_arch=$OPTARG ;;
    t ) build_type=$OPTARG ;;
    e ) build_extras=$OPTARG ;;
    n ) nuopc=true ;;
    i ) nuopc_install=true nuopc_install_dir=$OPTARG ;;
    c ) clean=true build=false ;;
    v ) verbosity=$((verbosity+1)) ;;
    h ) print $help; exit 0;;
    * ) print >&2 "See Help: $0 -h"; exit 1;;
  esac
done

if $nuopc; then
  if $nuopc_install; then
    mkdir -p $nuopc_install_dir
  fi
  build_extras=${build_extras}" -DEOS_SIG2 -DEOS_17T -DESPC_COUPLE"
else
  build_extras=${build_extras}" -DEOS_SIG2 -DEOS_17T"
fi

nuopc_dir=$(readlink -f $nuopc_dir)
nuopc_install_dir=$(readlink -f $nuopc_install_dir)

# Print Command Line Options
if [ $verbosity -ge 1 ]; then
  print "$HR"
  print "Compile Script Settings"
  print "$HR"
  print "\tVerbosity:         $verbosity"
  print "\tBuild:             $build"
  print "\tBuild Arch:        $build_arch"
  print "\tBuild Type:        $build_type"
  print "\tBuild Extras:      $build_extras"
  print "\tClean:             $clean"
  print "\tNUOPC:             $nuopc"
  print "\tNUOPC Dir:         $nuopc_dir"
  print "\tNUOPC Install:     $nuopc_install"
  print "\tNUOPC Install Dir: $nuopc_install_dir"
  print "$HR"
fi

if $nuopc; then
  if $build; then
    (cd ${nuopc_dir} && make nuopc \
      ARCH="${build_arch}" TYPE="${build_type}" \
      CPP_EXTRAS="${build_extras}")
  fi
  if $clean; then
    (cd ${nuopc_dir} && make nuopcdistclean \
      ARCH="${build_arch}" TYPE="${build_type}" \
      CPP_EXTRAS="${build_extras}")
  fi
else
  if $build; then
    make -f Makefile hycom \
      ARCH="${build_arch}" TYPE="${build_type}" \
      CPP_EXTRAS="${build_extras}"
  fi
  if $clean; then
    make -f Makefile clean \
      ARCH="${build_arch}" TYPE="${build_type}" \
      CUSE_FLAG="${build_extras}"
  fi
fi

# Edit Executable for specified architectures
if [[ $build_arch == "Asp5" || $build_arch == "sp5" ]]; then
  ldedit -bdatapsize=64K -bstackpsize=64K hycom
fi
if [[ $build_arch == "Asp6" || $build_arch == "sp6" ]]; then
  ldedit -bdatapsize=64K -bstackpsize=64K hycom
fi
if [[ $build_arch == "Asp6-nofl" || $build_arch == "sp6-nofl" ]]; then
  ldedit -bdatapsize=64K -bstackpsize=64K hycom
fi
