#!/bin/sh
# set -x

if test "$1" == "--force"
then
    echo "Forcing rebuilding of all dependencies"
    unset BEAGLE_ROOT
    unset BEAGLE_PREFIX
    unset NCL_ROOT
    unset NCL_PREFIX
    unset GSL_ROOT
    unset GSL_PREFIX
fi

env_filename='treeshrew_deps_env.sh'
echo '#!/bin/sh' > "${env_filename}"
TREESHREW_ROOT="$PWD"
echo "export TREESHREW_ROOT=${TREESHREW_ROOT}" >> "${env_filename}"

LIB_ROOT="$TREESHREW_ROOT/ext"
echo "export LIB_ROOT=${LIB_ROOT}" >> "${env_filename}"
if ! test -d $LIB_ROOT
then
    mkdir $LIB_ROOT
fi

echo $OSTYPE | grep darwin >/dev/null
is_linux=$?

export MAKE="make"

################################################################################
# GSL
################################################################################
if test -z $GSL_PREFIX
then
    if test -z $GSL_ROOT
    then
        cd $LIB_ROOT
        if ! test -d gsl-1.15
        then
            tar xvf ../gsl/gsl-1.15.tar.gz
        fi
        GSL_ROOT="$LIB_ROOT/gsl-1.15"
    fi
    echo "export GSL_ROOT=${GSL_ROOT}" >> "${env_filename}"
    cd "$GSL_ROOT" || exit 1
    mkdir build
    autoreconf
    sh "./autogen.sh" || exit 1
    cd build
    export GSL_PREFIX="${GSL_ROOT}/installed"
    ../configure --prefix="${GSL_PREFIX}" || exit 1
    ${MAKE} || exit 1
    ${MAKE} install || exit 1
    if test ${is_linux} -eq 0
    then
        export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${GSL_PREFIX}/lib"
    else
        export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GSL_PREFIX}/lib"
    fi
    ${MAKE} check || exit
    cd $TREESHREW_ROOT
else
    if test ${is_linux} -eq 0
    then
        export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${GSL_PREFIX}/lib"
    else
        export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GSL_PREFIX}/lib"
    fi
fi
echo "export GSL_PREFIX=${GSL_PREFIX}" >> "${env_filename}"
if test ${is_linux} -eq 0
then
    echo 'export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${GSL_PREFIX}/lib"' >> "${env_filename}"
else
    echo 'export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GSL_PREFIX}/lib"' >> "${env_filename}"
fi
export PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:${GSL_PREFIX}/lib/pkgconfig"
echo 'export PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:${GSL_PREFIX}/lib/pkgconfig"' >> "${env_filename}"


################################################################################
# BEAGLE
################################################################################
if test -z $BEAGLE_PREFIX
then
    if test -z $BEAGLE_ROOT
    then
        cd $LIB_ROOT
        if ! test -d beagle-lib
        then
            svn checkout http://beagle-lib.googlecode.com/svn/trunk beagle-lib || exit 1
        fi
        BEAGLE_ROOT="$LIB_ROOT/beagle-lib"
    fi
    echo "export BEAGLE_ROOT=${BEAGLE_ROOT}" >> "${env_filename}"
    cd "$BEAGLE_ROOT" || exit 1
    mkdir build
    if test ${is_linux} -eq 0
    then
        perl -pi -e 's/\-march=native//g' ./configure.ac
    fi
    #rxtext.py -f '\-march=native' -r '' ./configure.ac
    sh "./autogen.sh" || exit 1
    cd build
    export BEAGLE_PREFIX="${BEAGLE_ROOT}/installed"
    ../configure --prefix="${BEAGLE_PREFIX}" || exit 1
    ${MAKE} || exit 1
    ${MAKE} install || exit 1
    if test ${is_linux} -eq 0
    then
        export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${BEAGLE_PREFIX}/lib"
    else
        export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${BEAGLE_PREFIX}/lib"
    fi
    ${MAKE} check || exit
    cd $TREESHREW_ROOT
else
    if test ${is_linux} -eq 0
    then
        export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${BEAGLE_PREFIX}/lib"
    else
        export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${BEAGLE_PREFIX}/lib"
    fi
fi
echo "export BEAGLE_PREFIX=${BEAGLE_PREFIX}" >> "${env_filename}"
if test ${is_linux} -eq 0
then
    echo 'export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${BEAGLE_PREFIX}/lib"' >> "${env_filename}"
else
    echo 'export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${BEAGLE_PREFIX}/lib"' >> "${env_filename}"
fi
export PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:${BEAGLE_PREFIX}/lib/pkgconfig"
echo 'export PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:${BEAGLE_PREFIX}/lib/pkgconfig"' >> "${env_filename}"



################################################################################
# NCL
################################################################################
if test -z $NCL_PREFIX
then
    if test -z $NCL_ROOT
    then
        cd $LIB_ROOT
        if ! test -d v2.1
        then
            svn checkout https://ncl.svn.sourceforge.net/svnroot/ncl/branches/v2.1 || exit 1
        fi
        NCL_ROOT="$LIB_ROOT/v2.1"
    fi
    echo "export NCL_ROOT=${NCL_ROOT}" >> "${env_filename}"
    cd "$NCL_ROOT" || exit 1
    mkdir build
    sh "./bootstrap.sh" || exit 1
    cd build
    export NCL_PREFIX="${NCL_ROOT}/build/installed"
    CXXFLAGS="-std=c++11" ../configure --prefix="${NCL_PREFIX}" || exit 1
    CXXFLAGS="-std=c++11" ${MAKE} || exit 1
    ${MAKE} install || exit 1
    ${MAKE} install check || exit 1
    if test ${is_linux} -eq 0
    then
        export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${NCL_PREFIX}/lib/ncl"
    else
        export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${NCL_PREFIX}/lib/ncl"
    fi
    ${MAKE} check || exit
    cd $TREESHREW_ROOT
else
    if test ${is_linux} -eq 0
    then
        export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${NCL_PREFIX}/lib/ncl"
    else
        export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${NCL_PREFIX}/lib/ncl"
    fi
fi
echo "export NCL_PREFIX=${NCL_PREFIX}" >> "${env_filename}"
if test ${is_linux} -eq 0
then
    echo 'export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${NCL_PREFIX}/lib/ncl"'>> "${env_filename}"
else
    echo 'export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${NCL_PREFIX}/lib/ncl"'>> "${env_filename}"
fi

################################################################################
# Convenience
################################################################################

cat << 'EOF' > cfgcommand.sh
#! /bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if test -z $NCL_DIR
then
    source ${SCRIPT_DIR}/treeshrew_deps_env.sh
fi
if [[ $1 == "release" ]]
then
    CXXFLAGS="-O3 -Wall" ${SCRIPT_DIR}/configure --prefix=`pwd`/installed "--with-ncl=${NCL_PREFIX}" "--with-beagle=${BEAGLE_PREFIX}"
    echo
    echo ---
    echo Release configuration complete.
elif [[ $1 == "profile" ]]
then
    CXXFLAGS="-O3 -pg -g -Wall" ${SCRIPT_DIR}/configure --prefix=`pwd`/installed "--with-ncl=${NCL_PREFIX}" "--with-beagle=${BEAGLE_PREFIX}"
    echo
    echo ---
    echo Release configuration complete.
else
    CXXFLAGS="-g -O0 -Wall -DTREESHREW_DEBUG_PRINTING -DDEBUGGING_MODE" ${SCRIPT_DIR}/configure --prefix=`pwd`/installed "--with-ncl=${NCL_PREFIX}" "--with-beagle=${BEAGLE_PREFIX}"
    echo
    echo ---
    echo Debugging configuration complete.
fi
EOF
chmod a+x cfgcommand.sh

echo
echo ---
echo "Run './bootstrap.sh' first if this is a clean checkout (or if you do not yet have a './configure' file in the project root directory."
echo "Configure by invoking 'cfgcommand.sh debug', 'cfgcommand.sh profile', or 'cfgcommand.sh release' (script should be in the project root directory, but invoked from the build directory)."

