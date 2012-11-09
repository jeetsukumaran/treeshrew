treeshrew is a C++ (2011) program for the Bayesian estimation of alignments,
trees, and evolutionary parameters from short-read sequence data.

Build and install dependencies locally (run in project root):

    $ ./build-deps.sh

Before first building (run in project root):

    $ ./bootstrap.sh

Configuration (assuming 'build_deps.sh' was run):

    $ mkdir -p build/debug
    $ cd build/debug
    $ ../../cfgcommand.sh

Building and installing (assuming 'build_deps.sh' was run, and in build subdirectory):

    $ make
    $ make install

Following this recipe will result in all products being installed in
subdirectory 'installed/' of the build directory. To run the tests:

    $ source ../../treeshrew_deps_env.sh
    $ installed/opt/treeshrew/test/run-tests.py

If you get an error message regarding missing libraries during the build or
run, just make sure that the environmental variables are correctly set:

    $ source ../../treeshrew_deps_env.sh

##############################################################################
## treeshrew
##
## Copyright 2012 Jeet Sukumaran.
## All rights reserved.
##
## With code contributions from: Mark T. Holder.
## See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
