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

    $ source ../../treeshrew_deps_env.sh
    $ make
    $ make install

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
