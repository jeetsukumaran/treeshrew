#ifndef TREESHREW_TESTUTILS_HPP
#define TREESHREW_TESTUTILS_HPP

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include "../../src/dataio.hpp"
#include "../../src/genetree.hpp"
#include "../../src/utility.hpp"

std::vector<treeshrew::GeneTree *> get_trees(int argc, char * argv[]);

#endif
