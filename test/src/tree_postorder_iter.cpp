#include <vector>
#include <iostream>
#include <string>
#include <iomanip>
#include "../../src/dataio.hpp"
#include "../../src/genetree.hpp"
#include "../../src/utility.hpp"
#include "testutils.hpp"

int main(int argc, char * argv[]) {
    auto trees = get_trees(argc, argv);
    for (auto & tree : trees) {
        int postorder_count = 0;
        for (treeshrew::GeneTreeNode::postorder_iterator ndi = tree->postorder_begin(); ndi != tree->postorder_end(); ++ndi, ++postorder_count) {
            treeshrew::GeneTreeNode * nd = *ndi;
            std::cout << nd->get_label() << "\t" << std::setprecision(8) << nd->get_edge_length() << std::endl;
        }
    }
}


