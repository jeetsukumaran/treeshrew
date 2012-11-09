#include <vector>
#include <iostream>
#include <string>
#include "../../src/dataio.hpp"
#include "../../src/genetree.hpp"
#include "../../src/utility.hpp"
#include "testutils.hpp"

int main(int argc, char * argv[]) {
    auto trees = get_trees(argc, argv);
    for (auto & tree : trees) {
        int leaf_count = 0;
        for (treeshrew::GeneTreeNode::leaf_iterator ndi = tree->leaf_begin(); ndi != tree->leaf_end(); ++ndi, ++leaf_count) {
            treeshrew::GeneTreeNode * node = *ndi;
            std::cout << leaf_count + 1 << "\t" << node->get_label() << std::endl;
            if (leaf_count > 40) {
                std::cerr << "terminating: too many nodes" << std::endl;
                exit(1);
            }
        }
    }
}
