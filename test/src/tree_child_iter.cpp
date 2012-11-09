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
        int postorder_count = 0;
        for (treeshrew::GeneTreeNode::postorder_iterator ndi = tree->postorder_begin(); ndi != tree->postorder_end(); ++ndi, ++postorder_count) {
            auto node = *ndi;
            // if (node->is_leaf()) {
            //     continue;
            // }
            std::cout << node->get_label();
            for (treeshrew::GeneTreeNode::child_iterator chi = node->children_begin();
                    chi != node->children_end();
                    ++chi) {
                std::cout << "\t" << (*chi)->get_label();
            }
            std::cout << std::endl;
        }
    }
}
