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
        treeshrew::treeio::write_newick(tree, std::cout);
        std::cerr << "---" << std::endl;
        treeshrew::GeneTreeNode * root = tree->get_seed_node();
        std::cerr << "root is " <<  (void*)root << ": " << root->get_index() << std::endl;
        int postorder_count = 0;
        for (treeshrew::GeneTreeNode::postorder_iterator ndi = tree->postorder_begin(); ndi != tree->postorder_end(); ++ndi, ++postorder_count) {
            treeshrew::GeneTreeNode * nd = *ndi;
            std::cerr << (void*)nd << ": " << postorder_count << ": " << nd->get_index() << ": " << nd->get_label() << std::endl;
        }
    }
}
