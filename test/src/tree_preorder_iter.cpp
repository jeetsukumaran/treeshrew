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
        int preorder_count = 0;
        for (auto ndi = tree->preorder_begin(); ndi != tree->preorder_end(); ++ndi, ++preorder_count) {
            std::cerr << (void*)ndi.node() << ": " << preorder_count << ":" << ndi->get_index() << ": " << ndi->get_label() << std::endl;
        }
    }
}


