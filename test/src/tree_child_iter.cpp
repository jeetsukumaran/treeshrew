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
        for (auto ndi = tree->postorder_begin(); ndi != tree->postorder_end(); ++ndi, ++postorder_count) {
            std::cout << ndi->get_label();
            for (auto chi = tree->children_begin(ndi);
                    chi != tree->children_end(ndi);
                    ++chi) {
                std::cout << "\t" << chi->get_label();
            }
            std::cout << std::endl;
        }
    }
}
