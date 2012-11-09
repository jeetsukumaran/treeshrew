#include "testutils.hpp"

std::vector<treeshrew::GeneTree *> get_trees(int argc, char * argv[]) {
    if (argc == 1) {
        std::cerr << "Usage: " << argv[0] << " <TREEFILE> [nexus|newick]" << std::endl;
        exit(1);
    }
    if (argc > 3) {
        std::cerr << "Expecting at most two arguments" << std::endl;
        exit(1);
    }
    std::string format;
    if (argc == 3) {
        format = argv[2];
    } else {
        format = "nexus";
    }
    std::string filepath(argv[1]);
    std::vector<treeshrew::GeneTree *> trees;
    treeshrew::treeio::read_from_filepath(trees, filepath, format);
    return trees;
}
