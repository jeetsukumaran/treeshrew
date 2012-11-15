#include <vector>
#include <iomanip>
#include <iostream>
#include <string>
#include "../src/statespace.hpp"


int main(int argc, char * argv[]) {
    // treeshrew::GeneTree g(10);
    if (argc < 3) {
        std::cerr << "Usage: read_tree <NEWICK-TREEFILE> <FASTA-DATAFILE>" << std::endl;
        exit(1);
    }
    std::ifstream tree_src(argv[1]);
    std::ifstream data_src(argv[2]);
    treeshrew::StateSpace state_space(100, 50000);
    state_space.initialize_with_tree_and_alignment(tree_src, data_src);
    state_space.write_phylogenetic_data(std::cerr);
    std::cout << std::setprecision(12) << state_space.get_gene_tree()->calc_ln_probability() << std::endl;
}
