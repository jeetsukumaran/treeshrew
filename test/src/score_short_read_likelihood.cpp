#include <vector>
#include <iomanip>
#include <iostream>
#include <string>
#include "../src/statespace.hpp"


int main(int argc, char * argv[]) {
    // treeshrew::GeneTree g(10);
    if (argc < 4) {
        std::cerr << "Usage: read_tree <SHORT-READ-DATAFILE> <ALIGNMENT-DATAFILE> <TREE-FILE>" << std::endl;
        exit(1);
    }
    std::ifstream short_read_src(argv[1]);
    std::ifstream data_src(argv[2]);
    std::ifstream tree_src(argv[3]);
    treeshrew::StateSpace state_space(100, 50000);
    state_space.load_short_reads(short_read_src);
    state_space.initialize_with_tree_and_alignment(tree_src, data_src);
    state_space.write_phylogenetic_data(std::cerr);
    std::cout << "Probability of Alignment | Short Reads: " << std::setprecision(12) << state_space.calc_ln_probability_of_short_reads() << std::endl;
    std::cout << "  Probability of Gene Tree | Alignment: " << std::setprecision(12) << state_space.get_gene_tree()->calc_ln_probability() << std::endl;
}
