#include <vector>
#include <iostream>
#include <string>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include "../../src/character.hpp"
#include "../../src/dataio.hpp"
#include "../../src/utility.hpp"

int main(int argc, char * argv[]) {
    if (argc == 1) {
        std::cerr << "Usage: read_tree <DATAFILE> [nexus|phylip|fasta]" << std::endl;
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
        format = "fasta";
    }
    std::string filepath(argv[1]);
    // treeshrew::NucleotideSequences * dna = treeshrew::create_sequences_from_filepath<treeshrew::NucleotideSequences>(filepath, format);
    // dna->write_fasta(std::cout);
    treeshrew::NucleotideSequences dna;
    treeshrew::sequenceio::read_from_filepath(dna, filepath, format);
    for (auto & seq : dna) {
        std::cout << seq->get_label() << ":";
        std::copy(seq->partials_cbegin(), seq->partials_cend(), std::ostream_iterator<double>(std::cout, ";"));
        std::cout << std::endl;
    }
}
