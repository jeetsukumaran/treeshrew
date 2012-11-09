#include <vector>
#include <iomanip>
#include <iostream>
#include <string>
#include "../../src/dataio.hpp"
#include "../../src/character.hpp"
#include "../../src/genetree.hpp"
#include "../../src/utility.hpp"

void write_paup_calculator(treeshrew::GeneTree * tree, treeshrew::NucleotideSequences& data, std::ostream& out) {
    unsigned long num_seqs = data.get_num_sequences();
    out << "#NEXUS\n";
    out << "Begin Paup;\n    set storebr;\n    set warnreset = no warnroot = no; End;" << std::endl;
    out << "begin data;\n";
    out << "    dimensions ntax=" << num_seqs << " nchar=" << data.get_num_sites() << ";\n";
    out << "    format datatype=dna gap=- missing=? matchchar=.;\n";
    out << "    matrix\n";
    for (unsigned long idx = 0; idx < num_seqs; ++idx) {
        treeshrew::NucleotideSequence * seq = data.get_sequence(idx);
        out << "        " << seq->get_label() << "    ";
        seq->write_states_as_symbols(out);
        out << "\n";
    }
    out << "    ;\nend;\n\n";
    out << "begin trees;\n";
    out << "    tree 1 = ";
    treeshrew::treeio::write_newick(tree, out);
    out << "end;\n";
    out << "begin paup;\n";
    out << "    set crit=likelihood;\n";
    out << "    lset userbr nst=1 rmatrix=estimate basefreq=equal rates=equal pinvar=0;\n";
    out << "    lscore;\n";
    out << "end;\n";
}

int main(int argc, char * argv[]) {
    // treeshrew::GeneTree g(10);
    if (argc < 3) {
        std::cerr << "Usage: read_tree <NEWICK-TREEFILE> <FASTA-DATAFILE>" << std::endl;
        exit(1);
    }
    std::string tree_filepath(argv[1]);
    std::string data_filepath(argv[2]);
    treeshrew::NucleotideSequences data;
    treeshrew::sequenceio::read_from_filepath(data, data_filepath, "fasta");
    std::vector<treeshrew::GeneTree *> trees;
    treeshrew::treeio::read_from_filepath(trees, tree_filepath, "newick");
    for (auto & tree : trees) {
        tree->create_beagle_instance(data.get_num_sites());
        data.set_tip_data(tree);
        write_paup_calculator(tree, data, std::cerr);
        std::cout << std::setprecision(12) << tree->calc_ln_probability() << std::endl;
    }
}
