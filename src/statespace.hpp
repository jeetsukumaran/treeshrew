#ifndef TREESHREW_STATESPACE_HPP
#define TREESHREW_STATESPACE_HPP

#include <iostream>
#include "character.hpp"
#include "genetree.hpp"

namespace treeshrew {

class StateSpace {

    public:
        StateSpace(unsigned long max_sequences, unsigned long max_sites);
        ~StateSpace();
        void initialize_with_tree_and_alignment(
                std::istream& tree_src,
                std::istream& alignment_src,
                const std::string& tree_format="newick",
                const std::string& alignment_format="fasta"
                );
        void dispose_gene_tree();
        void dispose_alignment();
        void write_phylogenetic_data(std::ostream&);

    private:
        NucleotideAlignment                 alignment_;
        GeneTree *                          gene_tree_;


}; // StateSpace

} // namespace treeshrew

#endif
