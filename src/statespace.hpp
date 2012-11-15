#ifndef TREESHREW_STATESPACE_HPP
#define TREESHREW_STATESPACE_HPP

#include "character.hpp"
#include "genetree.hpp"

namespace treeshrew {

class StateSpace {

    public:
        StateSpace(unsigned long max_sequences, unsigned long max_sites);
        ~StateSpace();
        void set_gene_tree(std::istream& src, const std::string& format="newick");
        void set_gene_tree(GeneTree * gene_tree);
        inline GeneTree * get_gene_tree() const {
            return this->gene_tree_;
        }
        void dispose_gene_tree();
        void set_alignment(std::istream& src, const std::string& format="fasta");
        void dispose_alignment();

    private:
        NucleotideAlignment                 alignment_;
        GeneTree *                          gene_tree_;


}; // StateSpace

} // namespace treeshrew

#endif
