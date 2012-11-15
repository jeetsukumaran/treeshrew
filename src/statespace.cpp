#include "statespace.hpp"
#include "dataio.hpp"
#include "utility.hpp"

namespace treeshrew {

StateSpace::StateSpace(unsigned long max_sequences,
        unsigned long max_sites)
    : alignment_(max_sequences, max_sites) {
}

StateSpace::~StateSpace() {
    this->dispose_gene_tree();
    this->dispose_alignment();
}

void StateSpace::set_gene_tree(std::istream& src, const std::string& format) {
    std::vector<GeneTree *> trees;
    treeio::read_from_stream(trees,
            src,
            format,
            this->alignment_.get_max_sequences());
    if (trees.size() == 0) {
        treeshrew_abort("No trees found in data source");
    } else if (trees.size() > 1) {
        treeshrew_abort("Multiple trees found in data source");
    } else {
        this->set_gene_tree(trees[0]);
    }

}

void StateSpace::set_gene_tree(GeneTree * gene_tree) {
    this->dispose_gene_tree();
    this->gene_tree_ = gene_tree;
}

void StateSpace::dispose_gene_tree() {
    if (this->gene_tree_) {
        this->gene_tree_->clear();
        delete this->gene_tree_;
        this->gene_tree_ = nullptr;
    }
}

void StateSpace::set_alignment(std::istream& src, const std::string& format) {
    this->dispose_alignment();
    NucleotideSequences dna;
    sequenceio::read_from_stream(dna, src, format);
    this->alignment_.set_num_active_sites(dna.get_num_sites());
    for (auto leaf_iter = this->gene_tree_->leaf_begin(); leaf_iter != this->gene_tree_->leaf_end(); ++leaf_iter) {
        NucleotideSequence * dseq = dna.get_sequence(leaf_iter->get_label());
        NucleotideSequence * aseq = this->alignment_.new_sequence(&(*leaf_iter), dseq);
    }
}

void StateSpace::dispose_alignment() {
    this->alignment_.clear();
}

} // treeshrew
