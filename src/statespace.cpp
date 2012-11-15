#include "statespace.hpp"
#include "dataio.hpp"
#include "utility.hpp"

namespace treeshrew {

StateSpace::StateSpace(unsigned long max_sequences,
        unsigned long max_sites)
    : alignment_(max_sequences, max_sites) {
}

StateSpace::~StateSpace() {
}

void StateSpace::dispose_gene_tree() {
    if (this->gene_tree_) {
        this->gene_tree_->clear();
        delete this->gene_tree_;
        this->gene_tree_ = nullptr;
    }
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
    this->alignment_.set_gene_tree(this->gene_tree_);
}

} // treeshrew
