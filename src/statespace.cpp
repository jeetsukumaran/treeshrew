#include "statespace.hpp"
#include "dataio.hpp"
#include "utility.hpp"

namespace treeshrew {

StateSpace::StateSpace(unsigned long max_sequences,
        unsigned long max_sites)
    : alignment_(max_sequences, max_sites)
    , gene_tree_(nullptr) {
}

StateSpace::~StateSpace() {
    this->dispose_gene_tree();
    this->dispose_alignment();
}

void StateSpace::load_short_reads(std::istream& src,
        const std::string& format) {
    NucleotideSequences dna;
    sequenceio::read_from_stream(dna, src, format);
    for (auto & seq : dna) {
        this->short_reads_.add(*seq);
    }
}

void StateSpace::initialize_with_tree_and_alignment(
        std::istream& tree_src,
        std::istream& alignment_src,
        const std::string& tree_format,
        const std::string& alignment_format
        ) {

    // gene tree
    std::vector<GeneTree *> trees;
    treeio::read_from_stream(trees,
            tree_src,
            tree_format,
            this->alignment_.get_max_sequences());
    if (trees.size() == 0) {
        treeshrew_abort("No trees found in data source");
    } else if (trees.size() > 1) {
        treeshrew_abort("Multiple trees found in data source");
    } else {
        this->gene_tree_ = trees[0];
    }
    this->gene_tree_->create_beagle_instance(this->alignment_.get_max_sites());

    // alignment
    NucleotideSequences dna;
    sequenceio::read_from_stream(dna, alignment_src, alignment_format);
    this->alignment_.set_num_active_sites(dna.get_num_sites());
    for (auto leaf_iter = this->gene_tree_->leaf_begin(); leaf_iter != this->gene_tree_->leaf_end(); ++leaf_iter) {
        NucleotideSequence * dseq = dna.get_sequence(leaf_iter->get_label());
        this->alignment_.new_sequence(&(*leaf_iter), dseq);
        this->gene_tree_->set_tip_partials(*leaf_iter,
                this->alignment_.get_partials_data(&(*leaf_iter)));
    }
}

void StateSpace::dispose_gene_tree() {
    if (this->gene_tree_) {
        this->gene_tree_->clear();
        delete this->gene_tree_;
        this->gene_tree_ = nullptr;
    }
}

void StateSpace::dispose_alignment() {
    this->alignment_.clear();
}

double StateSpace::calc_ln_probability_of_short_reads() {
    double ln_prob = 0.0;
    for (auto sri = this->short_reads_.cbegin(); sri != this->short_reads_.cend(); ++sri) {
        const ShortReadSequence& short_read = *sri;
        double sub_prob = 0.0;
        for (auto ndi = this->gene_tree_->leaf_begin(); ndi != this->gene_tree_->leaf_end(); ++ndi) {
            GeneNodeData& gnd = *ndi;
            sub_prob += short_read.calc_probability_of_sequence(
                    this->alignment_.sequence_states_cbegin(&gnd),
                    this->alignment_.sequence_states_cend(&gnd),
                    0.5);
        }
        std::cerr << "*** " << sub_prob << std::endl;
        ln_prob += std::log(sub_prob);
    }
    return ln_prob;
}

void StateSpace::write_phylogenetic_data(std::ostream& out) {
    out << "#NEXUS\n\n";

    // make sure branch lengths are stored when tree is read
    out << "begin Paup;\n    set storebr;\n    set warnreset = no warnroot = no;\nend;\n\n" << std::endl;

    // character matrix
    out << "begin data;\n";
    out << "    dimensions ntax=" << this->gene_tree_->get_num_leaves() << " nchar=" << this->alignment_.get_num_active_sites() << ";\n";
    out << "    format datatype=dna gap=- missing=? matchchar=.;\n";
    out << "    matrix\n";
    for (auto ndi = this->gene_tree_->leaf_begin(); ndi != this->gene_tree_->leaf_end(); ++ndi) {
        out << "    " << ndi->get_label() << "        ";
        this->alignment_.write_states_as_symbols(&(*ndi), out);
        out << "\n";
    }
    out << "    ;\nend;\n";

    // tree
    out << "begin trees;\n";
    out << "    tree 1 = ";
    treeio::write_newick(this->gene_tree_, out);
    out << "end;\n";

    // calculate the likelihood
    out << "begin paup;\n    set crit=likelihood;\n    lset userbr nst=1 rmatrix=estimate basefreq=equal rates=equal pinvar=0;\n    lscore;\nend;\n";
}

} // treeshrew
