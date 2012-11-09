#include <config.h>
#include <libhmsbeagle/beagle.h>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <array>

static const unsigned int NUM_STATES = 4;

class PhylogeneticNode {
    public:
        PhylogeneticNode(PhylogeneticNode * parent=nullptr,
                double edge_length=0.0,
                int index=0) :
            parent_(parent),
            edge_length_(edge_length),
            index_(index) {
        }
        ~PhylogeneticNode() {
            for (auto & ch : this->child_nodes_) {
                delete ch;
            }
        }
        PhylogeneticNode * new_child(double edge_length, int index=0) {
            PhylogeneticNode * ch = new PhylogeneticNode(this, edge_length, index);
            this->child_nodes_.push_back(ch);
            return ch;
        }
        std::vector<PhylogeneticNode *>& get_leaf_nodes(std::vector<PhylogeneticNode *>& nodes) {
            if (this->child_nodes_.begin() == this->child_nodes_.end()) {
                nodes.push_back(this);
            } else {
                for (auto & ni : this->child_nodes_) {
                    ni->get_leaf_nodes(nodes);
                }
            }
            return nodes;
        }
        std::vector<PhylogeneticNode *>& get_nodes_postorder(std::vector<PhylogeneticNode *>& nodes) {
            for (auto & ni : this->child_nodes_) {
                ni->get_nodes_postorder(nodes);
            }
            nodes.push_back(this);
            return nodes;
        }
        std::vector<PhylogeneticNode *>& get_internal_nodes_postorder(std::vector<PhylogeneticNode *>& nodes, bool include_root=false) {
            if (this->is_leaf_node() || (!include_root && this->parent_ == nullptr)) {
                return nodes;
            }
            for (auto & ni : this->child_nodes_) {
                ni->get_internal_nodes_postorder(nodes);
            }
            nodes.push_back(this);
            return nodes;
        }
        int get_child_node_index(unsigned long child_index) {
            return this->child_nodes_.at(child_index)->index_;
        }
        void set_state_vector_by_symbol(const std::string& seq_str) {
            this->state_vector_.clear();
            for (auto & ch_str : seq_str) {
                if (ch_str == 'A' || ch_str == 'a') {
                    this->state_vector_.push_back(0);
                    this->partials_.push_back(1.0);
                    this->partials_.push_back(0.0);
                    this->partials_.push_back(0.0);
                    this->partials_.push_back(0.0);
                } else if (ch_str == 'C' || ch_str == 'c') {
                    this->state_vector_.push_back(1);
                    this->partials_.push_back(0.0);
                    this->partials_.push_back(1.0);
                    this->partials_.push_back(0.0);
                    this->partials_.push_back(0.0);
                } else if (ch_str == 'G' || ch_str == 'g') {
                    this->state_vector_.push_back(2);
                    this->partials_.push_back(0.0);
                    this->partials_.push_back(0.0);
                    this->partials_.push_back(1.0);
                    this->partials_.push_back(0.0);
                } else if (ch_str == 'T' || ch_str == 't') {
                    this->state_vector_.push_back(3);
                    this->partials_.push_back(0.0);
                    this->partials_.push_back(0.0);
                    this->partials_.push_back(0.0);
                    this->partials_.push_back(1.0);
                } else {
                    this->state_vector_.push_back(4);
                    this->partials_.push_back(1.0);
                    this->partials_.push_back(1.0);
                    this->partials_.push_back(1.0);
                    this->partials_.push_back(1.0);
                }
            }
        }
        void set_label(const std::string& label) {
            this->label_ = label;
        }
        const std::string& get_label() {
            return this->label_;
        }
        int get_index() {
            return this->index_;
        }
        void set_index(int index) {
            this->index_ = index;
        }
        double get_edge_len() {
            return this->edge_length_;
        }
        unsigned long get_state_vector_len() {
            return this->state_vector_.size();
        }
        int get_state(unsigned long index) {
            return this->state_vector_[index];
        }
        int * get_state_vector_data() {
            return this->state_vector_.data();
        }
        bool is_leaf_node() {
            return this->child_nodes_.size() == 0;
        }
        double * get_partials_data() {
            return this->partials_.data();
        }
        double get_partial(unsigned long index) {
            return this->partials_[index];
        }
        void write_newick(std::ostream& out) {
            if (this->child_nodes_.begin() == this->child_nodes_.end()) {
                out << this->label_;
            } else {
                out << "(";
                for (std::vector<PhylogeneticNode *>::const_iterator ni = this->child_nodes_.begin();
                        ni != this->child_nodes_.end();
                        ++ni) {
                    if (ni != this->child_nodes_.begin()) {
                        out << ", ";
                    }
                    (*ni)->write_newick(out);
                }
                out << ")";
            }
            out << ":" << this->edge_length_;
        }
        void write_states_as_symbols(std::ostream& out) {
            for (auto &s : this->state_vector_) {
                if (s == 0) {
                    out << "A";
                } else if (s == 1) {
                    out << "C";
                } else if (s == 2) {
                    out << "G";
                } else if (s == 3) {
                    out << "T";
                } else {
                    out << "-";
                }
            }
        }
        void write_data(std::ostream& out) {
            out << this->label_ << "        ";
            this->write_states_as_symbols(out);
        }
    private:
        PhylogeneticNode *              parent_;
        double                          edge_length_;
        int                             index_;
        std::string                     label_;
        std::vector<int>                state_vector_;
        std::vector<double>             partials_;
        std::vector<PhylogeneticNode *>   child_nodes_;

}; // PhylogeneticNode

class PhylogeneticTree {
    public:
        PhylogeneticTree() {
            this->seed_node_ = new PhylogeneticNode(nullptr);
        }
        ~PhylogeneticTree() {
            delete this->seed_node_;
        }
        PhylogeneticNode * get_seed_node() {
            return this->seed_node_;
        }
        void read_data() {
            this->seed_node_->set_index(6);
            PhylogeneticNode * i1 = this->seed_node_->new_child(0.02, 5);
            i1->new_child(0.01, 0);
            i1->new_child(0.02, 1);
            PhylogeneticNode * i2 = this->seed_node_->new_child(0.05, 4);
            i2->new_child(0.03, 2);
            i2->new_child(0.04, 3);
            std::vector<PhylogeneticNode *> leaves;
            this->get_leaf_nodes(leaves);
            std::string label;
            std::string sequence;
            unsigned long index = 0;
            std::ifstream infile("/Users/jeet/Sandbox/prototypes/treeshrew/treeshrew3/test/data/d20121031.txt");
            while (infile >> label >> sequence) {
                std::cerr << "Reading " << index+1 << "/" << leaves.size() << ": " << label << std::endl;
                leaves.at(index)->set_state_vector_by_symbol(sequence);
                leaves.at(index)->set_label(label);
                ++index;
            }
            for (auto & nd : leaves) {
                std::cerr << "Tip '" << nd->get_label() << "': " << nd->get_state_vector_len() << std::endl;
            }
        }
        std::vector<PhylogeneticNode *>& get_leaf_nodes(std::vector<PhylogeneticNode *>& nodes) {
            return this->seed_node_->get_leaf_nodes(nodes);
        }

        double calc_ln_likelihood() {

            //////////////////////////////////////////////////////////////////////////////
            // Setup and initialze Beagle

            std::vector<PhylogeneticNode *> nodes;
            this->seed_node_->get_nodes_postorder(nodes);
            std::vector<PhylogeneticNode *> leaf_nodes;
            this->seed_node_->get_leaf_nodes(leaf_nodes);
            std::vector<PhylogeneticNode *> internal_nodes;
            this->seed_node_->get_internal_nodes_postorder(internal_nodes, true);
            int num_sites = leaf_nodes[0]->get_state_vector_len();
            int num_nodes = nodes.size();
            int num_tip_nodes = leaf_nodes.size();
            int num_int_nodes = internal_nodes.size() ;
            int num_partials = num_nodes;
            BeagleInstanceDetails * return_info = new BeagleInstanceDetails();
            int beagle_instance = beagleCreateInstance(
                num_tip_nodes,  // Number of tip data elements (input)
                num_partials,   // Number of partials buffers to create (input) -- internal node count
                num_tip_nodes,  // Number of compact state representation buffers to create -- for use with setTipStates (input)
                4,              // Number of states in the continuous-time Markov chain (input) -- DNA
                num_sites,      // Number of site patterns to be handled by the instance (input) -- not compressed in this case
                1,              // Number of eigen-decomposition buffers to allocate (input)
                num_nodes,      // Number of transition matrix buffers (input) -- one per edge
                1,              // Number of rate categories
                0,              // Number of scaling buffers -- can be zero if scaling is not needed
                NULL,           // List of potential resource on which this instance is allowed (input, NULL implies no restriction
                0,              // Length of resourceList list (input) -- not needed to use the default hardware config
                0,              // Bit-flags indicating preferred implementation charactertistics, see BeagleFlags (input)
                0,              // Bit-flags indicating required implementation characteristics, see BeagleFlags (input)
                return_info
                );
            if (beagle_instance < 0) {
                std::cerr << "\n\n***ERROR*** Failed to obtain BEAGLE instance\n" << std::endl;
                exit(1);
            }
            int ret_code = 0;
            for (int i; i < num_tip_nodes; ++i) {
                // ret_code = beagleSetTipPartials(
                //             beagle_instance,            // instance
                //             leaf_nodes[i]->get_index(),                          // bufferIndex
                //             leaf_nodes[i]->get_partials_data()            // inPartials
                // );
                ret_code = beagleSetTipStates(
                        beagle_instance,
                        leaf_nodes[i]->get_index(),
                        leaf_nodes[i]->get_state_vector_data()
                        );
                if (ret_code != 0) {
                    std::cerr << "\n\n***ERROR*** Failed to set tip data\n" << std::endl;
                    exit(1);
                }
            }

            // let all sites have equal weight
            std::vector<double> pattern_weights( num_sites, 1 );
            beagleSetPatternWeights(beagle_instance, pattern_weights.data());

            // create array of state background frequencies
            double freqs[4] = { 0.25, 0.25, 0.25, 0.25 };
            beagleSetStateFrequencies(beagle_instance, 0, freqs);

            // create an array containing site category weights and rates
            const double weights[1] = { 1.0 };
            const double rates[1] = { 1.0 };
            beagleSetCategoryWeights(beagle_instance, 0, weights);
            beagleSetCategoryRates(beagle_instance, rates);


            double evec[4 * 4] = {
                1.0,  2.0,  0.0,  0.5,
                1.0,  -2.0,  0.5,  0.0,
                1.0,  2.0, 0.0,  -0.5,
                1.0,  -2.0,  -0.5,  0.0
            };

            // JC69 model inverse eigenvector matrix
            double ivec[4 * 4] = {
                0.25,  0.25,  0.25,  0.25,
                0.125,  -0.125,  0.125,  -0.125,
                0.0,  1.0,  0.0,  -1.0,
                1.0,  0.0,  -1.0,  0.0
            };

            // JC69 model eigenvalues
            double eval[4] = { 0.0, -1.3333333333333333, -1.3333333333333333, -1.3333333333333333 };

            ret_code = beagleSetEigenDecomposition(
                    beagle_instance,                 // instance
                    0,                               // eigenIndex,
                    (const double *)evec,            // inEigenVectors,
                    (const double *)ivec,            // inInverseEigenVectors,
                    eval);                           // inEigenValues

            if (ret_code != 0) {
                std::cerr << "\n\n***ERROR*** Failed to set eigen decomposition\n" << std::endl;
                exit(1);
            }

            //////////////////////////////////////////////////////////////////////////////
            // Calculate log-likelihood

            // a list of indices and edge lengths
            // these get used to tell beagle which edge length goes with which node
            std::vector<int> node_indices;
            std::vector<double> edge_lens;
            for (auto &nd : nodes) {
                node_indices.push_back(nd->get_index());
                edge_lens.push_back(nd->get_edge_len());
            }
            // tell BEAGLE to populate the transition matrices for the above edge lengthss
            beagleUpdateTransitionMatrices(beagle_instance,     // instance
                    0,             // eigenIndex
                    node_indices.data(),   // probabilityIndices
                    NULL,          // firstDerivativeIndices
                    NULL,          // secondDervativeIndices
                    edge_lens.data(),   // edgeLengths
                    node_indices.size());            // count

            // create a list of partial likelihood update operations
            // the order is [dest, sourceScaling, destScaling, source1, matrix1, source2, matrix2]
            // these operations say: first peel node 0 and 1 to calculate the per-site partial likelihoods, and store them
            // in buffer 3.  Then peel node 2 and buffer 3 and store the per-site partial likelihoods in buffer 4.
            // BeagleOperation operations[2] = {
            //     {3, BEAGLE_OP_NONE, BEAGLE_OP_NONE, 0, 0, 1, 1},
            //     {4, BEAGLE_OP_NONE, BEAGLE_OP_NONE, 2, 2, 3, 3}
            // };
            std::vector<BeagleOperation> beagle_operations;
            int ch1_idx = 0;
            int ch2_idx = 0;
            for (auto &nd : internal_nodes) {
                ch1_idx = nd->get_child_node_index(0);
                ch2_idx = nd->get_child_node_index(1);
                std::cerr << nd->get_index() << ": " << ch1_idx << ", " << ch2_idx << std::endl;
                beagle_operations.push_back(
                        {nd->get_index(), BEAGLE_OP_NONE, BEAGLE_OP_NONE, ch1_idx, ch1_idx, ch2_idx, ch2_idx}
                        );
            }

            // this invokes all the math to carry out the likelihood calculation
            beagleUpdatePartials( beagle_instance,      // instance
                    beagle_operations.data(),     // eigenIndex
                    beagle_operations.size(),              // operationCount
                    BEAGLE_OP_NONE);             // cumulative scale index

            // for (auto &nd : nodes) {
            //     std::cerr << nd->get_index() << ":";
            //     for (unsigned int i = 0; i < 4; ++i) {
            //         std::cerr << "     " << nd->get_partial(i);
            //     }
            //     std::cerr << std::endl;
            // }

            double logL = 0;
            int root_index[1] = {this->seed_node_->get_index()};
            int category_weight_index[1] = {0};
            int state_freq_index[1] = {0};
            int cumulative_scale_index[1] = {BEAGLE_OP_NONE};

            // calculate the site likelihoods at the root node
            // this integrates the per-site root partial likelihoods across sites, background state frequencies, and rate categories
            // results in a single log likelihood, output here into logL
            beagleCalculateRootLogLikelihoods(beagle_instance,               // instance
                    root_index,// bufferIndices
                    category_weight_index,                // weights
                    state_freq_index,                 // stateFrequencies
                    cumulative_scale_index,     // scaleBuffer to use
                    1,                      // count
                    &logL);         // outLogLikelihoods

            return logL;
        }

        void write_newick(std::ostream& out) {
            out << "[&R] ";
            this->seed_node_->write_newick(out);
            out << ";";
        }

        void write_sequences(std::ostream& out) {
            std::vector<PhylogeneticNode *> leaf_nodes;
            this->seed_node_->get_leaf_nodes(leaf_nodes);
            for (auto & nd : leaf_nodes) {
                nd->write_data(out);
                out << std::endl;
            }
            out << ";" << std::endl;
        }

    private:
        PhylogeneticNode *                  seed_node_;
}; // PhylogeneticTree


int main() {
    PhylogeneticTree tree;
    tree.read_data();
    double lnlike = tree.calc_ln_likelihood();
    std::ostream& out = std::cout;
    out << "#NEXUS\n";
    out << "Begin Paup;\n    set storebr;\n    set warnreset = no warnroot = no; End;" << std::endl;
    out << "begin data;\n";
    out << "    dimensions ntax=4 nchar=1314;\n";
    out << "    format datatype=dna gap=- missing=? matchchar=.;\n";
    out << "    matrix\n";
    tree.write_sequences(out);
    out << ";\n";
    out << "end;\n";
    out << "begin trees;\n";
    out << "    tree 1 = ";
    tree.write_newick(out);
    out << "\n;end;\n";
    out << "begin paup;\n";
    out << "    set crit=likelihood;\n";
    out << "    lset userbr nst=1 rmatrix=estimate basefreq=equal rates=equal pinvar=0;\n";
    out << "    lscore;\n";
    out << "end;\n";
    std::cout << "[! beagle = " << lnlike << "]" << std::endl;
}
