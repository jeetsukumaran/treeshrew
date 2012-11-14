#include <iostream>
#include "utility.hpp"
#include "genetree.hpp"

namespace treeshrew {

////////////////////////////////////////////////////////////////////////////////
// GeneTreeNode

GeneTree::GeneTree(unsigned long max_tips):
        Tree<GeneNodeData>(false),
        max_tips_(max_tips),
        leaf_node_allocator_(max_tips * 2, 0),
        internal_node_allocator_(max_tips * 2 + 1, max_tips * 2),
        beagle_instance_(-1),
        beagle_return_info_(nullptr) {
    this->create(this->allocate_internal_node(),
        this->allocate_internal_node());
}

GeneTree::~GeneTree() {
    this->free_beagle_instance();
    this->clear();
}

void GeneTree::clear() {
}

int GeneTree::create_beagle_instance(int num_sites) {
    this->free_beagle_instance();
    this->beagle_return_info_ = new BeagleInstanceDetails();
    int num_tip_nodes = this->max_tips_ * 2;
    int num_internal_nodes = num_tip_nodes + 1;
    int total_nodes = num_tip_nodes + num_internal_nodes;
    this->beagle_instance_ = beagleCreateInstance(
        num_tip_nodes,          // Number of tip data elements (input)
        num_internal_nodes,     // Number of partials buffers to create (input) -- internal node count
        num_tip_nodes,          // Number of compact state representation buffers to create -- for use with setTipStates (input)
        4,                      // Number of states in the continuous-time Markov chain (input) -- DNA
        num_sites,              // Number of site patterns to be handled by the instance (input) -- not compressed in this case
        1,                      // Number of eigen-decomposition buffers to allocate (input)
        total_nodes,            // Number of transition matrix buffers (input) -- one per edge
        1,                      // Number of rate categories
        0,                      // Number of scaling buffers -- can be zero if scaling is not needed
        NULL,                   // List of potential resource on which this instance is allowed (input, NULL implies no restriction
        0,                      // Length of resourceList list (input) -- not needed to use the default hardware config
        0,                      // Bit-flags indicating preferred implementation charactertistics, see BeagleFlags (input)
        0,                      // Bit-flags indicating required implementation characteristics, see BeagleFlags (input)
        this->beagle_return_info_
        );
    if (this->beagle_instance_ < 0) {
        treeshrew_abort("Failed to obtain BEAGLE instance");
    }

    // let all sites have equal weight
    std::vector<double> pattern_weights( num_sites, 1 );
    beagleSetPatternWeights(this->beagle_instance_, pattern_weights.data());

    // create array of state background frequencies
    double freqs[4] = { 0.25, 0.25, 0.25, 0.25 };
    beagleSetStateFrequencies(this->beagle_instance_, 0, freqs);

    // create an array containing site category weights and rates
    const double weights[1] = { 1.0 };
    const double rates[1] = { 1.0 };
    beagleSetCategoryWeights(this->beagle_instance_, 0, weights);
    beagleSetCategoryRates(this->beagle_instance_, rates);


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

    int ret_code = beagleSetEigenDecomposition(
            this->beagle_instance_,                 // instance
            0,                               // eigenIndex,
            (const double *)evec,            // inEigenVectors,
            (const double *)ivec,            // inInverseEigenVectors,
            eval);                           // inEigenValues

    if (ret_code != 0) {
        treeshrew_abort("Failed to set eigen decomposition");
    }

    return this->beagle_instance_;
}

int GeneTree::set_tip_states(const GeneNodeData& tip, const int * data) {
    int beagle_index = tip.get_index() ;
    int ret_code = beagleSetTipStates(
            this->beagle_instance_,
            beagle_index,
            data
            );
    if (ret_code != 0) {
        treeshrew_abort("Failed to set tip states for node index ", beagle_index);
    }
    return ret_code;
}

int GeneTree::set_tip_partials(const GeneNodeData& tip, const double * data) {
    int beagle_index = tip.get_index() ;
    int ret_code = beagleSetTipPartials(
            this->beagle_instance_,
            beagle_index,
            data
            );
    if (ret_code != 0) {
        treeshrew_abort("Failed to set tip partials for node index ", beagle_index);
    }
    return ret_code;
}

double GeneTree::calc_ln_probability() {
    std::vector<int> node_indices;
    std::vector<double> edge_lens;
    for (GeneTree::postorder_iterator ndi = this->postorder_begin(); ndi != this->postorder_end(); ++ndi) {
        node_indices.push_back(ndi->get_index());
        edge_lens.push_back(ndi->get_edge_length());
    }
    // tell BEAGLE to populate the transition matrices for the above edge lengthss
    beagleUpdateTransitionMatrices(this->beagle_instance_,     // instance
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
    for (GeneTree::postorder_iterator ndi = this->postorder_begin(); ndi != this->postorder_end(); ++ndi) {
        if (ndi.is_leaf()) {
            continue;
        }
        ch1_idx = ndi.first_child().get_index();
        ch2_idx = ndi.last_child().get_index();
        // std::cerr << nd->get_index() << ": " << ch1_idx << ", " << ch2_idx << std::endl;
        beagle_operations.push_back(
                {ndi->get_index(), BEAGLE_OP_NONE, BEAGLE_OP_NONE, ch1_idx, ch1_idx, ch2_idx, ch2_idx}
                );
    }

    // this invokes all the math to carry out the likelihood calculation
    beagleUpdatePartials( this->beagle_instance_,      // instance
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
    int root_index[1] = {this->head_node_->data().get_index()};
    int category_weight_index[1] = {0};
    int state_freq_index[1] = {0};
    int cumulative_scale_index[1] = {BEAGLE_OP_NONE};

    // calculate the site likelihoods at the root node
    // this integrates the per-site root partial likelihoods across sites, background state frequencies, and rate categories
    // results in a single log likelihood, output here into logL
    beagleCalculateRootLogLikelihoods(this->beagle_instance_,               // instance
            root_index,// bufferIndices
            category_weight_index,                // weights
            state_freq_index,                 // stateFrequencies
            cumulative_scale_index,     // scaleBuffer to use
            1,                      // count
            &logL);         // outLogLikelihoods

    return logL;
    return 0.0;
}

void GeneTree::free_beagle_instance() {
    if (this->beagle_return_info_) {
        delete this->beagle_return_info_;
    }
    this->beagle_return_info_ = nullptr;
    if (this->beagle_instance_ >= 0) {
        beagleFinalizeInstance(this->beagle_instance_);
    }
    this->beagle_instance_ = -1;
}

} // treeshrew
