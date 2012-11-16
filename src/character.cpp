#include <algorithm>
#include <iterator>
#include "character.hpp"

namespace treeshrew {

//////////////////////////////////////////////////////////////////////////////
// Utility Functions

unsigned long hamming_distance(const CharacterStateVectorType::const_iterator& short_read_begin,
        const CharacterStateVectorType::const_iterator& short_read_end,
        const CharacterStateVectorType::const_iterator& long_read_begin) {
    return std::inner_product(
            short_read_begin, short_read_end, long_read_begin,
            0, std::plus<unsigned int>(),
            std::not2(std::equal_to<CharacterStateVectorType::value_type>()));
}

unsigned long sliding_hamming_distance(const CharacterStateVectorType& short_read,
        const CharacterStateVectorType& long_read) {
    CharacterStateVectorType::const_iterator start_pos = long_read.begin();
    CharacterStateVectorType::const_iterator stop_pos = long_read.end() - short_read.size() + 1;
    assert(stop_pos >= start_pos);
    unsigned long d = 0;
    while (start_pos < stop_pos) {
        d += hamming_distance(short_read.begin(), short_read.end(), start_pos);
        ++start_pos;
    }
    return d;
}


//////////////////////////////////////////////////////////////////////////////
// NucleotideSequence

const std::map<char, CharacterStateType> NucleotideSequence::symbol_to_state_map_ {
    {'A', 0},
    {'C', 1},
    {'G', 2},
    {'T', 3},
    {'U', 3},
    {'N', 4},
    {'X', 4},
    {'-', 4},
    {'?', 4},
    {'R', 5},
    {'Y', 6},
    {'M', 7},
    {'W', 8},
    {'S', 9},
    {'K', 10},
    {'V', 11},
    {'H', 12},
    {'D', 13},
    {'B', 14}
};

const std::map<CharacterStateType, char> NucleotideSequence::state_to_symbol_map_ {
    { 0, 'A'},
    { 1, 'C'},
    { 2, 'G'},
    { 3, 'T'},
    // { 4, 'N'},
    // { 4, 'X'},
    { 4, '-'},
    // { 4, '?'},
    { 5, 'R'},
    { 6, 'Y'},
    { 7, 'M'},
    { 8, 'W'},
    { 9, 'S'},
    {10, 'K'},
    {12, 'H'},
    {13, 'D'},
    {14, 'B'},
    {11, 'V'}
};

const std::map<CharacterStateType, std::array<double, 4>> NucleotideSequence::state_to_partials_map_ {
    { 0, {{1.0, 0.0, 0.0, 0.0}}},   // A
    { 1, {{0.0, 1.0, 0.0, 0.0}}},   // C
    { 2, {{0.0, 0.0, 1.0, 0.0}}},   // G
    { 3, {{0.0, 0.0, 0.0, 1.0}}},   // T, U
    { 4, {{1.0, 1.0, 1.0, 1.0}}},   // N, X, -, ?
    { 5, {{1.0, 0.0, 1.0, 0.0}}},   // R
    { 6, {{0.0, 1.0, 0.0, 1.0}}},   // Y
    { 7, {{1.0, 1.0, 0.0, 0.0}}},   // M
    { 8, {{1.0, 0.0, 0.0, 1.0}}},   // W
    { 9, {{0.0, 1.0, 1.0, 0.0}}},   // S
    {10, {{0.0, 0.0, 1.0, 1.0}}},   // K
    {11, {{1.0, 1.0, 1.0, 0.0}}},   // V
    {12, {{1.0, 1.0, 0.0, 1.0}}},   // H
    {13, {{1.0, 0.0, 1.0, 1.0}}},   // D
    {14, {{0.0, 1.0, 1.0, 1.0}}}    // B
};

const CharacterStateType NucleotideSequence::missing_data_state(NucleotideSequence::symbol_to_state_map_.find('?')->second);
const std::array<double, 4> NucleotideSequence::missing_data_partials = NucleotideSequence::state_to_partials_map_.find(NucleotideSequence::missing_data_state)->second;

NucleotideSequence::NucleotideSequence() {
}

NucleotideSequence::NucleotideSequence(const std::string& label) :
    label_(label) {
}

NucleotideSequence::~NucleotideSequence() {
}

void NucleotideSequence::write_states_as_symbols(std::ostream& out) const {
    for (auto & s : this->sequence_) {
        out << this->state_to_symbol_map_.find(s)->second;
    }
}

void NucleotideSequence::write_states_as_symbols(std::ostream& out,
    const CharacterStateVectorType::const_iterator& begin,
    const CharacterStateVectorType::const_iterator& end) const {
    for (auto iter = begin; iter != end; ++iter) {
        out << this->state_to_symbol_map_.find(*iter)->second;
    }
}

std::string NucleotideSequence::get_states_as_symbols() const {
    std::ostringstream out;
    this->write_states_as_symbols(out);
    return out.str();
}

//////////////////////////////////////////////////////////////////////////////
// NucleotideSequences

NucleotideSequences::NucleotideSequences() {
}

NucleotideSequences::~NucleotideSequences() {
    this->clear();
}

void NucleotideSequences::clear() {
    for (auto & v : this->sequences_) {
        if (v) {
            delete v;
        }
    }
    this->sequences_.clear();
    this->label_sequence_map_.clear();
}

void NucleotideSequences::set_tip_data(GeneTree * gene_tree) {
    unsigned long idx=0;
    for (auto leaf_iter = gene_tree->leaf_begin(); leaf_iter != gene_tree->leaf_end(); ++leaf_iter, ++idx) {
        const std::string& label = leaf_iter->get_label();
        NucleotideSequence * seq = this->label_sequence_map_[label];
        if (!seq) {
            treeshrew_abort("Null sequence for taxon '", label, "'");
        }
        gene_tree->set_tip_partials(*leaf_iter, seq->partials_data());
    }
}

//////////////////////////////////////////////////////////////////////////////
// NucleotideAlignment

NucleotideAlignment::NucleotideAlignment(unsigned long max_sequences,
        unsigned long max_sites)
        : max_sequences_(max_sequences)
        , max_sites_(max_sites)
        , num_active_sites_(0) {
    this->create();
}

NucleotideAlignment::~NucleotideAlignment() {
    this->clear();
}

void NucleotideAlignment::create() {
    this->sequence_storage_.reserve(this->max_sequences_);
    for (unsigned long row = 0; row < this->max_sequences_; ++row) {
        NucleotideSequence * s = new NucleotideSequence();
        for (unsigned long col = 0; col < this->max_sites_; ++col) {
            s->append_state(NucleotideSequence::missing_data_state);
        }
        this->sequence_storage_.push_back(s);
        this->available_sequences_.push(s);
    }
}

void NucleotideAlignment::clear() {
    for (auto & v : this->sequence_storage_) {
        if (v) {
            delete v;
        }
    }
    while (this->available_sequences_.size() > 0) {
        this->available_sequences_.pop();
    }
    this->sequence_storage_.clear();
    this->sequence_node_data_map_.clear();
    this->node_data_sequence_map_.clear();
}

void NucleotideAlignment::write_states_as_symbols(
        GeneNodeData * gene_node_data,
        std::ostream& out) const {
    auto si = this->node_data_sequence_map_.find(gene_node_data);
    TREESHREW_NDEBUG_ASSERT(si != this->node_data_sequence_map_.end());
    NucleotideSequence * seq = si->second;
    seq->write_states_as_symbols(out, seq->cbegin(), seq->cbegin()+this->num_active_sites_);
    // std::copy(seq->cbegin(),
    //         seq->cbegin() + this->num_active_sites_,
    //         std::ostream_iterator<double>(out, ""));
}

//////////////////////////////////////////////////////////////////////////////
// ShortReadSequences

ShortReadSequences::ShortReadSequences()
    : max_short_read_length_(0) {
}

ShortReadSequences::~ShortReadSequences() {
}

void ShortReadSequences::set(const NucleotideSequences& data) {
    this->short_reads_.reserve(data.size());
    for (std::vector<NucleotideSequence *>::const_iterator si = data.cbegin();
            si != data.cend();
            ++si) {
        NucleotideSequence * seq = *si;
    // for (auto seq : data) {
        this->short_reads_.emplace_back(*seq);
        if (seq->size() > this->max_short_read_length_) {
            this->max_short_read_length_ = seq->size();
        }
    }
}

} // namespace treeshrew


