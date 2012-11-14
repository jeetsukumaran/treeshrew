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

std::map<char, CharacterStateType> NucleotideSequence::symbol_to_state_map_ {
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

std::map<CharacterStateType, char> NucleotideSequence::state_to_symbol_map_ {
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

std::map<CharacterStateType, std::array<double, 4>> NucleotideSequence::state_to_partials_map_ {
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

NucleotideSequence::NucleotideSequence() {
}

NucleotideSequence::NucleotideSequence(const std::string& label) :
    label_(label) {
}

NucleotideSequence::~NucleotideSequence() {
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

} // namespace treeshrew


