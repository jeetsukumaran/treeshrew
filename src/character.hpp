#ifndef TREESHREW_CHARACTER_HPP
#define TREESHREW_CHARACTER_HPP

#include <array>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <numeric>    //inner_product
#include <functional> //plus, equal_to, not2
#include <gsl/gsl_randist.h>
#include "genetree.hpp"
#include "utility.hpp"

namespace treeshrew {

//////////////////////////////////////////////////////////////////////////////
// Typedefs

typedef int CharacterStateType;
typedef std::vector<CharacterStateType> CharacterStateVectorType;

//////////////////////////////////////////////////////////////////////////////
// Utility Functions

// Assumes that ``long_read_begin`` sequence is at least as long as ``short_read_begin`` to
// ``short_read_end``.
unsigned long hamming_distance(const CharacterStateVectorType::const_iterator& short_read_begin,
        const CharacterStateVectorType::const_iterator& short_read_end,
        const CharacterStateVectorType::const_iterator& long_read_begin);

// Assumes that ``long_read`` sequence is at least as long as ``short_read``
unsigned long sliding_hamming_distance(const CharacterStateVectorType& short_read,
        const CharacterStateVectorType& long_read);

//////////////////////////////////////////////////////////////////////////////
// NucleotideSequence

class NucleotideSequence {

    public:
        NucleotideSequence();
        NucleotideSequence(const std::string& label);
        ~NucleotideSequence();
        void reserve(unsigned long size) {
            this->sequence_.reserve(size);
        }
        inline unsigned long size() const {
            return this->sequence_.size();
        }
        inline CharacterStateVectorType::iterator begin() {
            return this->sequence_.begin();
        }
        inline CharacterStateVectorType::iterator end() {
            return this->sequence_.end();
        }
        inline const CharacterStateVectorType::const_iterator cbegin() const {
            return this->sequence_.cbegin();
        }
        inline const CharacterStateVectorType::const_iterator cend() const {
            return this->sequence_.cend();
        }
        inline std::vector<double>::iterator partials_begin() {
            return this->partials_.begin();
        }
        inline std::vector<double>::iterator partials_end() {
            return this->partials_.end();
        }
        inline const std::vector<double>::const_iterator partials_cbegin() const {
            return this->partials_.cbegin();
        }
        inline const std::vector<double>::const_iterator partials_cend() const {
            return this->partials_.cend();
        }
        inline void append_state(CharacterStateType c) {
            this->sequence_.push_back(c);
        }
        inline void append_state_by_symbol(char s) {
            auto state = NucleotideSequence::get_state_from_symbol(s);
            this->sequence_.push_back(state);
            auto state_partials = this->state_to_partials_map_.find(state)->second;
            this->partials_.insert(this->partials_.end(), state_partials.begin(), state_partials.end());
        }
        inline void append_states_by_symbols(const std::string& s) {
            this->sequence_.reserve(this->sequence_.size() + s.size());
            for (auto & c : s) {
                this->append_state_by_symbol(c);
            }
        }
        inline const CharacterStateType * state_data() const {
            return this->sequence_.data();
        }
        inline const double * partials_data() const {
            return this->partials_.data();
        }
        inline const std::string& get_label() const {
            return this->label_;
        }
        inline void set_label(const std::string& label) {
            this->label_ = label;
        }

        void write_states_as_symbols(std::ostream& out) const;
        std::string get_states_as_symbols() const;

    protected:
        std::string                 label_;
        CharacterStateVectorType    sequence_;
        std::vector<double>         partials_;

    public:
        static const std::map<char, CharacterStateType>                     symbol_to_state_map_;
        static const std::map<CharacterStateType, char>                     state_to_symbol_map_;
        static const std::map<CharacterStateType, std::array<double, 4>>    state_to_partials_map_;
        static const CharacterStateType                                     missing_data_state;
        static const std::array<double, 4>                                  missing_data_partials;

    public:
        inline static const CharacterStateType get_state_from_symbol(char s) {
            auto state_lookup = NucleotideSequence::symbol_to_state_map_.find(s);
            if (state_lookup == NucleotideSequence::symbol_to_state_map_.end()) {
                treeshrew_abort("Invalid state symbol '", s, "'");
                return 4;
            } else {
                return state_lookup->second;
            }
        }
        inline static const char get_symbol_from_state(CharacterStateType s) {
            auto symbol_lookup = NucleotideSequence::state_to_symbol_map_.find(s);
            if (symbol_lookup == NucleotideSequence::state_to_symbol_map_.end()) {
                treeshrew_abort("Invalid state: ", s);
            }
            return symbol_lookup->second;
        }

}; // NucleotideSequence

//////////////////////////////////////////////////////////////////////////////
// NucleotideSequences

class NucleotideSequences {

    public:
        NucleotideSequences();
        ~NucleotideSequences();
        void clear();
        inline std::vector<NucleotideSequence *>::iterator begin() {
            return this->sequences_.begin();
        }
        inline std::vector<NucleotideSequence *>::iterator end() {
            return this->sequences_.end();
        }
        inline const std::vector<NucleotideSequence *>::const_iterator cbegin() const {
            return this->sequences_.cbegin();
        }
        inline const std::vector<NucleotideSequence *>::const_iterator cend() const {
            return this->sequences_.cend();
        }
        inline NucleotideSequence * new_sequence(const std::string& label) {
            NucleotideSequence * v = new NucleotideSequence(label);
            this->sequences_.push_back(v);
            this->label_sequence_map_[label] = v;
            return v;
        }
        NucleotideSequence * get_sequence(unsigned long index) const {
            TREESHREW_ASSERT(index < this->sequences_.size());
            return this->sequences_[index];
        }
        NucleotideSequence * get_sequence(const std::string& label) const {
            auto label_sequence = this->label_sequence_map_.find(label);
            TREESHREW_ASSERT(label_sequence != this->label_sequence_map_.end());
            return label_sequence->second;
        }
        unsigned long get_num_sequences() {
            return this->sequences_.size();
        }
        unsigned long get_num_sites() {
            if (this->sequences_.size() > 0) {
                return this->sequences_[0]->size();
            } else {
                return 0;
            }
        }
        void set_tip_data(GeneTree * gene_tree);

    protected:
        std::vector<NucleotideSequence *>               sequences_;
        std::map<std::string, NucleotideSequence *>     label_sequence_map_;


}; // NucleotideSequences

//////////////////////////////////////////////////////////////////////////////
// ShortReadSequence

class ShortReadSequence {

    public:
        ShortReadSequence(const NucleotideSequence& seq)
            : label_(seq.get_label())
            , sequence_(seq.cbegin(), seq.cend()) {
            this->begin_ = this->sequence_.begin();
            this->end_ = this->sequence_.end();
            this->size_ = this->sequence_.size();
        }

        inline double calc_probability_of_sequence(
                const CharacterStateVectorType::const_iterator& long_read_begin,
                const CharacterStateVectorType::const_iterator& long_read_end,
                const NucleotideSequence * long_read, double error_probability) {
            CharacterStateVectorType::const_iterator start_pos = long_read_begin;
            CharacterStateVectorType::const_iterator stop_pos = long_read_end - this->size_ + 1;
            assert(stop_pos >= start_pos);
            unsigned long num_mismatches = 0;
            double prob = 0.0;
            while (start_pos < stop_pos) {
                num_mismatches = std::inner_product(
                        this->begin_, this->end_, start_pos,
                        0, std::plus<unsigned int>(),
                        std::not2(std::equal_to<CharacterStateVectorType::value_type>()));
                prob += gsl_ran_binomial_pdf(num_mismatches, error_probability, this->size_);
                ++start_pos;
            }
            return std::log(prob);
        }

    private:
        std::string                                 label_;
        CharacterStateVectorType                    sequence_;
        CharacterStateVectorType::const_iterator    begin_;
        CharacterStateVectorType::const_iterator    end_;
        unsigned long                               size_;

}; // ShortReadSequence

//////////////////////////////////////////////////////////////////////////////
// NucleotideAlignment

class NucleotideAlignment {

    public:
        NucleotideAlignment(unsigned long max_sequences,
                unsigned long max_sites);
        ~NucleotideAlignment();
        void create();
        void clear();
        void set_alignment(const NucleotideSequences& sequences);
        void set_gene_tree(GeneTree * gene_tree);
        inline unsigned long get_max_sequences() const {
            return this->max_sequences_;
        }
        inline unsigned long get_max_sites() const {
            return this->max_sites_;
        }
        inline unsigned long get_num_active_sites() const {
            return this->num_active_sites_;
        }
        inline void set_num_active_sites(unsigned long num) {
            this->num_active_sites_ = num;
        }
        inline void new_sequence(
                GeneNodeData * gene_node_data,
                const NucleotideSequence * src_seq=nullptr) {
            TREESHREW_NDEBUG_ASSERT(gene_node_data);
            if (this->available_sequences_.size() == 0) {
                treeshrew_abort("Maximum number of sequences exceeded");
            }
            NucleotideSequence * seq = this->available_sequences_.top();
            this->available_sequences_.pop();
            this->sequence_node_data_map_[seq] = gene_node_data;
            this->node_data_sequence_map_[gene_node_data] = seq;
            seq->set_label(gene_node_data->get_label());
            if (src_seq) {
                this->set_sequence_states(seq, src_seq);
            }
        }

    protected:
        void set_sequence_states(NucleotideSequence * seq, const NucleotideSequence * src_seq) {
            unsigned long len = src_seq->size();
            if (len > this->num_active_sites_) {
                this->num_active_sites_ = len;
            }
            std::copy(src_seq->cbegin(), src_seq->cend(), seq->begin());
            std::fill(seq->begin() + len, seq->end(), NucleotideSequence::missing_data_state);
            std::copy(src_seq->partials_cbegin(), src_seq->partials_cend(), seq->partials_begin());
            std::fill(seq->partials_begin() + len, seq->partials_end(), 1.0);
        }

    protected:
        unsigned long                                           max_sequences_;
        unsigned long                                           max_sites_;
        unsigned long                                           num_active_sites_;
        std::vector<NucleotideSequence *>                       sequence_storage_;
        std::stack<NucleotideSequence *>                        available_sequences_;
        std::map<NucleotideSequence *, GeneNodeData *>          sequence_node_data_map_;
        std::map<GeneNodeData *, NucleotideSequence *>          node_data_sequence_map_;

}; // NucleotideAlignment


} // namespace treeshrew

#endif
