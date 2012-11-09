#ifndef TREESHREW_CHARACTER_HPP
#define TREESHREW_CHARACTER_HPP

#include <array>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include "genetree.hpp"
#include "utility.hpp"

namespace treeshrew {

typedef int CharacterStateType;
typedef std::vector<CharacterStateType> CharacterStateVectorType;

class NucleotideSequence {

    public:
        NucleotideSequence();
        NucleotideSequence(const std::string& label);
        ~NucleotideSequence();
        void reserve(unsigned long size) {
            this->sequence_.reserve(size);
        }
        unsigned long size() {
            return this->sequence_.size();
        }
        void append_state(CharacterStateType c) {
            this->sequence_.push_back(c);
        }
        void append_state_by_symbol(char s) {
            auto state = NucleotideSequence::get_state_from_symbol(s);
            this->sequence_.push_back(state);
            auto state_partials = this->state_to_partials_map_[state];
            this->partials_.insert(this->partials_.end(), state_partials.begin(), state_partials.end());
        }
        void append_states_by_symbols(const std::string& s) {
            this->sequence_.reserve(this->sequence_.size() + s.size());
            for (auto & c : s) {
                this->append_state_by_symbol(c);
            }
        }
        CharacterStateType * state_data() {
            return this->sequence_.data();
        }
        double * partials_data() {
            return this->partials_.data();
        }
        const std::string& get_label() const {
            return this->label_;
        }
        void write_states_as_symbols(std::ostream& out) const {
            for (auto & s : this->sequence_) {
                out << this->state_to_symbol_map_[s];
            }
        }
        std::string get_states_as_symbols() const {
            std::ostringstream out;
            this->write_states_as_symbols(out);
            return out.str();
        }

    private:
        std::string                 label_;
        CharacterStateVectorType    sequence_;
        std::vector<double>         partials_;

    protected:
        static std::map<char, CharacterStateType>                     symbol_to_state_map_;
        static std::map<CharacterStateType, char>                     state_to_symbol_map_;
        static std::map<CharacterStateType, std::array<double, 4>>    state_to_partials_map_;

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


class NucleotideSequences {

    public:
        NucleotideSequences();
        ~NucleotideSequences();
        void clear();
        inline NucleotideSequence * new_sequence(const std::string& label) {
            NucleotideSequence * v = new NucleotideSequence(label);
            this->sequences_.push_back(v);
            this->label_sequence_map_[label] = v;
            return v;
        }
        NucleotideSequence * get_sequence(unsigned long index) {
            TREESHREW_ASSERT(index < this->sequences_.size());
            return this->sequences_[index];
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


} // namespace treeshrew

#endif
