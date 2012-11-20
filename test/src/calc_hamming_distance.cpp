#include <iostream>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include <vector>
#include "../../src/character.hpp"

using namespace treeshrew;

// int main() {
//     GeneNodeData gnd;
//     NucleotideSequence sr_seq;
//     sr_seq.append_states_by_symbols("AACCGGTT");
//     ShortReadSequence short_read(sr_seq);
//     NucleotideSequence lr_seq;
//     lr_seq.append_states_by_symbols("AACCGGTT");
//     NucleotideAlignment a(8, 8);
//     a.new_sequence(&gnd, &lr_seq);
//     std::cout << std::log(a.calc_probability_of_sequence(&gnd, short_read, 0.5)) << std::endl;
// }

void get_sliding_windows(
        std::vector<CharacterStateVectorType>& window_sets,
        const CharacterStateVectorType& short_read,
        const CharacterStateVectorType& long_read) {
    CharacterStateVectorType::const_iterator start_pos = long_read.begin();
    CharacterStateVectorType::const_iterator stop_pos = long_read.end() - short_read.size() + 1;
    assert(stop_pos >= start_pos);
    while (start_pos < stop_pos) {
        window_sets.emplace_back(start_pos, start_pos + short_read.size());
        ++start_pos;
    }
}

int main() {
    CharacterStateVectorType  short_read{1,1,2};
    CharacterStateVectorType  long_read{1,1,2,2,1,1,2,3,1,1,2,3};
    std::vector<CharacterStateVectorType> window_sets;
    get_sliding_windows(window_sets, short_read, long_read);
    unsigned long total_diff_count = 0;
    for (auto & wset : window_sets) {
        unsigned long diff_count = 0;
        for (unsigned int i = 0; i < short_read.size(); ++ i) {
            if (short_read.at(i) != wset.at(i)) {
                diff_count += 1;
            }
            std::cerr << wset.at(i);
        }
        total_diff_count += diff_count;
        std::cerr << "   " << diff_count;
        // std::copy(wset.begin(), wset.end(), std::ostream_iterator<CharacterStateType>(std::cerr, ""));
        std::cerr << std::endl;
    }
    unsigned long hd = sliding_hamming_distance(short_read, long_read);
    std::cerr << "--" << std::endl;
    std::cerr << "Test target returned: " << hd << std::endl;
    std::cerr << "Test program returned: " << total_diff_count << std::endl;
    if (total_diff_count != hd) {
        exit(1);
    } else {
        exit(0);
    }
}
