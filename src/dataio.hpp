#ifndef TREESHREW_DATAIO_HPP
#define TREESHREW_DATAIO_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <ncl/nxsmultiformat.h>
#include "utility.hpp"

namespace treeshrew {

namespace sequenceio {

////////////////////////////////////////////////////////////////////////////////
// Sequence Reading Utility Functions

template <class SequencesType>
int read_from_stream(
        SequencesType& sequences,
        std::istream& src,
        const std::string& format) {
    MultiFormatReader reader(-1, NxsReader::IGNORE_WARNINGS);
    reader.SetWarningOutputLevel(NxsReader::AMBIGUOUS_CONTENT_WARNING);
    reader.SetCoerceUnderscoresToSpaces(true);
    const char * format_cstr = nullptr;
    if (format == "phylip" || format == "dnaphylip") {
        format_cstr = "dnarelaxedphylip";
    } else if (format == "rnaphylip") {
        format_cstr = "rnarelaxedphylip";
    } else if (format == "fasta") {
        format_cstr = "dnafasta";
    } else {
        format_cstr = format.c_str();
    }
    reader.ReadStream(src, format_cstr);
    unsigned num_taxa_blocks = reader.GetNumTaxaBlocks();
    NxsTaxaBlock *  taxa_block = reader.GetTaxaBlock(num_taxa_blocks-1);
    if (!taxa_block) {
        treeshrew_abort("No taxon definitions were parsed (invalid file format?)");
    }
    NxsCharactersBlock * chars_block = reader.GetCharactersBlock(taxa_block, 0);
    if (!chars_block) {
        treeshrew_abort("No character states were parsed (invalid file format?)");
    }
    unsigned int ntax = taxa_block->GetNTaxTotal();
    for (unsigned int taxon_idx = 0; taxon_idx < ntax; ++taxon_idx) {
        // const char * label = NxsString::GetEscaped(taxa_block->GetTaxonLabel(taxon_idx)).c_str();
        const char * label = taxa_block->GetTaxonLabel(taxon_idx).c_str();
        NxsDiscreteStateRow row = chars_block->GetDiscreteMatrixRow(taxon_idx);
        // seqs->reserve(row.size());
        std::ostringstream o;
        for (unsigned int col = 0; col < row.size(); ++col) {
            chars_block->ShowStateLabels(o, taxon_idx, col, 0);
        }
        auto * seq = sequences.new_sequence(label);
        seq->append_states_by_symbols(o.str());
        // for (auto & ncl_state : row) {
        //     seqs->append_state_by_symbol(symbols[ncl_state]);
        // }
    }
    return ntax;
}

template <class SequencesType>
int read_from_filepath(
        SequencesType& sequences,
        const std::string& filepath,
        const std::string& format) {
    std::ifstream f(filepath);
    if (!f.good()) {
        treeshrew_abort("Error opening file for input");
    }
    return read_from_stream(sequences, f, format);
}

template <class SequencesType>
int read_from_string(
        SequencesType& sequences,
        const std::string& str,
        const std::string& format) {
    std::istringstream s(str);
    return read_from_stream(sequences, s, format);
}

////////////////////////////////////////////////////////////////////////////////
// Sequence Writing Utility Functions

template <class SequencesType>
int write_fasta(
        SequencesType& sequences,
        std::ostream& out) {
    // TREESHREW_ASSERT(this->sequences_.size() == this->labels_.size());
    // for (unsigned long idx = 0; idx < this->sequences_.size(); ++idx) {
    //     out << ">" << this->labels_[idx] << std::endl;
    //     unsigned int col_count = 0;
    //     for (auto & ch : (*this->sequences_[idx]) ) {
    //         if (col_count == 70) {
    //             out << std::endl;
    //             col_count = 0;
    //         }
    //         out << NucleotideSequences::get_symbol_from_state(ch);
    //         ++col_count;
    //     }
    //     out << "\n" << std::endl;
    // }
    return 1;
}

} // namespace sequenceio

namespace treeio {

////////////////////////////////////////////////////////////////////////////////
// Tree Reading Utility Functions

template <class TreeType>
int construct_tree(
        TreeType * ttree,
        const NxsTaxaBlock * tb,
        const NxsFullTreeDescription & ftd) {
    NxsSimpleTree ncl_tree(ftd, -1, -1.0);
    auto * root = ttree->head_node();
    decltype(root) node_parent = nullptr;
    decltype(root) new_node = nullptr;
    std::vector<const NxsSimpleNode *> ncl_nodes = ncl_tree.GetPreorderTraversal();
    std::map<const NxsSimpleNode *, decltype(root)> ncl_to_native;
    int size = 0;
    for (auto & ncl_node : ncl_nodes) {
        const NxsSimpleEdge & ncl_edge = ncl_node->GetEdgeToParentRef();
        const NxsSimpleNode * ncl_par = ncl_edge.GetParent();
        std::vector<NxsSimpleNode *> ncl_child_nodes = ncl_node->GetChildren();
        unsigned int nchildren = ncl_child_nodes.size();
        double edge_len = ncl_edge.GetDblEdgeLen();
        if (edge_len < 0) {
            edge_len = 0.0;
        }
        std::string label;

        if (nchildren > 2) {
            treeshrew_abort("Tree has node with more than 2 children");
        } else if (nchildren == 1) {
            treeshrew_abort("Tree has node with only 1 child");
        } else if (nchildren == 0) {
            unsigned int ncl_taxon_idx = ncl_node->GetTaxonIndex();
            // label = NxsString::GetEscaped(tb->GetTaxonLabel(ncl_taxon_idx)).c_str();
            label = tb->GetTaxonLabel(ncl_taxon_idx).c_str();
            new_node = ttree->allocate_leaf_node();
        } else if (nchildren == 2) {
            label = NxsString::GetEscaped(ncl_node->GetName()).c_str();
            if (!ncl_par) {
                new_node = ttree->head_node();
            } else {
                new_node = ttree->allocate_internal_node();
            }
        }
        new_node->data().set_label(label);
        new_node->data().set_edge_length(edge_len);
        ncl_to_native[ncl_node] = new_node;
        if (ncl_par) {
            if (ncl_to_native.find(ncl_par) == ncl_to_native.end()) {
                treeshrew_abort("Parent node not visited in preorder traversal");
            }
            node_parent = ncl_to_native[ncl_par];
            if (!node_parent) {
                treeshrew_abort("Null parent node");
            }
            node_parent->add_child(new_node);
        }
        ++size;
    }
    return size;
}

template <class TreeType>
int read_from_stream(
        std::vector<TreeType *>& trees,
        std::istream& src,
        const std::string& format="nexus",
        unsigned long max_tips=0) {
    MultiFormatReader reader(-1, NxsReader::IGNORE_WARNINGS);
    reader.SetWarningOutputLevel(NxsReader::AMBIGUOUS_CONTENT_WARNING);
    reader.SetCoerceUnderscoresToSpaces(false);
    const char * format_cstr = nullptr;
    if (format == "newick") {
        format_cstr = "relaxedphyliptree";
    } else {
        format_cstr = format.c_str();
    }
    reader.ReadStream(src, format_cstr);
    unsigned num_taxa_blocks = reader.GetNumTaxaBlocks();
    NxsTaxaBlock *  taxa_block = reader.GetTaxaBlock(num_taxa_blocks-1);
    if (!taxa_block) {
        treeshrew_abort("No taxon definitions were parsed (invalid file format?)");
    }
    NxsTreesBlock * trees_block = reader.GetTreesBlock(taxa_block, 0);
    if (!trees_block) {
        return 0;
    }
    if (max_tips == 0) {
        max_tips = taxa_block->GetNTaxTotal();
    }
    unsigned int num_trees = trees_block->GetNumTrees();
    int tree_count = 0;
    for (unsigned int tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
        auto * tree = new TreeType(max_tips);
        const NxsFullTreeDescription & ftd = trees_block->GetFullTreeDescription(tree_idx);
        construct_tree<TreeType>(tree, taxa_block, ftd);
        trees.push_back(tree);
        ++tree_count;
    }
    return tree_count;
}

template <class TreeType>
int read_from_filepath(
        std::vector<TreeType *>& trees,
        const std::string& filepath,
        const std::string& format="nexus",
        unsigned long max_tips=0) {
    std::ifstream f(filepath);
    if (!f.good()) {
        treeshrew_abort("Error opening file for input");
    }
    return read_from_stream<TreeType>(trees, f, format, max_tips);
}

template <class TreeType>
int read_from_string(
        std::vector<TreeType *>& trees,
        const std::string& str,
        const std::string& format="nexus",
        unsigned long max_tips=0) {
    std::istringstream s(str);
    return read_from_stream<TreeType>(trees, s, format, max_tips);
}

////////////////////////////////////////////////////////////////////////////////
// Tree Writing Utility Functions

template <class TreeType, class iter>
void write_newick_node(TreeType * tree, const iter& tree_iter, std::ostream& out) {
    TREESHREW_NDEBUG_ASSERT(tree_iter);
    if (!tree_iter.is_leaf()) {
        out << "(";
        int ch_count = 0;
        for (typename TreeType::child_iterator chi = tree->children_begin(tree_iter);
                chi != tree->children_end(tree_iter);
                ++chi, ++ch_count) {
            if (ch_count > 0) {
                out << ", ";
            }
            write_newick_node(tree, chi, out);
        }
        out << ")";
    } else {
    }
    // label = NxsString::GetEscaped(tb->GetTaxonLabel(ncl_taxon_idx)).c_str();
    // out << tree_iter->get_label();
    out << NxsString::GetEscaped(tree_iter->get_label()).c_str();
    out << ":" << tree_iter->get_edge_length();
}

template <class TreeType>
void write_newick(TreeType * tree, std::ostream& out) {
    out << "[&R] ";
    write_newick_node(tree, tree->begin(), out);
    out << ";" << std::endl;
}


} // namespace treeio

} // namespace treeshrew

#endif
