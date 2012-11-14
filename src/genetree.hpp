#ifndef TREESHREW_GENETREE_HPP
#define TREESHREW_GENETREE_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <stack>
#include <map>
#include <libhmsbeagle/beagle.h>
#include "utility.hpp"
#include "tree.hpp"

namespace treeshrew {

////////////////////////////////////////////////////////////////////////////////
// GeneTreeNode
class GeneNodeData {

    public:
        GeneNodeData(int index = -1) :
            index_(index),
            edge_length_(0.0),
            is_dirty_(true),
            ln_likelihood_(0.0) { }

        inline void clear() {
            this->edge_length_ = 0.0;
            this->is_dirty_ = true;
            this->ln_likelihood_ = 0.0;
        }

        inline void set(int index, double edge_length_, const std::string& label, bool is_dirty) {
            this->index_ = index;
            this->edge_length_ = edge_length_;
            this->label_ = label;
            this->is_dirty_ = is_dirty;
        }
        inline void set_index(int index) {
            this->index_ = index;
        }
        inline void set_edge_length(double edge_length) {
            this->edge_length_ = edge_length;
        }
        inline void set_label(const std::string& label) {
            this->label_ = label;
        }
        inline void set_dirty(bool dirty=true) {
            this->is_dirty_ = dirty;
        }
        inline void flag_as_dirty() {
            this->is_dirty_ = true;
        }
        inline int get_index() const {
            return this->index_;
        }
        inline double get_edge_length() const {
            return this->edge_length_;
        }
        inline const std::string& get_label() const {
            return this->label_;
        }
        inline bool is_dirty() const {
            return this->is_dirty_;
        }

    private:
        int             index_;
        double          edge_length_;
        std::string     label_;
        bool            is_dirty_;
        double          ln_likelihood_;

}; // GeneNodeData

////////////////////////////////////////////////////////////////////////////////
// RestrictedResourceAllocator
template <class T>
class RestrictedResourceAllocator {
    public:
        RestrictedResourceAllocator(unsigned long max_size, int index_offset=0)
                : index_offset_(index_offset) {
            this->reserve(max_size);
        }

        ~RestrictedResourceAllocator() {
            this->clear();
        }
        void clear(){
            while (this->available_.size() > 0) {
                this->available_.pop();
            }
            for (auto nd : this->storage_) {
                if (nd) {
                    delete nd;
                }
            }
            this->storage_.clear();
        }
        unsigned long add_storage(unsigned long count) {
            for (unsigned long index = 0; index < count; ++index) {
                int item_index = this->index_offset_ + this->storage_.size();
                T * item = new T();
                item->data().set_index(item_index);
                this->storage_.push_back(item);
                this->available_.push(item);
            }
            return this->storage_.size();
        }
        unsigned long reserve(unsigned long total){
            if (total <= this->storage_.size()) {
                return 0;
            }
            unsigned long storage_to_add = total - this->storage_.size();
            this->add_storage(storage_to_add);
            return storage_to_add;
        }
        inline T * allocate() {
            if (this->available_.size() == 0) {
                treeshrew_abort("Node reserves exhausted");
            }
            T * nd = this->available_.top();
            this->available_.pop();
            return nd;
        }
        inline void deallocate(T * nd) {
            this->available_.push(nd);
        }
    private:
        int               index_offset_;
        std::vector<T *>  storage_;
        std::stack<T *>   available_;
}; // RestrictedResourceAllocator

////////////////////////////////////////////////////////////////////////////////
// GeneTree

class GeneTree : public Tree<GeneNodeData> {

    public:
        typedef TreeNode<GeneNodeData> GeneTreeNode;

    public:
        GeneTree(unsigned long max_tips);
        ~GeneTree();

        inline GeneTreeNode * allocate_leaf_node() {
            return this->leaf_node_allocator_.allocate();
        }

        inline void deallocate_leaf_node(GeneTreeNode * nd) {
            return this->leaf_node_allocator_.deallocate(nd);
        }

        inline GeneTreeNode * allocate_internal_node() {
            return this->internal_node_allocator_.allocate();
        }

        inline void deallocate_internal_node(GeneTreeNode * nd) {
            return this->internal_node_allocator_.deallocate(nd);
        }

        void clear();

        int create_beagle_instance(int num_sites);
        int set_tip_states(const GeneNodeData& tip, const int * data);
        int set_tip_partials(const GeneNodeData& tip, const double * data);
        double calc_ln_probability();
        void free_beagle_instance();

    private:
        unsigned long                              max_tips_;
        RestrictedResourceAllocator<GeneTreeNode>  leaf_node_allocator_;
        RestrictedResourceAllocator<GeneTreeNode>  internal_node_allocator_;
        int                                        beagle_instance_;
        BeagleInstanceDetails *                    beagle_return_info_;


}; // GeneTree


} // namespace treeshrew

#endif
