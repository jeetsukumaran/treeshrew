#ifndef TREESHREW_TREE_HPP
#define TREESHREW_TREE_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <stack>
#include <map>
#include <libhmsbeagle/beagle.h>
#include "utility.hpp"

namespace treeshrew {

////////////////////////////////////////////////////////////////////////////////
// GeneTreeNode
class GeneTreeNode {

    public:
        GeneTreeNode(GeneTreeNode * parent = nullptr, int index = -1) :
            parent_(parent),
            first_child_(nullptr),
            last_child_(nullptr),
            next_sibling_(nullptr),
            index_(index),
            edge_length_(0.0),
            is_dirty_(true),
            ln_likelihood_(0.0) { }

        inline void add_child(GeneTreeNode * ch) {
            if (this->first_child_ == nullptr) {
                this->first_child_ = ch;
                this->last_child_ = ch;
            } else {
                this->last_child_->next_sibling_ = ch;
                this->last_child_ = ch;
            }
            ch->parent_ = this;
            ch->next_sibling_ = nullptr;
        }

        inline GeneTreeNode * get_parent() const {
            return this->parent_;
        }

        inline GeneTreeNode * get_first_child() const {
            return this->first_child_;
        }

        inline GeneTreeNode * get_last_child() const {
            return this->last_child_;
        }

        inline GeneTreeNode * get_next_sibling() const {
            return this->next_sibling_;
        }
        inline void set_next_sibling(GeneTreeNode * nd) {
            this->next_sibling_ = nd;
        }

        bool is_leaf() const {
            return this->first_child_ == nullptr;
        }

        inline void clear() {
            this->parent_ = nullptr;
            this->first_child_ = nullptr;
            this->last_child_ = nullptr;
            this->next_sibling_ = nullptr;
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

        void get_children(std::vector<GeneTreeNode *>& children) {
            GeneTreeNode * nd = this->first_child_;
            while (nd) {
                children.push_back(nd);
                nd = nd->next_sibling_;
            }
        }

        class base_iterator {

            public:
                typedef base_iterator self_type;
                typedef GeneTreeNode * value_type;
                typedef std::forward_iterator_tag iterator_category;
                typedef int difference_type;
                base_iterator(value_type node)
                    : node_(node) { }
                virtual ~base_iterator() {}
                value_type& operator*() {
                    return this->node_;
                }
                value_type* operator->() {
                    return &(this->node_);
                }
                bool operator==(const self_type& rhs) {
                    return this->node_ == rhs.node_;
                }
                bool operator!=(const self_type& rhs) {
                    return !(*this == rhs);
                }
            protected:
                value_type node_;
        }; // base_iterator

        class child_iterator : public base_iterator {

            public:
                typedef child_iterator self_type;
                typedef GeneTreeNode * value_type;
                typedef std::forward_iterator_tag iterator_category;
                typedef int difference_type;
                child_iterator(value_type node)
                    : base_iterator(node) { }
                virtual ~child_iterator() {}
                const self_type& operator++() {
                    this->node_ = this->node_->get_next_sibling();
                    return *this;
                }
                self_type operator++(int) {
                    ++(*this);
                    return *this;
                }
        }; // child_iterator

        child_iterator children_begin() {
            return child_iterator(this->first_child_);
        }

        child_iterator children_end() {
            // return child_iterator(this);
            return child_iterator(nullptr);
        }

        class postorder_iterator : public base_iterator {

            public:
                typedef postorder_iterator self_type;
                typedef GeneTreeNode * value_type;
                typedef std::forward_iterator_tag iterator_category;
                typedef int difference_type;
                postorder_iterator(value_type node)
                    : base_iterator(node) { }
                virtual ~postorder_iterator() {}
                const self_type& operator++() {
                    GeneTreeNode * nd = this->node_->get_next_sibling();
                    if (nd == nullptr) {
                        // if (!this->node_->get_parent()) {
                        //     std::cerr << "** " << this->node_->get_label() << " " << (void*)this->node_ << ", " << (void*)this->node_->get_parent() << std::endl;
                        // }
                        this->node_ = this->node_->get_parent();
                    } else {
                        this->node_ = nd;
                        nd = this->node_->get_first_child();
                        while (nd) {
                            this->node_= nd;
                            nd = this->node_->get_first_child();
                        }
                    }
                    return *this;
                }
                self_type operator++(int) {
                    ++(*this);
                    return *this;
                }
        }; // postorder_iterator

        postorder_iterator postorder_begin() {
            GeneTreeNode * nd = this;
            while (nd->first_child_ != nullptr) {
                nd = nd->first_child_;
            }
            return postorder_iterator(nd);
        }

        postorder_iterator postorder_end() {
            return postorder_iterator(this->next_sibling_);
        }

        class leaf_iterator : public base_iterator {

            public:
                typedef leaf_iterator self_type;
                typedef GeneTreeNode * value_type;
                typedef std::forward_iterator_tag iterator_category;
                typedef int difference_type;
                leaf_iterator(value_type node, value_type top_node = nullptr)
                    : base_iterator(node),
                      top_node_(top_node) { }
                virtual ~leaf_iterator() {}

                // assert(this->node!=0);
                // if(this->node->first_child!=0) { // current node is no longer leaf (children got added)
                //     while(this->node->first_child)
                //         this->node=this->node->first_child;
                // }
                // else {
                //     while(this->node->next_sibling==0) {
                //         if (this->node->parent==0) return *this;
                //         this->node=this->node->parent;
                //         if (top_node != 0 && this->node==top_node) return *this;
                //     }
                //     this->node=this->node->next_sibling;
                //     while(this->node->first_child)
                //         this->node=this->node->first_child;
                // }
                // return *this;

                const self_type& operator++() {
                    if (this->node_->get_first_child() != nullptr) {
                        while (this->node_->get_first_child()) {
                            this->node_ = this->node_->get_first_child();
                        }
                    } else {
                        while (this->node_->get_next_sibling() == nullptr) {
                            if (this->node_->get_parent() == nullptr) {
                                std::cerr << "***** " << __LINE__ << std::endl;
                                return *this;
                            }
                            this->node_ = this->node_->get_parent();
                            if (this->top_node_ != nullptr && this->node_ == this->top_node_) {
                                return *this;
                            }
                        }
                        this->node_ = this->node_->get_next_sibling();
                        while (this->node_->get_first_child()) {
                            this->node_ = this->node_->get_first_child();
                        }
                    }
                    return *this;
                }
                self_type operator++(int) {
                    ++(*this);
                    return *this;
                }
            private:
                value_type   top_node_;
        }; // leaf_iterator

        leaf_iterator leaf_begin() {
            GeneTreeNode * nd = this;
            while (nd->first_child_ != nullptr) {
                nd = nd->first_child_;
            }
            return leaf_iterator(nd, this);
        }

        leaf_iterator leaf_end() {
            return leaf_iterator(this);
        }

    private:
        GeneTreeNode *      parent_;
        GeneTreeNode *      first_child_;
        GeneTreeNode *      last_child_;
        GeneTreeNode *      next_sibling_;
        int             index_;
        double          edge_length_;
        std::string     label_;
        bool            is_dirty_;
        double          ln_likelihood_;

}; // GeneTreeNode

////////////////////////////////////////////////////////////////////////////////
// GeneTreeNodeAllocator
class GeneTreeNodeAllocator {
    public:
        GeneTreeNodeAllocator(unsigned long max_nodes, int index_offset=0);
        ~GeneTreeNodeAllocator();
        void clear();
        unsigned long add_storage(unsigned long count);
        unsigned long reserve(unsigned long total);
        inline GeneTreeNode * allocate() {
            if (this->available_nodes_.size() == 0) {
                treeshrew_abort("Node reserves exhausted");
            }
            GeneTreeNode * nd = this->available_nodes_.top();
            this->available_nodes_.pop();
            return nd;
        }
        inline void deallocate(GeneTreeNode * nd) {
            nd->clear();
            this->available_nodes_.push(nd);
        }
    private:
        int                             index_offset_;
        std::vector<GeneTreeNode *>     node_storage_;
        std::stack<GeneTreeNode *>      available_nodes_;
}; // GeneTreeNodeAllocator

////////////////////////////////////////////////////////////////////////////////
// GeneTree

class GeneTree {

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

        inline GeneTreeNode * get_seed_node() {
            return this->seed_node_;
        }

        GeneTreeNode::postorder_iterator postorder_begin() {
            return this->seed_node_->postorder_begin();
        }

        GeneTreeNode::postorder_iterator postorder_end() {
            return this->seed_node_->postorder_end();
        }

        GeneTreeNode::leaf_iterator leaf_begin() {
            return this->seed_node_->leaf_begin();
        }

        GeneTreeNode::leaf_iterator leaf_end() {
            return this->seed_node_->leaf_end();
        }

        int create_beagle_instance(int num_sites);
        int set_tip_states(GeneTreeNode * tip, int * data);
        int set_tip_partials(GeneTreeNode * tip, double * data);
        double calc_ln_probability();
        void free_beagle_instance();

    private:
        unsigned long               max_tips_;
        GeneTreeNode *              seed_node_;
        GeneTreeNode *              foot_node_;
        GeneTreeNodeAllocator       leaf_node_allocator_;
        GeneTreeNodeAllocator       internal_node_allocator_;
        int                         beagle_instance_;
        BeagleInstanceDetails *     beagle_return_info_;


}; // GeneTree


} // namespace treeshrew

#endif
