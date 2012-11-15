#ifndef TREESHREW_TREE_HPP
#define TREESHREW_TREE_HPP

#include <cassert>
#include <limits>
#include <stack>
#include <iostream>
#include <stdexcept>

namespace treeshrew {

////////////////////////////////////////////////////////////////////////////////
// TreeNode

template<class NodeData>
class TreeNode {

    public:

        /////////////////////////////////////////////////////////////////////////
        // Lifecycle

        TreeNode(size_t owner = 0)
            : parent_(nullptr)
            , first_child_(nullptr)
            , last_child_(nullptr)
            , next_sibling_(nullptr) { }

        TreeNode(const NodeData& data, size_t owner = 0)
            : parent_(nullptr)
            , first_child_(nullptr)
            , last_child_(nullptr)
            , next_sibling_(nullptr)
            , data_(data) { }

        ~TreeNode() {
            this->clear_links();
        }

        /////////////////////////////////////////////////////////////////////////
        // Structure

        inline void add_child(TreeNode<NodeData> * ch) {
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

        inline TreeNode<NodeData> * parent_node() const {
            return this->parent_;
        }

        inline void set_parent(TreeNode<NodeData> * parent) {
            this->parent_ = parent;
        }

        inline TreeNode<NodeData> * first_child_node() const {
            return this->first_child_;
        }

        inline TreeNode<NodeData> * last_child_node() const {
            return this->last_child_;
        }

        inline TreeNode<NodeData> * next_sibling_node() const {
            return this->next_sibling_;
        }
        inline void set_next_sibling_node(TreeNode<NodeData> * nd) {
            this->next_sibling_ = nd;
        }

        bool is_leaf() const {
            return this->first_child_ == nullptr;
        }

        inline void clear_links() {
            this->parent_ = nullptr;
            this->first_child_ = nullptr;
            this->last_child_ = nullptr;
            this->next_sibling_ = nullptr;
        }

        /////////////////////////////////////////////////////////////////////////
        // Data

        const NodeData& data() const {
            return this->data_;
        }

        NodeData& data() {
            return this->data_;
        }

        void set_data(const NodeData& data) {
            this->data_ = data;
        }

        inline NodeData& parent() const {
            assert(this->parent_);
            return this->parent_->data;
        }

        inline NodeData& first_child() const {
            assert(this->first_child_);
            return this->first_child_->data_;
        }

        inline NodeData& last_child() const {
            assert(this->last_child_);
            return this->last_child_->data_;
        }

        inline NodeData& next_sibling() const {
            assert(this->next_sibling_);
            return this->next_sibling_->data;
        }

    private:
        TreeNode<NodeData> *        parent_;
        TreeNode<NodeData> *        first_child_;
        TreeNode<NodeData> *        last_child_;
        TreeNode<NodeData> *        next_sibling_;
        NodeData                    data_;

}; // TreeNode

////////////////////////////////////////////////////////////////////////////////
// Tree

template<class NodeData, class TreeNodeAllocator = std::allocator<TreeNode<NodeData>>>
class Tree {

    public:
        typedef TreeNode<NodeData> node_type;
        typedef NodeData value_type;

    public:

        /////////////////////////////////////////////////////////////////////////
        // Lifecycle

        // If ``manage_node_allocation`` is false, then client code is
        // responsible for creation (= allocation of memory + construction) +
        // disposal (= destruction + deallocation) of *all* node objects.
        //
        // All structural methods that require the base class to create nodes
        // will be unavailable (std::logic_error will be thrown).
        //
        // In addition, client code *must* call Tree::create(), passing in two
        // constructed node objects to service as the tree before any
        // operations.
        Tree(bool manage_node_allocation=true)
                : manage_node_allocation_(manage_node_allocation) {
            if (this->manage_node_allocation_) {
                this->create(this->acquire_node(), this->acquire_node());
            }
        }

        virtual ~Tree() {
            this->release_node(this->head_node_);
            this->release_node(this->stop_node_);
        }

        // If node allocation is non-managed (``manage_node_allocation ==
        // false``), then client code *must* call this method, passing in two
        // node objects before this class is used.
        void create(node_type * head_node, node_type * stop_node) {
            this->head_node_ = head_node;
            this->stop_node_ = stop_node;
            this->head_node_->set_next_sibling_node(this->stop_node_);
        }

        void clear();

        /////////////////////////////////////////////////////////////////////////
        // Allocators

        inline node_type * acquire_node() {
            if (this->manage_node_allocation_) {
                node_type * nd = this->tree_node_allocator_.allocate(1, 0);
                this->tree_node_allocator_.construct(nd, node_type());
                return nd;
            } else {
                throw std::logic_error("Tree::acquire_node(): Request for node allocation but resource is not managed");
            }
        }

        inline node_type * acquire_node(const value_type& data) {
            if (this->manage_node_allocation_) {
                node_type * nd = this->tree_node_allocator_.allocate(1, 0);
                this->tree_node_allocator_.construct(nd, node_type(data));
                return nd;
            } else {
                throw std::logic_error("Tree::acquire_node(const value_type& data): Request for node allocation but resource is not managed");
            }
        }

        template <typename... Types>
        inline node_type * acquire_emplaced_data_node(const Types&... args) {
            if (this->manage_node_allocation_) {
                node_type * nd = this->tree_node_allocator_.allocate(1, 0);
                this->tree_node_allocator_.construct(nd, node_type(value_type (args...)));
                return nd;
            } else {
                throw std::logic_error("Tree::acquire_emplaced_data_node(const Types&... args): Request for node allocation but resource is not managed");
            }
        }

        inline void release_node(node_type * nd) {
            if (this->manage_node_allocation_) {
                this->tree_node_allocator_.destroy(nd);
                this->tree_node_allocator_.deallocate(nd, 1);
            }
            // } else {
            //     throw std::logic_error("Tree::release_node(node_type * nd): Request for node allocation but resource is not managed");
            // }
        }

        /////////////////////////////////////////////////////////////////////////
        // Iterators

        class base_iterator {
            public:
				typedef base_iterator               self_type;
				typedef NodeData                    value_type;
				typedef NodeData *                  pointer;
				typedef NodeData &                  reference;
				typedef unsigned long               size_type;
				typedef int                         difference_type;
				typedef std::forward_iterator_tag   iterator_category;

                base_iterator(node_type * node)
                    : node_(node) { }
                virtual ~base_iterator() {}
                reference operator*() const {
                    return this->node_->data();
                }
                pointer operator->() const {
                    return &(this->node_->data());
                }
                bool operator==(const self_type& rhs) const {
                    return this->node_ == rhs.node_;
                }
                bool operator!=(const self_type& rhs) const {
                    return !(*this == rhs);
                }
                node_type * node() const {
                    return this->node_;
                }
                bool is_leaf() const {
                    assert(this->node_);
                    return this->node_->is_leaf();
                }
                reference parent() const {
                    assert(this->node_);
                    return this->node_->parent();
                }
                reference first_child() const {
                    assert(this->node_);
                    return this->node_->first_child();
                }
                reference last_child() const {
                    assert(this->node_);
                    return this->node_->last_child();
                }
                reference next_sibling() const {
                    assert(this->node_);
                    return this->node_->next_sibling();
                }
                node_type * parent_node() const {
                    assert(this->node_);
                    return this->node_->parent_node();
                }
                node_type * first_child_node() const {
                    assert(this->node_);
                    return this->node_->first_child_node();
                }
                node_type * last_child_node() const {
                    assert(this->node_);
                    return this->node_->last_child_node();
                }
                node_type * next_sibling_node() const {
                    assert(this->node_);
                    return this->node_->next_sibling_node();
                }
                bool operator!() const {
                    return this->node_ == nullptr;
                }
            protected:
                node_type * node_;
        }; // base_iterator
        typedef base_iterator iterator;

        class preorder_iterator : public base_iterator {
            public:
                typedef preorder_iterator self_type;
                preorder_iterator(node_type * node)
                    : base_iterator(node)
                    , skip_current_children_(false) { }
                virtual ~preorder_iterator() {}
                const self_type& operator++() {
                    assert(this->node_ != nullptr);
                    if (this->skip_current_children_ || this->node_->first_child_node() == nullptr) {
                        this->skip_current_children_ = false;
                        while (this->node_->next_sibling_node() == nullptr) {
                            this->node_ = this->node_->parent_node();
                            if (this->node_ == nullptr) {
                                return *this;
                            }
                        }
                        this->node_ = this->node_->next_sibling_node();
                    } else {
                        this->node_ = this->node_->first_child_node();
                    }
                    return *this;
                }
                self_type operator++(int) {
                    self_type i = *this;
                    ++(*this);
                    return i;
                }
            private:
                bool skip_current_children_;
        }; // preorder_iterator

        class postorder_iterator : public base_iterator {
            public:
                typedef postorder_iterator self_type;
                postorder_iterator(node_type * node)
                    : base_iterator(node) { }
                virtual ~postorder_iterator() {}
                const self_type& operator++() {
                    assert(this->node_ != nullptr);
                    node_type * nd = this->node_->next_sibling_node();
                    if (nd == nullptr) {
                        this->node_ = this->node_->parent_node();
                    } else {
                        this->node_ = nd;
                        nd = this->node_->first_child_node();
                        while (nd) {
                            this->node_= nd;
                            nd = this->node_->first_child_node();
                        }
                    }
                    return *this;
                }
                self_type operator++(int) {
                    self_type i = *this;
                    ++(*this);
                    return i;
                }
        }; // postorder_iterator

        class leaf_iterator : public base_iterator {

            public:
                typedef leaf_iterator self_type;
                leaf_iterator(node_type * node, node_type * top_node = nullptr)
                    : base_iterator(node),
                      top_node_(top_node) { }
                virtual ~leaf_iterator() {}
                const self_type& operator++() {
                    if (this->node_->first_child_node() != nullptr) {
                        while (this->node_->first_child_node()) {
                            this->node_ = this->node_->first_child_node();
                        }
                    } else {
                        while (this->node_->next_sibling_node() == nullptr) {
                            if (this->node_->parent_node() == nullptr) {
                                return *this;
                            }
                            this->node_ = this->node_->parent_node();
                            if (this->top_node_ != nullptr && this->node_ == this->top_node_) {
                                return *this;
                            }
                        }
                        this->node_ = this->node_->next_sibling_node();
                        while (this->node_->first_child_node()) {
                            this->node_ = this->node_->first_child_node();
                        }
                    }
                    return *this;
                }
                self_type operator++(int) {
                    self_type i = *this;
                    ++(*this);
                    return i;
                }
            private:
                node_type * top_node_;
        }; // leaf_iterator

        class child_iterator : public base_iterator {

            public:
                typedef child_iterator self_type;
                child_iterator(node_type * node)
                    : base_iterator(node) {
                    }
                virtual ~child_iterator() {}
                const self_type& operator++() {
                    if (this->node_ != nullptr) {
                        this->node_ = this->node_->next_sibling_node();
                    }
                    return *this;
                }
                self_type operator++(int) {
                    self_type i = *this;
                    ++(*this);
                    return i;
                }
        }; // child_iterator

        /////////////////////////////////////////////////////////////////////////
        // Iteration Control


        // -- preorder iterator --

        preorder_iterator begin() const {
            return preorder_iterator(this->head_node_);
        }

        preorder_iterator end() const {
            return preorder_iterator(this->stop_node_);
        }

        preorder_iterator preorder_begin() const {
            return preorder_iterator(this->head_node_);
        }

        preorder_iterator preorder_end() const {
            return preorder_iterator(this->stop_node_);
        }

        // -- postorder iterator --

        postorder_iterator postorder_begin() const {
            node_type * nd = this->head_node_;
            while (nd->first_child_node() != nullptr) {
                nd = nd->first_child_node();
            }
            return postorder_iterator(nd);
        }

        postorder_iterator postorder_end() const {
            return postorder_iterator(this->head_node_->next_sibling_node());
        }

        // -- leaf iterator --

        leaf_iterator leaf_begin() const {
            return this->leaf_begin(this->head_node_);
        }

        leaf_iterator leaf_begin(node_type * nd) const {
            node_type * first = nd;
            while (first->first_child_node() != nullptr) {
                first = first->first_child_node();
            }
            return leaf_iterator(first, nd);
        }

        template <typename iter>
        leaf_iterator leaf_begin(iter& pos) const {
            return leaf_begin(pos.node());
        }

        leaf_iterator leaf_end() const {
            return this->leaf_end(this->head_node_);
        }

        leaf_iterator leaf_end(node_type * nd) const {
            return leaf_iterator(nd);
        }

        template <typename iter>
        leaf_iterator leaf_end(iter& pos) const {
            return leaf_iterator(pos->node());
        }

        // -- children iterator --

        child_iterator children_begin() const {
            return child_iterator(this->head_node_->first_child_node());
        }

        child_iterator children_begin(node_type * nd) const {
            return child_iterator(nd->first_child_node());
        }

        template <typename iter>
        child_iterator children_begin(iter& pos) const {
            return child_iterator(pos.node()->first_child_node());
        }

        child_iterator children_end() const {
            return child_iterator(nullptr);
        }

        child_iterator children_end(node_type * nd) const {
            return child_iterator(nullptr);
        }

        template <typename iter>
        child_iterator children_end(iter& pos) const {
            return child_iterator(nullptr);
        }

        /////////////////////////////////////////////////////////////////////////
        // Metrics

        inline unsigned long get_num_leaves() const {
            unsigned long leaf_count = 0;
            for (auto ndi = this->leaf_begin(); ndi != this->leaf_end(); ++ndi, ++leaf_count) {
            }
            return leaf_count;
        }


        /////////////////////////////////////////////////////////////////////////
        // Structure Access

        node_type * head_node() const {
            return this->head_node_;
        }

        node_type * stop_node() const {
            return this->stop_node_;
        }

        /////////////////////////////////////////////////////////////////////////
        // Structure Manipulation

        template<typename iter> iter add_child(iter& pos, node_type * node) {
            pos.node()->add_child(node);
            return iter(node);
        }

        template<typename iter> iter add_child(iter& pos) {
            return this->add_child(pos, this->acquire_node());
        }

        template<typename iter> iter add_child(iter& pos, const value_type& data) {
            return this->add_child(pos, this->acquire_node(data));
        }

        template<typename iter, typename... Types> iter emplace_child(iter& pos, const Types&... args) {
            return this->add_child(pos, this->acquire_emplaced_data_node(args...));
        }

    protected:
        TreeNodeAllocator           tree_node_allocator_;
        bool                        manage_node_allocation_;
        node_type *                 head_node_;
        node_type *                 stop_node_;

}; // Tree

} // namespace treeshrew

#endif
