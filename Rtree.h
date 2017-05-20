#ifndef PCM_R_TREE_RTREE_H
#define PCM_R_TREE_RTREE_H

#include <ostream>
#include <algorithm>
#include <functional>

#include "MBR.h"
#include "RTObject.h"

namespace rtree
{
    namespace details
    {
        template<class ObjectT, class BoundValueT, uint16_t dimension, uint16_t max_childs_number, uint16_t min_child_number>
        struct Node
        {
            typedef RTObject<ObjectT, BoundValueT, dimension> RTObjectT;
            typedef MBR<BoundValueT, dimension> MBRT;

            explicit Node(uint16_t level = 0, const MBRT &mbr = MBRT(), Node *parent = nullptr)
                    : level_(level),
                      childs_number_(0),
                      parent_(parent),
                      mbr_(mbr)
            {
                static_assert(max_childs_number >= 2 * min_child_number,
                              "Minimum number of childs should be at most a half of maximum number of childs. Otherwise split is impossible");

                if (is_leaf())
                {
                    data_ = new RTObjectT[max_childs_number];
                } else
                {
                    children_ = new Node[max_childs_number];
                }

                ++this->stats_.writes_number;
                this->stats_.bytes_written += sizeof(*this);
            }

            Node(Node const &other) = delete;

            Node &operator=(Node const &other) = delete;

            Node(Node &&other)
                    : level_(other.level_),
                      childs_number_(other.childs_number_),
                      parent_(other.parent_),
                      mbr_(other.mbr_)
            {
                if (is_leaf())
                {
                    data_ = other.data_;
                } else
                {
                    children_ = other.children_;
                    std::for_each(children_, children_ + childs_number_, [this](Node &node) { node.parent_ = this; });
                }
                other.data_ = nullptr;
                other.children_ = nullptr;
            }

            virtual ~Node()
            {
                clear();
            }

            void clear()
            {
                if (is_leaf() && data_)
                {
                    delete[] data_;
                    data_ = nullptr;
                } else if (children_)
                {
                    delete[] children_;
                    children_ = nullptr;
                }
                parent_ = nullptr;
            }

            Node &operator=(Node &&other)
            {
                clear();

                childs_number_ = other.childs_number_;
                parent_ = other.parent_;
                mbr_ = other.mbr_;
                level_ = other.level_;


                if (other.is_leaf())
                {
                    data_ = other.data_;
                } else
                {
                    children_ = other.children_;
                    std::for_each(children_, children_ + childs_number_, [this](Node &node) { node.parent_ = this; });
                }

                other.data_ = nullptr;
                other.children_ = nullptr;

                return *this;
            }

            bool is_leaf() const { return level_ == 0; }

            bool is_root() const { return parent_ == nullptr; }

            Node *get_parent() const { return parent_; }

            Node *get_root()
            {
                Node* result = this;
                while (result->parent_ != nullptr)
                {
                    result = result->parent_;
                }
                return result;
            }

            uint16_t get_level() const { return level_; }

            const MBRT &get_mbr() const { return mbr_; }

            size_t get_childs_count() const { return childs_number_; }

            Node *choose_subtree(MBRT const &mbr)
            {
                Node *result = this;

                while (!result->is_leaf())
                {
                    result = result->level_ == 1 ? result->choose_leaf(mbr) : result->choose_internal(mbr);
                }
                return result;
            }

            RTObjectT *insert(RTObjectT &&object)
            {
                if (!is_leaf() || childs_number_ >= max_childs_number)
                {
                    return nullptr;
                }

                data_[childs_number_++] = std::move(object);
                expand_mbr(data_[childs_number_ - 1].mbr_);
                return &data_[childs_number_ - 1];
            }

            Node *insert(Node &&node)
            {
                if (level_ - node.level_ != 1 || childs_number_ >= max_childs_number) { return nullptr; }

                node.parent_ = this;
                children_[childs_number_++] = std::move(node);
                expand_mbr(children_[childs_number_ - 1].mbr_);
                return &children_[childs_number_ - 1];
            }

            Node *split()
            {
                if (is_root() || childs_number_ < 2 * min_child_number)
                {
                    // if you want to split the root, you should create a new root, insert into it the current root and only then call split
                    return nullptr;
                }
                uint16_t split_axis = choose_split_axis();
                size_t separator = choose_separating_index(split_axis);

                Node second_node(level_);
                for (size_t i = separator; i < childs_number_; ++i)
                {
                    if (is_leaf())
                    {
                        second_node.insert(std::move(data_[i]));
                    } else
                    {
                        second_node.insert(std::move(children_[i]));
                    }
                }

                // cut current node.
                childs_number_ = separator;
                recalc_mbr();

                parent_->insert(std::move(second_node));
                return parent_;
            }

            Node *choose_internal(MBRT const &mbr)
            {
                if (!childs_number_) { return this; }

                return std::min_element(children_, children_ + childs_number_,
                                        std::bind(choosing_internal_less, std::placeholders::_1, std::placeholders::_2,
                                                  mbr));
            }

            Node *choose_leaf(MBRT const &mbr)
            {
                Node *result = this;

                BoundValueT min_overlap = result->mbr_.area() * childs_number_;

                for (size_t i = 0; i < childs_number_; ++i)
                {
                    MBRT expanded_child_mbr = children_[i].mbr_.expanded_mbr(mbr);
                    // sum of overlaps excluding overlap with itself
                    BoundValueT overlap_sum = std::accumulate(children_, children_ + childs_number_, BoundValueT(),
                                                              [&expanded_child_mbr](BoundValueT partial_resutlt,
                                                                                    Node &child)
                                                              {
                                                                  return partial_resutlt +
                                                                         expanded_child_mbr.overlap_area(child.mbr_);
                                                              }) - expanded_child_mbr.area();
                    if (overlap_sum < min_overlap ||
                        overlap_sum == min_overlap && choosing_internal_less(*result, children_[i], mbr))
                    {
                        result = &children_[i];
                        min_overlap = overlap_sum;
                    }
                }

                return result;
            }

            void find(const MBRT &mbr, std::vector<RTObjectT> &vector)
            {
                if (is_leaf())
                {
                    std::for_each(data_, data_ + childs_number_, [&vector, &mbr](RTObjectT &object)
                    {
                        if (object.mbr_.overlap_area(mbr) > 0)
                        {
                            vector.push_back(object);
                        }
                    });
                } else
                {
                    std::for_each(children_, children_ + childs_number_, [&vector, &mbr](Node &node)
                    {
                        if (node.mbr_.overlap_area(mbr) > 0)
                        {
                            node.find(mbr, vector);
                        }
                    });
                }
            }

        protected:
            uint16_t level_;
            size_t childs_number_;
            Node *parent_;
            union
            {
                Node *children_;
                RTObjectT *data_;
            };
            MBRT mbr_;

            struct
            {
                size_t bytes_read;
                size_t bytes_written;
                size_t reads_number;
                size_t writes_number;
            } stats_ = {0, 0, 0, 0};

            static bool choosing_internal_less(Node const &first, Node const &second, MBRT const &mbr)
            {
                BoundValueT first_expantion_area = first.mbr_.expantion_area(mbr);
                BoundValueT second_expantion_area = second.mbr_.expantion_area(mbr);
                return first_expantion_area < second_expantion_area ||
                       first_expantion_area == second_expantion_area && first.mbr_.area() < second.mbr_.area();
            }

            void expand_mbr(MBRT const &mbr)
            {
                Node *current = this;
                while (current && current->mbr_.expand(mbr))
                {
                    current = current->parent_;
                }
            }

            void sort_children_by_min_bounds(uint16_t axis)
            {
                std::sort(children_, children_ + childs_number_,
                          [axis](Node &first, Node &second)
                          {
                              return first.mbr_.min[axis] < second.mbr_.min[axis];
                          });
            }

            void sort_children_by_max_bounds(uint16_t axis)
            {
                std::sort(children_, children_ + childs_number_,
                          [axis](Node &first, Node &second)
                          {
                              return first.mbr_.max[axis] < second.mbr_.max[axis];
                          });
            }

            void sort_data_by_min_bounds(uint16_t axis)
            {
                std::sort(data_, data_ + childs_number_,
                          [axis](RTObjectT &first, RTObjectT &second)
                          {
                              return first.mbr_.min[axis] < second.mbr_.min[axis];
                          });
            }

            void sort_data_by_max_bounds(uint16_t axis)
            {
                std::sort(data_, data_ + childs_number_,
                          [axis](RTObjectT &first, RTObjectT &second)
                          {
                              return first.mbr_.max[axis] < second.mbr_.max[axis];
                          });
            }

            uint16_t choose_split_axis()
            {
                uint16_t best_axis = 0;
                BoundValueT min_margins_sum = -1;

                for (uint16_t i = 0; i < dimension; ++i)
                {
                    BoundValueT margins_sum = 0;

                    if (is_leaf())
                    {
                        sort_data_by_max_bounds(i);
                        margins_sum += calculate_distribution_margin_sum();

                        sort_data_by_min_bounds(i);
                        margins_sum += calculate_distribution_margin_sum();

                    } else
                    {
                        sort_children_by_max_bounds(i);
                        margins_sum += calculate_distribution_margin_sum();

                        sort_children_by_min_bounds(i);
                        margins_sum += calculate_distribution_margin_sum();
                    }

                    if (margins_sum < min_margins_sum || min_margins_sum < 0)
                    {
                        min_margins_sum = margins_sum;
                        best_axis = i;
                    }
                }

                return best_axis;
            }

            BoundValueT calculate_distribution_margin_sum()
            {
                BoundValueT result = BoundValueT();

                if (childs_number_ < 2 * min_child_number) { return result; }

                size_t distribution_range = childs_number_ - 2 * min_child_number;

                MBRT first_node_mbr = MBRT();
                first_node_mbr.clear();

                MBRT second_node_mbr = MBRT();
                second_node_mbr.clear();

                if (is_leaf())
                {
                    // pointers to data objects that can be distributed to any of new nodes
                    RTObjectT *unstable_data_begin = data_ + min_child_number;
                    RTObjectT *unstable_data_end = data_ + childs_number_ - min_child_number;


                    // place minimum number of data nodes to new nodes
                    first_node_mbr.expand(data_, unstable_data_begin);
                    second_node_mbr.expand(unstable_data_end, data_ + childs_number_);

                    for (size_t k = 0; k <= distribution_range; ++k)
                    {
                        result +=
                                first_node_mbr.expanded_mbr(unstable_data_begin, unstable_data_begin + k).perimeter() +
                                second_node_mbr.expanded_mbr(unstable_data_begin + k, unstable_data_end).perimeter();
                    }

                } else
                {
                    // pointers to children that can be distributed to any of new nodes
                    Node *unstable_children_begin = children_ + min_child_number;
                    Node *unstable_children_end = children_ + childs_number_ - min_child_number;


                    // place minimum number of childs to new nodes
                    first_node_mbr.expand(children_, unstable_children_begin);
                    second_node_mbr.expand(unstable_children_end, children_ + childs_number_);

                    for (size_t k = 0; k <= distribution_range; ++k)
                    {
                        result += first_node_mbr.expanded_mbr(unstable_children_begin,
                                                              unstable_children_begin + k).perimeter() +
                                  second_node_mbr.expanded_mbr(unstable_children_begin + k,
                                                               unstable_children_end).perimeter();
                    }
                }


                return result;
            }

            /**
             * @return k in [min_childs_count; childs_number - min_childs_count] that is the dilemeter of childs for first and second splitted nodes
             *
             */
            size_t choose_separating_index(uint16_t split_axis)
            {

                BoundValueT min_overlap_value = mbr_.area();
                BoundValueT current_area_value = 2 * mbr_.area();
                size_t result_index = min_child_number;

                if (is_leaf())
                {
                    sort_data_by_min_bounds(split_axis);
                } else
                {
                    sort_children_by_min_bounds(split_axis);
                }

                update_distribution_values(result_index, min_overlap_value, current_area_value);

                if (is_leaf())
                {
                    sort_data_by_max_bounds(split_axis);
                } else
                {
                    sort_children_by_max_bounds(split_axis);
                }
                bool is_result_for_sort_by_max = update_distribution_values(result_index, min_overlap_value,
                                                                            current_area_value);

                if (!is_result_for_sort_by_max)
                {
                    if (is_leaf())
                    {
                        sort_data_by_min_bounds(split_axis);
                    } else
                    {
                        sort_children_by_min_bounds(split_axis);
                    }
                }

                return result_index;
            }

            /**
             * @return true if was at least one update, and false if there were no updates
             */
            bool
            update_distribution_values(size_t &index, BoundValueT &min_overlap_value, BoundValueT &current_area_value)
            {
                bool result = false;

                if (childs_number_ < 2 * min_child_number) { return result; }

                size_t distribution_range = childs_number_ - 2 * min_child_number;

                MBRT first_node_mbr = MBRT();
                first_node_mbr.clear();

                MBRT second_node_mbr = MBRT();
                second_node_mbr.clear();

                // pointers to children that can be distributed to any of new nodes
                Node *unstable_children_begin = children_ + min_child_number;
                Node *unstable_children_end = children_ + childs_number_ - min_child_number;

                RTObjectT *unstable_data_begin = data_ + min_child_number;
                RTObjectT *unstable_data_end = data_ + childs_number_ - min_child_number;


                if (is_leaf())
                {
                    // place minimum number of data objects to new nodes
                    first_node_mbr.expand(data_, unstable_data_begin);
                    second_node_mbr.expand(unstable_data_end, data_ + childs_number_);
                } else
                {
                    // place minimum number of childs to new nodes
                    first_node_mbr.expand(children_, unstable_children_begin);
                    second_node_mbr.expand(unstable_children_end, children_ + childs_number_);
                }

                for (size_t k = 0; k <= distribution_range; ++k)
                {
                    MBRT second_extended_mbr;

                    if (is_leaf())
                    {
                        first_node_mbr.expand((unstable_data_begin + k - 1)->mbr_);
                        second_extended_mbr = second_node_mbr.expanded_mbr(unstable_data_begin + k, unstable_data_end);
                    } else
                    {
                        first_node_mbr.expand((unstable_children_begin + k - 1)->mbr_);
                        second_extended_mbr = second_node_mbr.expanded_mbr(unstable_children_begin + k,
                                                                           unstable_children_end);
                    }

                    BoundValueT overlap_value = first_node_mbr.overlap_area(second_extended_mbr);
                    BoundValueT area_value = first_node_mbr.area() + second_extended_mbr.area();

                    if (overlap_value < min_overlap_value ||
                        overlap_value == min_overlap_value && area_value < current_area_value)
                    {
                        index = min_child_number + k;
                        min_overlap_value = overlap_value;
                        current_area_value = area_value;
                        result = true;
                    }
                }

                return result;
            }

            void recalc_mbr()
            {
                mbr_.clear();
                if (is_leaf())
                {
                    mbr_.expand(data_, data_ + childs_number_);
                } else
                {
                    mbr_.expand(children_, children_ + childs_number_);
                }
            }
        };
    }

    template<class ObjectT, class BoundValueT, uint16_t dimension, uint16_t max_childs_number, uint16_t min_child_number>
    struct Rtree
    {
        typedef details::Node<ObjectT, BoundValueT, dimension, max_childs_number, min_child_number> NodeT;
        typedef details::RTObject<ObjectT, BoundValueT, dimension> RTObjectT;
        typedef details::MBR<BoundValueT, dimension> MBRT;

        Rtree() : root_(new NodeT()) {}

        virtual ~Rtree()
        {
            if (root_)
            {
                delete root_;
            }
        }

        void insert(RTObjectT &&object)
        {
            NodeT *node = root_->choose_subtree(object.mbr_);
            if (node->get_childs_count() >= max_childs_number)
            {
                split(node, object.mbr_);
                node = root_->choose_subtree(object.mbr_);
            }

            node->insert(std::move(object));
        }

        std::vector<RTObjectT> find(const MBRT &mbr)
        {
            std::vector<RTObjectT> result;
            root_->find(mbr, result);
            return result;
        }

    private:

        NodeT *root_;

        /**
         * split up to root if needed
         * @param node is a splitted node
         * @param mbr is an mbr for choosing right node after splitting
         * @return pointer to the parent node, where results of the split were inserted
         */
        NodeT *split(NodeT *node, MBRT &mbr)
        {
            NodeT *parent;
            // we must be able to insert a new node into the parent after splitting. Otherwise, the parent needs to be splited too
            if ((parent = node->get_parent()) && parent->get_childs_count() >= max_childs_number)
            {
                node = split(parent, mbr)->choose_internal(mbr);
                if (node->get_childs_count() < max_childs_number)
                {
                    return node;
                }
            }

            // if we want to split the root, then we should create a new root that will contain 2 children after splitting
            if (node == root_)
            {
                NodeT *new_root = new NodeT(root_->get_level() + 1);
                node = new_root->insert(std::move(*root_));
                delete root_;
                root_ = new_root;
            }

            return node->split()->choose_internal(mbr);
        }
    };
}
#endif //PCM_R_TREE_RTREE_H
