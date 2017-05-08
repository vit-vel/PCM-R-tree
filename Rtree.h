#ifndef PCM_R_TREE_RTREE_H
#define PCM_R_TREE_RTREE_H

#include <ostream>
#include <algorithm>
#include <functional>

/**
 * MBR --- Minimal Bounding Rectangle for R-tree node which cover all MBR of children
 * @tparam BoundValueT is a type of the variable storing a coordinate of the rectangle side
 * @tparam dimension is a dimension of R-tree
 */
template<class BoundValueT, uint16_t dimension>
struct MBR
{
    BoundValueT min[dimension];
    BoundValueT max[dimension];

    void clear()
    {
        std::fill_n(min, dimension, BoundValueT());
        std::fill_n(max, dimension, BoundValueT());
    }

    BoundValueT area() const {
        BoundValueT result = 1;
        for (uint16_t i = 0; i < dimension; ++i) {
            result *= max[i] - min[i];
        }

        return result;
    }

    BoundValueT perimeter() const {
        BoundValueT result = 0;
        for (uint16_t i = 0; i < dimension; ++i) {
            result += max[i] - min[i];
        }

        return result * 2;
    }

    bool expand(MBR &other) {
        for (uint16_t i = 0; i < dimension; ++i)
        {
            min[i] = std::min(min[i], other.min[i]);
            max[i] = std::max(max[i], other.max[i]);
        }
    }

    template<class Iterator>
    typename std::enable_if<std::is_same<typename std::iterator_traits<Iterator>::value_type,MBR>::value, size_t>::type
    expand(Iterator *begin, Iterator *end)
    {
        size_t count = 0;
        while (begin != end)
        {
            expand(*begin++);
            ++count;
        }
        return count;
    }

    MBR expanded_mbr(MBR &mbr) const {
        MBR result = MBR();
        for (uint16_t i = 0; i < dimension; ++i) {
            result.max[i] = std::max(max[i], mbr.max[i]);
            result.min[i] = std::min(min[i], mbr.min[i]);
        }
        return result;
    }

    template<class Iterator>
    typename std::enable_if<std::is_same<typename std::iterator_traits<Iterator>::value_type,MBR>::value, MBR>::type
    expanded_mbr(Iterator *begin, Iterator *end)
    {
        MBR result = MBR();
        while (begin != end)
        {
            result.expand(*begin);
            ++begin;
        }
        return result;
    }

    BoundValueT expantion_area(MBR &mbr) const {
        BoundValueT expanded_mbr_area = 1;
        for (uint16_t i = 0; i < dimension && expanded_mbr_area; ++i) {
            expanded_mbr_area *= std::max(max[i], mbr.max[i]) - std::min(min[i], mbr.min[i]);
        }
        return expanded_mbr_area - area();
    }

    BoundValueT overlap_area(MBR &mbr) const {
        BoundValueT overlap_area = 1;
        // if overlap_area < 0 than there is no overlap
        for (uint16_t i = 0; i < dimension && overlap_area > 0; ++i) {
            overlap_area *= std::min(max[i], mbr.max[i]) - std::max(min[i], mbr.min[i]);
        }
        return overlap_area > 0 ? overlap_area : 0;
    }

    friend std::ostream &operator<<(std::ostream &os, const MBR &mbr)
    {
        os << "min: ";
        for (auto min_value : mbr.min) {
            os << min_value << " ";
        }
        os << std::endl << "max: ";
        for (auto min_value : mbr.max) {
            os << min_value << " ";
        }
        os << std::endl << "area: " << mbr.area();
        return os;
    }
};

/*
 * The R-tree object that contains:
 *   1. Id of a real object or pointer to an object or db_record (pointer in this case)
 *   2. Minimum bounding rectangle of a real object
 */
template<class ObjectT, class BoundValueT, uint16_t dimension>
struct RTObject
{
    typedef MBR<BoundValueT, dimension> MBRT;

    explicit RTObject(ObjectT *data = nullptr, const MBRT &mbr = MBRT())
            : mbr_(mbr), data_(data) {}

    virtual ~RTObject() { if (data_) { delete data_; }}

    ObjectT *data_;
    MBRT mbr_;
};


template<class ObjectT, class BoundValueT, uint16_t dimension, uint16_t max_childs_number, uint16_t min_child_number>
struct Node
{
    typedef Node<ObjectT, BoundValueT, dimension, max_childs_number, min_child_number> NodeT;
    typedef RTObject<ObjectT, BoundValueT, dimension> RTObjectT;
    typedef MBR<BoundValueT, dimension> MBRT;

    explicit Node(Node *parent = nullptr, const MBRT &mbr = MBRT(), uint16_t level = 0)
            : level_(level),
              childs_number_(0),
              parent_(parent),
              mbr_(mbr)
    {
        static_assert(max_childs_number <= 2 * min_child_number, "Minimum number of childs should be at most a half of maximum number of childs. Otherwise split is impossible");
        static_assert(std::is_arithmetic<BoundValueT>::value, "Bound value type should be arithmetic");

        if (is_leaf())
        {
            data_ = new RTObjectT[max_childs_number];
        } else
        {
            children_ = new NodeT[max_childs_number];
        }

        ++this->stats_.writes_number;
        this->stats_.bytes_written += sizeof(*this);
    }

    virtual ~Node()
    {
        if (is_leaf())
        {
            delete[] data_;
            data_ = nullptr;
        } else
        {
            delete[] children_;
            children_ = nullptr;
        }
    }

    bool is_leaf()
    {
        return level_ == 0;
    }

    bool is_root()
    {
        return parent_ == nullptr;
    }

    uint16_t get_level() const {
        return level_;
    }

    const MBRT &get_mbr() const {
        return mbr_;
    }

    size_t get_childs_count() const
    {
        return childs_number_;
    }

    NodeT* choose_subtree(MBRT &mbr)
    {
        NodeT* result = this;

        while (!result->is_leaf())
        {
            result = result->level_ == 1 ? choose_leaf(mbr) : choose_internal(mbr);
        }
        return result;
    }

    bool insert(RTObjectT &object)
    {
        if (!is_leaf())
        {
            return false;
        }
        if (childs_number_ < max_childs_number)
        {
            data_[childs_number_++] = object;
            expand_mbr(object.mbr_);
        } else
        {
            split();
        }
    }

    void split()
    {
        if (childs_number_ < max_childs_number) {
            return;
        }
        size_t split_axis = choose_split_axis();
    }

protected:
    uint16_t level_;
    size_t childs_number_;
    NodeT* parent_;
    union
    {
        NodeT* children_;
        RTObjectT *data_;
    };
    MBRT mbr_;

    struct {
        size_t bytes_read;
        size_t bytes_written;
        size_t reads_number;
        size_t writes_number;
    } stats_ = {0, 0, 0, 0};

    NodeT* choose_internal(MBRT &mbr) {
        if (!childs_number_) {
            return this;
        }

        return std::min_element(children_, children_ + childs_number_,
                                std::bind(choosing_internal_less, std::placeholders::_1, std::placeholders::_2, mbr));
    }

    NodeT* choose_leaf(MBRT &mbr) {
        NodeT* result = this;

        BoundValueT min_overlap = result->mbr_.area() * childs_number_;

        for (size_t i = 0; i < childs_number_; ++i) {
            MBRT expanded_child_mbr = children_[i].mbr_.expanded_mbr(mbr);
            // sum of overlaps excluding overlap with itself
            BoundValueT overlap_sum = std::accumulate(children_, children_ + childs_number_, BoundValueT(),
                                                      [&expanded_child_mbr](BoundValueT partial_resutlt, NodeT child) {
                                                          return partial_resutlt + expanded_child_mbr.overlap_area(child.mbr_);
                                                      }) - expanded_child_mbr.area();
            if (overlap_sum < min_overlap ||
                overlap_sum == min_overlap && choosing_internal_less(*result, children_[i], mbr)) {
                result = &children_[i];
                min_overlap = overlap_sum;
            }
        }

        return result;
    }

    bool choosing_internal_less(NodeT &first, NodeT &second, MBRT &mbr) {
        BoundValueT first_expantion_area = first.mbr_.expantion_area(mbr);
        BoundValueT second_expantion_area = second.mbr_.expantion_area(mbr);
        return first_expantion_area < second_expantion_area ||
               first_expantion_area == second_expantion_area && first.mbr_.area() < second.mbr_.area();
    }

    void expand_mbr(MBRT &mbr)
    {
        NodeT* current = this;
        while (current->mbr_.expand(mbr))
        {
            current = current->parent_;
        }
    }

    uint16_t choose_split_axis()
    {
        uint16_t best_axis = 0;
        BoundValueT min_margins_sum = -1;

        for (uint16_t i = 0; i < dimension; ++i)
        {
            BoundValueT margins_sum = 0;

            // sort by maximum bounds
            std::sort(children_, children_ + childs_number_,
                      [i](NodeT &first, NodeT &second)
                      {
                          return first.mbr_.max[i] < second.mbr_.max[i];
                      });

            margins_sum += calculate_distribution_margin_sum();

            std::sort(children_, children_ + childs_number_,
                      [i](NodeT &first, NodeT &second)
                      {
                          return first.mbr_.min[i] < second.mbr_.min[i];
                      });

            margins_sum += calculate_distribution_margin_sum();

            if (margins_sum < min_margins_sum || min_margins_sum < 0) {
                min_margins_sum = margins_sum;
                best_axis = i;
            }
        }

        return best_axis;
    }

    BoundValueT calculate_distribution_margin_sum() {
        BoundValueT result = BoundValueT();
        int distribution_range = max_childs_number - 2 * min_child_number;

        if (distribution_range < 0)
        { return result; }

        MBRT first_node_mbr = MBRT();
        first_node_mbr.clear();

        MBRT second_node_mbr = MBRT();
        second_node_mbr.clear();

        // pointers to children that can be distributed to any of new nodes
        NodeT* unstable_children_begin = children_ + min_child_number;
        NodeT* unstable_children_end = children_ + childs_number_ - min_child_number;


        // place minimum number of childs to new nodes
        first_node_mbr.expand(children_, unstable_children_begin);
        second_node_mbr.expand(unstable_children_end, children_ + childs_number_);

        for (size_t k = 1; k <= distribution_range; ++k)
        {
            result += first_node_mbr.expanded_mbr(unstable_children_begin, unstable_children_begin + k).perimeter() +
                      second_node_mbr.expanded_mbr(unstable_children_begin + k, unstable_children_end).perimetr();
        }

        return result;
    }
};

template<class ObjectT, class BoundValueT, uint16_t dimension, uint16_t max_childs_number, uint16_t min_child_number>
struct Rtree
{
    typedef Node<ObjectT, BoundValueT, dimension, max_childs_number, min_child_number> NodeT;
    typedef RTObject<ObjectT, BoundValueT, dimension> RTObjectT;
    typedef MBR<BoundValueT, dimension> MBRT;

    Rtree(NodeT* root_ = nullptr) : root_(root_) {}

    void setRoot(NodeT* root) { root_ = root; }

    void insert(RTObjectT &object)
    {
        root_->choose_subtree(object.mbr_)->insert(object);
    }

private:

    NodeT* root_;
};

#endif //PCM_R_TREE_RTREE_H
