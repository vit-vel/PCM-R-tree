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

    MBR expanded_mbr(MBR<BoundValueT, dimension> &mbr) const {
        MBR<BoundValueT, dimension> result = MBR();
        for (uint16_t i = 0; i < dimension; ++i) {
            result.max[i] = std::max(max[i], mbr.max[i]);
            result.min[i] = std::min(min[i], mbr.min[i]);
        }
        return result;
    }

    BoundValueT expantion_area(MBR<BoundValueT, dimension> &mbr) const {
        BoundValueT expanded_mbr_area = 1;
        for (uint16_t i = 0; i < dimension && expanded_mbr_area; ++i) {
            expanded_mbr_area *= std::max(max[i], mbr.max[i]) - std::min(min[i], mbr.min[i]);
        }
        return expanded_mbr_area - area();
    }

    BoundValueT overlap_area(MBR<BoundValueT, dimension> &mbr) const {
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
    : mbr_(mbr), data_(data)
    {}

    virtual ~RTObject()
    { if (data_) {delete data_;} }

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
        if (isLeaf())
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
        if (isLeaf()) {
            delete [] data_;
            data_ = nullptr;
        } else
        {
            delete [] children_;
            children_ = nullptr;
        }
    }

    bool isLeaf() {
        return level_ == 0;
    }

    uint16_t getLevel() const {
        return level_;
    }

    const MBRT &getMbr() const {
        return mbr_;
    }

    NodeT* choose_subtree(MBRT &mbr) {
        NodeT *result = this;

        while (!result->isLeaf()) {
            result = result->level_ == 1 ? choose_leaf(mbr) : choose_internal(mbr);
        }
        return result;
    }

protected:
    uint16_t level_;
    size_t childs_number_;
    NodeT *parent_;
    union
    {
        NodeT *children_;
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
        NodeT *result = this;

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
};

template <class ObjectT, class BoundValueT, uint16_t dimension, uint16_t max_childs_number, uint16_t min_child_number>
struct Rtree
{
    typedef Node<ObjectT, BoundValueT, dimension, max_childs_number, min_child_number> NodeT;
    typedef MBR<BoundValueT, dimension> MBRT;

    Rtree(NodeT *root_ = nullptr) : root_(root_) {}

    void setRoot(NodeT *root)
    { root_ = root;}

    NodeT* choose_subtree(MBRT &mbr) {
        return root_->choose_subtree(mbr);
    }

private:

    NodeT *root_;
};

#endif //PCM_R_TREE_RTREE_H
