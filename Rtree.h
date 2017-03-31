#ifndef PCM_R_TREE_RTREE_H
#define PCM_R_TREE_RTREE_H

#include <ostream>

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
        return os;
    }
};

template<class ObjectT, class BoundValueT, uint16_t dimension>
struct Node {
    typedef MBR<BoundValueT, dimension> MBRT;

    Node() {}
    Node(const MBRT &mbr) : mbr_(mbr)
    {}

protected:
    MBRT mbr_;

    struct NodeStatistics {
        size_t bytes_read;
        size_t bytes_written;
        size_t reads_number;
        size_t writes_number;
    } stats_;
};

template <class ObjectT, class BoundValueT, uint16_t dimension>
struct LeafNode : Node<ObjectT, BoundValueT, dimension>
{
    typedef Node<ObjectT, BoundValueT, dimension> NodeT;

    explicit LeafNode(ObjectT *data = nullptr, const typename NodeT::MBRT &mbr = typename NodeT::MBRT())
            : NodeT(mbr), data_(data)
    {}

    virtual ~LeafNode()
    { if (data_) {delete data_;} }

private:
    ObjectT *data_;
};

template<class ObjectT, class BoundValueT, uint16_t dimension>
struct IntrenalNode : Node<ObjectT, BoundValueT, dimension>
{
    typedef Node<ObjectT, BoundValueT, dimension> NodeT;
    typedef LeafNode<ObjectT, BoundValueT, dimension> LeafNodeT;

    IntrenalNode(uint16_t max_childs_number, uint16_t min_child_number, const typename NodeT::MBRT &mbr = typename NodeT::MBRT())
            : NodeT(mbr),
              childs_number_(0),
              min_childs_number_(min_child_number),
              max_childs_number_(max_childs_number),
              children_(new LeafNodeT[max_childs_number])
    {
        ++this->stats_.writes_number;
        this->stats_.bytes_written += sizeof(*this);
    }

    virtual ~IntrenalNode()
    { clear(); }

    virtual void clear()
    {
        if (children_) {
            this->mbr_.clear();
            delete [] children_;
            children_ = nullptr;
            childs_number_ = 0;
        }
    }

private:
    size_t childs_number_;
    size_t min_childs_number_;
    size_t max_childs_number_;
    NodeT *children_;
};

template <class ObjectT, class BoundValueT, uint16_t dimension>
struct Rtree
{
    typedef Node<ObjectT, BoundValueT, dimension> NodeT;

    virtual ~Rtree()
    { if (root_) {delete root_;} }

    void setRoot(NodeT *root)
    { root_ = root;}

private:
    NodeT *root_;
};

#endif //PCM_R_TREE_RTREE_H
