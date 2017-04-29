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
};

template <class ObjectT, class BoundValueT, uint16_t dimension, uint16_t max_childs_number, uint16_t min_child_number>
struct Rtree
{
    typedef Node<ObjectT, BoundValueT, dimension, max_childs_number, min_child_number> NodeT;

    virtual ~Rtree()
    { if (root_) {delete root_;} }

    void setRoot(NodeT *root)
    { root_ = root;}

private:
    NodeT *root_;
};

#endif //PCM_R_TREE_RTREE_H
