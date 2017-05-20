#ifndef PCM_R_TREE_RTOBJECT_H
#define PCM_R_TREE_RTOBJECT_H

namespace rtree
{
    namespace details
    {
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

            RTObject(const RTObject &other) : RTObject(other.data_, other.mbr_)
            {}

            RTObject(RTObject && other)
                    : mbr_(other.mbr_), data_(other.data_)
            {
                other.data_ = nullptr;
            }

            RTObject& operator= (RTObject && other)
            {
                data_ = other.data_;
                mbr_ = other.mbr_;
                other.data_ = nullptr;
            }

            const MBRT& get_mbr()
            { return mbr_; }

            ObjectT *data_;
            MBRT mbr_;
        };
    }
}

#endif //PCM_R_TREE_RTOBJECT_H
