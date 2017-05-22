#ifndef PCM_R_TREE_MBR_H
#define PCM_R_TREE_MBR_H

namespace rtree
{
    namespace details
    {
        /**
         * MBR --- Minimal Bounding Rectangle for R-tree node which cover all MBR of children
         * @tparam BoundValueT is a type of the variable storing a coordinate of the rectangle side. Should be an arithmetic type!
         * @tparam dimension is a dimension of R-tree
         */

        template<class BoundValueT, uint16_t dimension, class SFINAE = void>
        struct MBR;

        template<class BoundValueT, uint16_t dimension>
        struct MBR<BoundValueT, dimension, typename std::enable_if<std::is_arithmetic<BoundValueT>::value>::type>
        {
            BoundValueT min[dimension];
            BoundValueT max[dimension];

            MBR()
            {
                clear();
            }

            /**
             * example: MBR with  BoundValueT=int, dimension=4 can be createad as follow:
             * {{1, 10}, {2, 4}, {-1, -1}, {0, 100}}
             * where i-th pair defines min and max bounds for i-th dimention
             */
            MBR(const std::initializer_list<std::pair<BoundValueT, BoundValueT>> &bounds_list)
            {
                uint16_t dim = 0;
                for (auto bounds : bounds_list)
                {
                    assert(bounds.first <= bounds.second);
                    min[dim] = bounds.first;
                    max[dim++] = bounds.second;
                }
            }

            void clear()
            {
                std::fill_n(min, dimension, std::numeric_limits<BoundValueT>::max());
                std::fill_n(max, dimension, std::numeric_limits<BoundValueT>::min());
            }

            BoundValueT area() const
            {
                BoundValueT result = 1;
                for (uint16_t i = 0; i < dimension; ++i)
                {
                    result *= max[i] - min[i];
                }

                return result;
            }

            BoundValueT perimeter() const
            {
                BoundValueT result = 0;
                for (uint16_t i = 0; i < dimension; ++i)
                {
                    result += max[i] - min[i];
                }

                return result * 2;
            }

            bool expand(const MBR &other)
            {
                bool result = false;
                for (uint16_t i = 0; i < dimension; ++i)
                {
                    if (other.min[i] < min[i])
                    {
                        min[i] = other.min[i];
                        result = true;
                    }
                    if (other.max[i] > max[i]) {
                        max[i] = other.max[i];
                        result = true;
                    }
                }
                return result;
            }

            template<class Iterator>
//    typename std::enable_if<std::is_same<decltype( std::declval<typename std::iterator_traits<Iterator>::value_type>().get_mbr() ), const MBR>::value, size_t>::type
            size_t
            expand(Iterator* begin, Iterator* end)
            {
                size_t count = 0;
                while (begin != end)
                {
                    expand(begin->get_mbr());
                    ++begin;
                    ++count;
                }
                return count;
            }

            MBR expanded_mbr(const MBR &mbr) const
            {
                MBR result = MBR();
                for (uint16_t i = 0; i < dimension; ++i)
                {
                    result.max[i] = std::max(max[i], mbr.max[i]);
                    result.min[i] = std::min(min[i], mbr.min[i]);
                }
                return result;
            }

            template<class Iterator>
//    typename std::enable_if<std::is_same<decltype( std::declval<typename std::iterator_traits<Iterator>::value_type>().get_mbr() ), const MBR>::value, MBR>::type
            MBR
            expanded_mbr(Iterator* begin, Iterator* end)
            {
                MBR result = MBR();
                while (begin != end)
                {
                    result.expand(begin->get_mbr());
                    ++begin;
                }
                return result;
            }

            BoundValueT expantion_area(const MBR &mbr) const
            {
                BoundValueT expanded_mbr_area = 1;
                for (uint16_t i = 0; i < dimension && expanded_mbr_area; ++i)
                {
                    expanded_mbr_area *= std::max(max[i], mbr.max[i]) - std::min(min[i], mbr.min[i]);
                }
                return expanded_mbr_area - area();
            }

            BoundValueT overlap_area(const MBR &mbr) const
            {
                BoundValueT overlap_area = 1;
                // if overlap_area < 0 than there is no overlap
                for (uint16_t i = 0; i < dimension && overlap_area > 0; ++i)
                {
                    overlap_area *= std::min(max[i], mbr.max[i]) - std::max(min[i], mbr.min[i]);
                }
                return overlap_area > 0 ? overlap_area : 0;
            }

            friend std::ostream &operator<<(std::ostream &os, const MBR &mbr)
            {
                os << "min: ";
                for (auto min_value : mbr.min)
                {
                    os << min_value << " ";
                }
                os << std::endl << "max: ";
                for (auto min_value : mbr.max)
                {
                    os << min_value << " ";
                }
                os << std::endl << "area: " << mbr.area();
                return os;
            }
        };
    }
}

#endif //PCM_R_TREE_MBR_H
