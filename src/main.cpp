#include <iostream>
#include "Rtree.h"
#include "PCMRtree.h"
#include "smoke_test.h"

using bounds_type = int32_t;
using value_type = int32_t;

constexpr uint16_t dimension = 2;
constexpr size_t min_elements = 2;
constexpr size_t max_elements = 4;

rtree::details::RTObject<value_type, bounds_type, dimension> random_rtobject()
{
    rtree::details::RTObject<value_type, bounds_type, dimension> object;
    for (size_t j = 0; j < dimension; ++j)
    {
        int min = random() % INT16_MAX;
        int max = random() % INT16_MAX;
        if (min > max)
        {std::swap(min, max);}

        object.mbr_.min[j] = min;
        object.mbr_.max[j] = max;
    }

    return std::move(object);
}

int main(int argc, char *argv[])
{

    unsigned long long elements_number = 10000;
    bool pcm = true;

    if (argc > 1)
    {
        elements_number = std::stoull(argv[1]);
    }
    if (argc > 2)
    {
        pcm = std::stoi(argv[2]) == 1;
    }

    if (pcm)
    {
        rtree::pcm::Rtree<value_type, bounds_type, dimension, max_elements, min_elements> tree;

        for (size_t i = 0; i < elements_number; ++i)
        {
            tree.insert(random_rtobject());
        }
    } else
    {
        rtree::Rtree<value_type, bounds_type, dimension, max_elements, min_elements> tree;

        for (size_t i = 0; i < elements_number; ++i)
        {
            tree.insert(random_rtobject());
        }
    }

//    run_tests();

    return 0;
}
