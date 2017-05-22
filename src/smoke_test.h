#ifndef PCM_R_TREE_SMOKE_TEST_H
#define PCM_R_TREE_SMOKE_TEST_H

#include <iostream>
#include <assert.h>
#include "Rtree.h"
#include "PCMRtree.h"

using namespace rtree::details;

void test_mbr_creation()
{
    MBR<int, 3> mbr1{{1, 4},{2, 5}, {3, 6}};
    MBR<int, 3> mbr2{{11, 44}, {22, 55}, {33, 66}};

//    std::cout << mbr1 << std::endl;
//    std::cout << mbr2 << std::endl;

    MBR<int, 3> mbr3(mbr1);
    mbr2 = mbr1;
    mbr1.max[0] = 100000;
    assert(mbr1.max[0] != mbr2.max[0]);
    assert(mbr1.max[0] != mbr3.max[0]);
    assert(mbr2.max[0] == mbr3.max[0]);

//    std::cout << mbr1 << std::endl;
//    std::cout << mbr2 << std::endl;
//    std::cout << mbr3 << std::endl;

//    MBR<int, 30> zero_test;
//    std::cout << zero_test << std::endl;

    MBR<int, 2> q;
    MBR<size_t, 2> w;
    MBR<float, 2> e;
    MBR<double, 2> r;
    MBR<char, 2> t;
    MBR<bool, 2> y;
    MBR<long, 2> u;
    MBR<unsigned long, 2> i;
//    MBR<std::string, 2> k; // shouldn't compile

}

void test_mbr_expantion_and_overlapping()
{
    MBR<float, 2> mbr = {{1,5},{1,4}};
    assert(mbr.area() == 12);
    assert(mbr.perimeter() == 14);

    MBR<float, 2> crossing_mbr = {{4, 6}, {3, 5}};
    assert(mbr.expantion_area(crossing_mbr) == 8);
    assert(mbr.overlap_area(crossing_mbr) == 1);
    assert(mbr.expantion_area(crossing_mbr) == (mbr.expanded_mbr(crossing_mbr).area() -
                                                mbr.area()));


    MBR<float, 2> inner_mbr = {{2, 4}, {1, 3}};
    assert(mbr.expantion_area(inner_mbr) == 0);
    assert(mbr.overlap_area(inner_mbr) == inner_mbr.area());
    assert(mbr.expantion_area(inner_mbr) == (mbr.expanded_mbr(inner_mbr).area() - mbr.area()));

    MBR<float, 2> outer_mbr = {{-1, 0}, {0, 3}};
    assert(mbr.expantion_area(outer_mbr) == 12);
    assert(mbr.overlap_area(outer_mbr) == 0);
    assert(mbr.expantion_area(outer_mbr) == (mbr.expanded_mbr(outer_mbr).area() - mbr.area()));
}

void test_find_object()
{
    rtree::Rtree<int, float, 2, 2, 1> tree;
    tree.insert(RTObject<int, float, 2>(nullptr, {{2, 3}, {5, 6}}));
    tree.insert(RTObject<int, float, 2>(nullptr, {{2, 3}, {2, 6}}));
    tree.insert(RTObject<int, float, 2>(nullptr, {{1, 6}, {4, 8}}));
    tree.insert(RTObject<int, float, 2>(nullptr, {{3, 5}, {7, 10}}));
    tree.insert(RTObject<int, float, 2>(nullptr, {{4, 7}, {6, 9}}));

    assert(tree.find({{0, 3}, {0, 3}}).size() == 1);
    assert(tree.find({{4, 5}, {7, 8}}).size() == 3);
    assert(tree.find({{3, 4}, {7, 8}}).size() == 2);
    assert(tree.find({{2, 3}, {7, 8}}).size() == 1);
    assert(tree.find({{3, 4}, {8, 9}}).size() == 1);
    assert(tree.find({{1, 6}, {4, 8}}).size() == 5);
    assert(tree.find({{1, 2}, {3, 4}}).size() == 0);
}

void run_tests()
{
    test_mbr_creation();
    test_mbr_expantion_and_overlapping();
    test_find_object();
}

#endif //PCM_R_TREE_SMOKE_TEST_H
