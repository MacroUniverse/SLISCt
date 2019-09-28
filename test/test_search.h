#pragma once
#include "../SLISC/search.h"
#include "../SLISC/slice_arith.h"

void test_search()
{
    { // test search_row
        using namespace slisc;
        Long N1 = 10, N2 = 10;
        CmatDoub a(N1, N2); rand(a);
        DvecDoub sli;
        for (Long i = 0; i < N1; ++i) {
            slice2(sli, a, i);
            if (search_row(sli, a) != i)
                SLS_ERR("failed!");
        }
    }
}
