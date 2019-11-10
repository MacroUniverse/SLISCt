#pragma once
#include "global.h"

namespace slisc {

// see if elm == vec[i], return i
// return -1 if not found
template <class T, class Tv>
inline Long search(const T &elm, const Tv &vec, Long start = 0)
{
    for (Long i = start; i < Size(vec); ++i) {
        if (elm == vec[i])
            return i;
    }
    return -1;
}

// search a row from a matrix
// starting from the row `start`
// return the row index, return -1 if not found
template <class Tv, class Tmat>
inline Long search_row(const Tv &v, const Tmat &a, Long_I start = 0)
{
#ifdef SLS_CHECK_SHAPE
    if (v.size() != a.n2())
        SLS_ERR("wrong shape");
#endif
    for (Long i = start; i < a.n1(); ++i) {
        if (v == slice2(a, i))
            return i;
    }
    return -1;
}

// check if v[i] == v[j] for any i and j
// return the index of one of the repeated elements
// return -1 if no reapeat
template <class Tv>
inline Long find_repeat(Tv v)
{
    for (Long i = 0; i < Size(v); ++i) {
        for (Long j = i + 1; j < Size(v); ++j) {
            if (v[i] == v[j])
                return i;
        }
    }
    return -1;
}

// return 0, output ind: if v[ind] == s
// return -1: if s is not found and s < v[0]
// return 1: if s is not found and v[end] < s
// return -2: if s is not found and v[0] < s < v[end], output ind so that v[ind] < s < v[ind+1]
// return -3: if s is not found because v.size() == 0
// operator-(v[i], s) must be implemented and only need to return the correct sign or 0
template <class T, class T1>
inline Int lookupInt(Long_O ind, const T1 &v, const T &s)
{
    Int diff;
    Long i, N = v.size(), ind1 = 0, ind2 = N - 1;
    if (N < 1) return -3;
    diff = v[0] - s;
    if (diff > 0) return -1;
    if (diff == 0) {
        ind = 0; return 0;
    }
    diff = v[ind2] - s;
    if (diff < 0) return 1;
    
    // N >= 2 from here
    if (diff == 0) {
        ind = ind2; return 0;
    }
    // N > 2 from here
    for (i = 0; i < N; ++i) {
        ind = (ind1 + ind2) / 2;
        diff = v[ind] - s;
        if (diff == 0) return 0;
        if (diff > 0) ind2 = ind;
        else ind1 = ind;
        if (ind2 - ind1 == 1) {
            ind = ind1; return -2;
        }
    }
    SLS_ERR("slisc::lookup(): unknown error!");
    return 0;
}
} // namespace slisc
