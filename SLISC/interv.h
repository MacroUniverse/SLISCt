// interval algorithms
// an interval is represented by two numbers e.g. Long interv[2]
// interv[0] is the left bound and interv[1] is the right bound
// an interval can represent, e.g. a sub vector of a vector,
// elements at the bounds are included.

#pragma once
#include "sort.h"

namespace slisc {

class Intvs;
typedef const Intvs &Intvs_I;
typedef Intvs &Intvs_O, &Intvs_IO;

class Intvs : public vector<Long>
{
public:
    typedef vector<Long> Base;

    void push(Long_I left, Long_I right);

    void pushL(Long_I i);

    void pushR(Long_I i);

    void push_back(Long_I i);

    void check_pair() const;

    Long size() const;

    const Long &L(Long_I i) const;

    Long &L(Long_I i);

    const Long &R(Long_I i) const;

    Long &R(Long_I i);

    const Long &operator[](Long_I i) const;

    Long &operator[](Long_I i);

    void erase(Long_I start, Long_I count);
};

inline void Intvs::push(Long_I left, Long_I right)
{
#ifdef SLS_CHECK_BOUNDS
    if (left < 0 || right < 0)
        SLS_ERR("must be non-negative!");
    if (isodd(Base::size()))
        SLS_ERR("last pair not finished!");
#endif
    push_back(left);
    push_back(right);
}

inline void Intvs::pushL(Long_I i)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0)
        SLS_ERR("must be non-negative!");
    if (isodd(Base::size()))
        SLS_ERR("last pair not finished!");
#endif
    push_back(i);
}

inline void Intvs::pushR(Long_I i)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0)
        SLS_ERR("must be non-negative!");
    if (!isodd(Base::size()))
        SLS_ERR("left not pushed!");
#endif
    push_back(i);
}

inline void Intvs::push_back(Long_I i)
{
    Base::push_back(i);
}

inline void Intvs::check_pair() const
{
    if (isodd(Base::size()))
        SLS_ERR("side is odd!");
}

inline Long Intvs::size() const
{
    if (isodd(Base::size()))
        SLS_ERR("last pair unfinished!");
    return Base::size() / 2;
}

inline const Long &Intvs::L(Long_I i) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i * 2 >= (Long)Base::size())
        SLS_ERR("out of bound!");
#endif
    return Base::operator[](i * 2);
}

inline Long &Intvs::L(Long_I i)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i * 2 >= (Long)Base::size())
        SLS_ERR("out of bound!");
#endif
    return Base::operator[](i * 2);
}

inline const Long &Intvs::R(Long_I i) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i * 2 + 1 >= (Long)Base::size())
        SLS_ERR("out of bound!");
#endif
    return Base::operator[](i * 2 + 1);
}

inline Long &Intvs::R(Long_I i)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i * 2 + 1 >= (Long)Base::size())
        SLS_ERR("out of bound!");
#endif
    return Base::operator[](i * 2 + 1);
}

inline const Long &Intvs::operator[](Long_I i) const
{
    SLS_ERR("use .L() or .R() instead!");
    return Base::operator[](i);
}

inline Long &Intvs::operator[](Long_I i)
{
    SLS_ERR("use .L() or .R() instead!");
    return Base::operator[](i);
}

inline void Intvs::erase(Long_I start, Long_I count)
{
    auto temp = Base::begin() + 2 * start;
    Base::erase(temp, temp + 2 * count);
}

// see if an index i falls into the scopes of ind
inline Bool is_in(Long_I i, Intvs_I intvs)
{
    for (Long j = 0; j < intvs.size(); ++j) {
        if (intvs.L(j) <= i && i <= intvs.R(j))
            return true;
    }
    return false;
}

// invert ranges in ind0, output to ind1
// [0, N-1] is the total domain
Long invert(Intvs_O ind, Intvs_I ind0, Long_I N)
{
    ind.clear();
    if (ind0.size() == 0) {
        ind.push(0, N - 1);
        return 1;
    }

    Long count{}; // total num of ranges output
    if (ind0.L(0) > 0) {
        ind.push(0, ind0.L(0) - 1); ++count;
    }
    for (Long i = 0; i < ind0.size()-1; ++i) {
        ind.push(ind0.R(i) + 1, ind0.L(i+1) - 1);
        ++count;
    }
    if (ind0.back() < N - 1) {
        ind.push(ind0.back() + 1, N - 1); ++count;
    }
    return count;
}

// combine ranges ind1 and ind2
// a range can contain another range, but not partial overlap
// return total range number, or -1 if failed.
Long combine(Intvs_O ind, Intvs_I ind1, Intvs_I ind2)
{
    Long i, N1 = ind1.size(), N2 = ind2.size();
    ind1.check_pair(); ind2.check_pair();
    if (&ind == &ind1 || &ind == &ind2) {
        SLS_ERR("aliasing is not allowed!");
    }
    if (N1 == 0) {
        ind = ind2; return N2;
    }
    else if (N2 == 0) {
        ind = ind1; return N1;
    }

    // load start and end, and sort
    vector<Long> start, end;
    for (i = 0; i < N1; ++i)
        start.push_back(ind1.L(i));
    for (i = 0; i < N1; ++i)
        end.push_back(ind1.R(i));
    for (i = 0; i < N2; ++i)
        start.push_back(ind2.L(i));
    for (i = 0; i < N2; ++i)
        end.push_back(ind2.R(i));
    sort2(start, end);

    // load ind
    ind.clear();
    ind.pushL(start[0]);
    i = 0;
    while (i < Size(start) - 1) {
        if (end[i] > start[i + 1]) {
            if (end[i] > end[i + 1]) {
                end[i + 1] = end[i]; ++i;
            }
            else {
                cout << "error! range overlap!" << endl;
                return -1;  // break point here
            }
        }
        else if (end[i] == start[i + 1] - 1)
            ++i;
        else {
            ind.pushR(end[i]);
            ind.pushL(start[i + 1]);
            ++i;
        }
    }
    ind.pushR(end.back());
    return ind.size();
}

// combine intervals ind and ind1, and assign the result to ind
Long combine(Intvs_O ind, Intvs_I ind1)
{
    Intvs temp;
    Long ret = combine(temp, ind, ind1);
    ind = temp;
    return ret;
}

} // namespace slisc
