// general parser utilities
#pragma once
#include "../SLISC/scalar_arith.h"
#include "../SLISC/unicode.h"
#include "../SLISC/interv.h"

namespace slisc {

// skip N contiguous scope
// skip one "{...}"
// return one index after '}', return -1 if failed
// might return str.size()
// type can only be '{' for now
inline Long skip_scope(Str32_I str, Long_I ind, Long_I N = 1, Char32_I type = U'{')
{
    Long ind0 = ind;
    for (Long i = 0; i < N; ++i) {
        ind0 = expect(str, U"{", ind0);
        if (ind0 < 0)
            return -1;
        ind0 = pair_brace(str, ind0 - 1) + 1;
    }
    return ind0;
}

// Find the next "key{...}" in "str"
// find "key{...}{...}" when option = '2'
// if option = 'i', range does not include {}, if 'o', range from first character of <key> to '}'
// return the left interval boundary, output the right boundary
// return -1 if not found.
inline Long find_scope(Long_O right, Str32_I key, Str32_I str, Long_I start, Char option = 'i')
{
    Long ind0 = start, ind1;
    Long left;
    while (true) {
        ind1 = str.find(/*U"\\" +*/ key, ind0);
        if (ind1 < 0) {
            right = -1; return -1;
        }
        ind0 = ind1 + key.size();
        ind0 = expect(str, U"{", ind0);
        if (ind0 < 0) {
            ind0 = ind1 + key.size(); continue;
        }
        left = (option == 'i' ? ind0 : ind1);
        ind0 = pair_brace(str, ind0 - 1);
        if (option != '2') {
            right = option == 'i' ? ind0 - 1 : ind0;
            break;
        }
        else {
            ind0 = expect(str, U"{", ind0 + 1);
            if (ind0 < 0)
                continue;
            ind0 = pair_brace(str, ind0 - 1);
            right = ind0;
            break;
        }
    }
    return left;
}

// Find all "key{...}" in "str"
// find "key{}{}" when option = '2'
// if option = 'i', intervals are the strings inside "{}" (not including)
// if option = 'o', range from first character of "key" to '}' (including)
// return number of scopes found
inline Long find_scopes(Intvs_O intv, Str32_I key, Str32_I str, Char option = 'i')
{
    intv.clear();
    Long ind0 = 0, right;
    while (true)
    {
        ind0 = find_scope(right, key, str, ind0, option);
        if (ind0 < 0)
            return intv.size();
        intv.push(ind0, right);
        ++ind0;
    }
}

// find comments
// e.g. key = "%" for tex, "//" for c++
// key is escaped and only escaped if preceded by backslash '\'
// return number of comments found. return -1 if failed
// interval is from key[0] to '\n' (including)
inline Long find_comments(Intvs_O intv, Str32_I str, Str32_I key)
{
    Long ind0{}, ind1{};
    Long N{}; // number of comments found
    intv.clear();
    while (true) {
        ind1 = str.find(key, ind0);
        if (ind1 < 0)
            return N;
        if (ind1 == 0 || (ind1 > 0 && str.at(ind1 - 1) != U'\\')) {
            intv.pushL(ind1); ++N;
        }
        else {
            ind0 = ind1 + 1;  continue;
        }
            

        ind1 = str.find(U'\n', ind1 + 1);
        if (ind1 < 0) {// not found
            intv.pushR(str.size() - 1);
            return N;
        }
        else
            intv.pushR(ind1);
        ind0 = ind1;
    }
}

} // namespace slisc
