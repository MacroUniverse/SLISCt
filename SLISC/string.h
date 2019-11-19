// string utilities

#pragma once
#include "global.h"

namespace slisc {

// string utilities

template <typename T>
inline void num2str(Str_O str, T s)
{
    str = to_string(s);
    // erase trailing zeros
    if (str.find('.') != Str::npos)
        str.erase(str.find_last_not_of('0') + 1);
}

template <typename T>
inline Str num2str(T s)
{
    Str str;
    num2str(str, s);
    return str;
}

// trim all occurance of key on the left
// return the number of charaters trimed
// e.g. key = "\n " to trim space and '\n'
inline Long trimL(Str_IO str, Str_I key = " ")
{
    Long N;
    Long ind = str.find_first_not_of(key);
    if (ind < 0) {
        N = str.size();
        str.clear();
        return N;
    }
    N = ind;
    str.erase(0, N);
    return N;
}

// trim all occurance of key on the right
// return the number of charaters trimed
// e.g. key = "\n " to trim space and '\n'
inline Long trimR(Str_IO str, Str_I key = " ")
{
    Long N;
    Long ind = str.find_last_not_of(key);
    if (ind < 0) {
        N = str.size();
        str.clear();
        return N;
    }
    str.erase(ind + 1);
    N = str.size() - ind;
    return N;
}

// trim both sides
// e.g. key = "\n " to trim space and '\n'
inline Long trim(Str_IO str, Str_I key = " ")
{
    return trimL(str, key) + trimR(str, key);
}

// replace every "\r\n" with "\n"
inline Long CRLF_to_LF(Str_IO str)
{
    Long ind0{}, N{};
    while (true) {
        ind0 = str.find("\r\n", ind0);
        if (ind0 < 0) return N;
        str.erase(ind0, 1);
    }
}

} // namespace slisc
