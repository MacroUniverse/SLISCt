#pragma once
#include "scalar_arith.h"
#include <sstream>
#include <fstream>
#include <locale>
#include <codecvt>
#ifdef _MSC_VER
#include <Windows.h> // for console unicode output
#undef max
#undef min
#endif
#include "string.h"
#include "file.h"
#include "utfcpp/utf8.h"

#ifdef _MSC_VER
#define SLS_USE_UTFCPP
#endif

namespace slisc {

using std::stringstream;

inline Long CRLF_to_LF(Str32_IO str);

// set windows console to display utf-8
#ifdef _MSC_VER
struct set_windows_console_utf8 {
    set_windows_console_utf8() {
        SetConsoleOutputCP(65001); // code page 65001 is UTF-8
    }
};
// in case of ODR error, put this in main function;
set_windows_console_utf8 yes_set_windows_console_utf8;
#endif

// write Str to file
inline void write_file(Str_I str, Str_I name)
{
    ofstream fout(name, std::ios::binary);
    fout << str;
    fout.close();
}

#ifdef SLS_USE_UTFCPP
// convert from UTF-8 Str to UTF-32 Str32
inline void utf8to32(Str32_O str32, Str_I str)
{
    str32.clear();
    utf8::utf8to32(str.begin(), str.end(), back_inserter(str32));
}

// convert from UTF-32 Str32 to UTF-8 Str
inline void utf32to8(Str_O str, Str32_I str32)
{
    str.clear();
    utf8::utf32to8(str32.begin(), str32.end(), back_inserter(str));
}
#else
// convert from UTF-8 Str to UTF-32 Str32
inline void utf8to32(Str32_O str32, Str_I str)
{
    std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> myconv;
    str32 = myconv.from_bytes(str);
}

// convert from UTF-32 Str32 to UTF-8 Str
inline void utf32to8(Str_O str, Str32_I str32)
{
    std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> myconv;
    str = myconv.to_bytes(str32);
}
#endif

// convert from UTF-8 Str32 to UTF-32 Str
inline Str32 utf8to32(Str_I str)
{
    Str32 str32;
    utf8to32(str32, str);
    return str32;
}

// convert from UTF-32 Str32 to UTF-8 Str
inline Str utf32to8(Str32_I str32)
{
    Str str;
    utf32to8(str, str32);
    return str;
}

// display Str32
inline std::ostream &operator<<(std::ostream &out, Str32_I str32)
{
    Str str;
    utf32to8(str, str32);
    out << str;
    return out;
}

// operator+ that converts Str to Str32
inline Str32 operator+(Str32_I str32, Str_I str)
{
    return str32 + utf8to32(str);
}

inline Str32 operator+(Str_I str, Str32_I str32)
{
    return utf8to32(str) + str32;
}

// not case sensitive on Windows, see file_exist_case()
inline Bool file_exist(Str32_I fname) {
    return file_exist(utf32to8(fname));
}

// read a line from a string, from str[start] to 1 char before '\n'
// return the start of the next line, return -1 if out of bound
// if the file ends with `\n`, then the line.back() is not empty
inline Long get_line(Str32_O line, Str32_I str, Long_I start)
{
    Long ind = str.find(U'\n', start);
    line = str.substr(start, ind - start);
    if (ind < 0 || ind == Size(str) - 1)
        return -1;
    return ind + 1;
}

// skip to the next line
// return the index after `\n`
// return -1 if `\n` not found
inline Long skip_line(Str32_I &str, Long_I start)
{
    Long ind = str.find(U'\n', start);
    if (ind < 0 || ind == Size(str) - 1)
        return -1;
    return ind + 1;
}

// read a UTF-8 file into UTF-32 Str32
inline void read_file(Str32_O str32, Str_I fname)
{
    Str str;
    read_file(str, fname);
    utf8to32(str32, str);
}

inline void read_file(Str32_O str32, Str32_I fname)
{
    read_file(str32, utf32to8(fname));
}

// write UTF-32 Str32 into a UTF-8 file
inline void write_file(Str32_I str32, Str_I fname)
{
    Str str;
    utf32to8(str, str32);
    write_file(str, fname);
}

inline void write_file(Str32_I str32, Str32_I fname)
{
    write_file(str32, utf32to8(fname));
}

// write a vector of strings to file
// no `\n` allowed in each string
// file will be ended by a return
inline void write_vec_str(vector_I<Str32> vec_str, Str32_I fname)
{
    Str32 str;
    for (Long i = 0; i < Size(vec_str); ++i) {
        str += vec_str[i] + U'\n';
    }
    write_file(str, fname);
}

// read the file written by `write_vec_str()`
// file must be ended by a return
inline void read_vec_str(vector_O<Str32> vec_str, Str32_I fname)
{
    Str32 str;
    vec_str.clear();
    read_file(str, fname);
    CRLF_to_LF(str);
    Long ind0 = 0;
    for (Long i = 0; ; ++i) {
        vec_str.emplace_back();
        ind0 = get_line(vec_str[i], str, ind0);
        if (ind0 < 0)
            return;
    }
}

template <class T>
inline void num2str(Str32_O str, const T &num)
{
    Str str0;
    num2str(str0, num);
    utf8to32(str, str0);
}

template <class T>
inline Str32 num2str32(const T &num)
{
    Str32 str;
    num2str(str, num);
    return str;
}

 // same as str.insert(), but return one index after key after insertion
inline Long insert(Str32_IO str, Str32_I key, Long start)
{
    str.insert(start, key);
    return start + key.size();
}

// replace every U"\r\n" with U"\n"
inline Long CRLF_to_LF(Str32_IO str)
{
    Long ind0{}, N{};
    while (true) {
        ind0 = str.find(U"\r\n", ind0);
        if (ind0 < 0) return N;
        str.erase(ind0, 1);
    }
}

// Find the next appearance of one of "key"
// output the ikey of key[ikey] found
// return the first index of key[ikey] found, return -1 if nothing found
inline Long find(Long_O ikey, Str32_I str, vector_I<Str32> keys, Long_I start)
{
    Long i{}, ind0{}, Nkey{}, imin;
    Nkey = keys.size();
    imin = str.size();
    for (i = 0; i < Nkey; ++i) {
         ind0 = str.find(keys[i], start);
         if (ind0 >= start && ind0 < imin) {
             imin = ind0; ikey = i;
         }
    }
    if (imin == Size(str)) imin = -1;
    return imin;
}

// Find the previous appearance of one of "key"
// output the ikey of key[ikey] found
// return the first index of key[ikey] found, return -1 if nothing found
// keyword will be found even if starting from the middle of it
inline Long rfind(Long_O ikey, Str32_I str, vector_I<Str32> key, Long_I start)
{
    Long i{}, ind0{}, Nkey{}, imax;
    Nkey = key.size();
    imax = -1;
    for (i = 0; i < Nkey; ++i) {
        ind0 = str.rfind(key[i], start);
        if (ind0 > imax && ind0 <= start) {
            imax = ind0; ikey = i;
        }
    }
    return imax;
}

// same as FindMultipleReverse, but able to deal with multiple match
// return the number of matches, return -1 if not found
inline Long rfind(vector_O<Long> ikey, Str32_I str, vector_I<Str32> key, Long_I start)
{
    Long i{}, ind0{}, Nkey{}, imax;
    Nkey = key.size();
    imax = -1;
    for (i = 0; i < Nkey; ++i) {
        ind0 = str.rfind(key[i], start);
        if (ind0 > imax && ind0 <= start) {
            imax = ind0; ikey.clear();
            ikey.push_back(i);
        }
        // found another
        else if (ind0 == imax)
            ikey.push_back(i);
    }
    return imax;
}

// see if a key appears followed only by only white space or '\n'
// return the index after the key found, return -1 if nothing found.
inline Long expect(Str32_I str, Str32_I key, Long_I start)
{
    Long ind = start;
    Long ind0 = 0;
    Long L = str.size();
    Long L0 = key.size();
    Char32 c0, c;
    while (true) {
         c0 = key.at(ind0);
         c = str.at(ind);
         if (c == c0) {
             ++ind0;
             if (ind0 == L0)
                 return ind + 1;
         }
         else if (c != U' ' && c != U'\n')
             return -1;
         ++ind;
         if (ind == L)
             return -1;
    }
}

// trim all occurance of key on the left
// return the number of charaters trimed
// e.g. key = "\n " to trim space and '\n'
inline Long trimL(Str32_IO str, Str32_I key = U" ")
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
inline Long trimR(Str32_IO str, Str32_I key = U" ")
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
inline Long trim(Str32_IO str, Str32_I key = U" ")
{
    return trimL(str, key) + trimR(str, key);
}

// check if a character is a letter (a-z, A-Z)
inline Bool is_letter(Char32_I c)
{
    if ((c >= U'a' && c <= U'z') || (c >= U'A' && c <= U'Z'))
        return true;
    return false;
}

inline Bool is_num(Char32_I c)
{
    if (c >= U'0' && c <= U'9')
        return true;
    return false;
}

// check if a character is alphanumeric (a-z, A-Z, 0-9)
inline Bool is_alphanum(Char32_I c)
{
    if (is_letter(c) || is_num(c))
        return true;
    return false;
}

// check if a character is alphanumeric or an underscore
inline Bool is_alphanum_(Char32_I c)
{
    if (is_alphanum(c) || c == U'_')
        return true;
    return false;
}

// check if a word is a whole word
inline Bool is_whole_word(Str32_I str, Long_I ind, Long_I count)
{
    // check left
    Long ind0 = ind - 1;
    if (ind0 >= 0 && is_alphanum_(str[ind0]))
        return false;
    // check right
    ind0 += 1 + count;
    if (ind0 < Size(str) && is_alphanum_(str[ind0]))
        return false;
    return true;
}

// find whole word, like in Visual Studio Code, begin from "str[start]"
// return the first index of key found, return -1 if not found
inline Long find_whole_word(Str32_I str, Str32_I key, Long_I start)
{
    Long ind0 = start;
    while (true) {
        ind0 = str.find(key, ind0);
        if (ind0 < 0)
            return -1;
        if (is_whole_word(str, ind0, key.size()))
            return ind0;
        ++ind0;
    }
}

// replace all occurance of "key" with "new_str"
// return the number of keys replaced
inline Long replace(Str32_IO str, Str32_I key, Str32_I new_str)
{
    Long ind0 = 0;
    Long Nkey = key.size();
    Long N = 0;
    while (true) {
         ind0 = str.find(key, ind0);
        if (ind0 < 0) break;
        str.erase(ind0, Nkey);
        str.insert(ind0, new_str);
        ++N; ind0 += new_str.size();
    }
    return N;
}

// find next character that is not a space
// output single character c, return the position of the c
// return -1 if not found
// TODO: replace this with basic_string::find_first_not_of
inline Long NextNoSpace(Str32_O c, Str32_I str, Long start)
{
    for (Long i = start; i < Size(str); ++i) {
        c = str.at(i);
        if (c == U" ")
            continue;
        else
            return i;
    }
    return -1;
}

// reverse version of Expect key
// return the previous index of the key found, return -2 if nothing found.
inline Long ExpectKeyReverse(Str32_I str, Str32_I key, Long start)
{
    if (start < 0)
        return -2;
    Long ind = start;
    Long L0 = key.size();
    Long ind0 = L0 - 1;
    Char32 c0, c;
    while (true) {
        c0 = key.at(ind0);
        c = str.at(ind);
        if (c == c0) {
            --ind0;
            if (ind0 < 0)
                return ind - 1;
        }
        else if (c != ' ' && c != '\n')
            return -2;
        --ind;
        if (ind < 0)
            return -2;
    }
}

// Find the next number
// return -1 if not found
inline Long find_num(Str32_I str, Long start)
{
    Long i{}, end = str.size() - 1;
    Char32 c;
    for (i = start; i <= end; ++i) {
        c = str.at(i);
        if (c >= '0' && c <= '9')
            return i;
    }
    return -1;
}

// get non-negative integer from string
// return the index after the last digit, return -1 if failed
// str[start] must be a number
inline Long str2int(Long_O num, Str32_I str, Long start)
{
    Long i{};
    Char32 c;
    c = str.at(start);
    if (c < '0' || c > '9') {
        SLS_ERR("not a number!"); return -1;  // break point here
    }
    num = c - '0';
    for (i = start + 1; i < Size(str); ++i) {
        c = str.at(i);
        if (c >= '0' && c <= '9')
            num = 10 * num + (Long)(c - '0');
        else
            return i;
    }
    return i;
}

// get non-negative double from string
// return the index after the last digit, return -1 if failed
// str[start] must be a number
inline Long str2double(Doub& num, Str32_I str, Long start)
{
    Long ind0{}, num1{}, num2{};
    ind0 = str2int(num1, str, start);
    if (str.at(ind0) != '.') {
        num = (Doub)num1;
        return ind0;
    }
    ind0 = str2int(num2, str, ind0 + 1);
    num = (Doub)num2;
    while (num >= 1)
        num /= 10;
    num += (Doub)num1;
    return ind0;
}

// delete any following ' ' or '\n' characters starting from "start" (including "start")
// return the number of characters deleted
inline Long eatR(Str32& str, Long start, Str32_I chars)
{
    Long N;
    Long ind0 = str.find_first_not_of(chars, start);
    if (ind0 < 0) {
        N = str.size() - start;
        if (N > 0)
            str.erase(start, N);
    }
    N = ind0 - start;
    if (N > 0)
        str.erase(start, N);
    return N;
}

// eat to left, see eatR
inline Long eatL(Str32& str, Long start, Str32_I chars)
{
    Long N;
    Long ind0 = str.find_last_not_of(chars, start);
    if (ind0 < 0) {
        N = start + 1;
        str.erase(0, N);
    }
    N = start - ind0;
    if (N > 0)
        str.erase(ind0 + 1, N);
    return N;
}

// Pair right brace to left one (default)
// or () or [] or anying single character
// ind is inddex of left brace
// return index of right brace, -1 if failed
inline Long pair_brace(Str32_I str, Long ind, Char32_I type = U'{')
{
    Char32 left, right;
    if (type == U'{' || type == U'}') {
        left = U'{'; right = U'}';
    }
    else if (type == U'(' || type == U')') {
        left = U'('; right = U')';
    }
    else if (type == U'[' || type == U']') {
        left = U'['; right = U']';
    }
    else {// anything else
        left = type; right = type;
    }

    Char32 c, c0 = ' ';
    Long Nleft = 1;
    for (Long i = ind + 1; i < Size(str); i++)
    {
        c = str.at(i);
        if (c == left && c0 != '\\')
            ++Nleft;
        else if (c == right && c0 != '\\')
        {
            --Nleft;
            if (Nleft == 0)
                return i;
        }
        c0 = c;
    }
    return -1;
}

// match braces
// return -1 means failure, otherwise return number of {} paired
// output ind_left, ind_right, ind_RmatchL
inline Long MatchBraces(vector_O<Long> ind_left, vector_O<Long> ind_right,
    vector_O<Long> ind_RmatchL, Str32_I str, Long start, Long end)
{
    ind_left.resize(0); ind_right.resize(0); ind_RmatchL.resize(0);
    Char32 c, c_last = ' ';
    Long Nleft = 0, Nright = 0;
    vector<Bool> Lmatched;
    Bool matched;
    for (Long i = start; i <= end; ++i)
    {
        c = str.at(i);
        if (c == '{' && c_last != '\\')
        {
            ++Nleft;
            ind_left.push_back(i);
            Lmatched.push_back(false);
        }
        else if (c == '}' && c_last != '\\')
        {
            ++Nright;
            ind_right.push_back(i);
            matched = false;
            for (Long j = Nleft - 1; j >= 0; --j)
                if (!Lmatched[j])
                {
                    ind_RmatchL.push_back(j);
                    Lmatched[j] = true;
                    matched = true;
                    break;
                }
            if (!matched)
                return -1; // unbalanced braces
        }
        c_last = c;
    }
    if (Nleft != Nright)
        return -1;
    return Nleft;
}

// list all files in current directory, with a given extension
inline void file_list_ext(vector_O<Str32> fnames, Str32_I path, Str32_I ext, Bool_I keep_ext = true)
{
    vector<Str> fnames8;
    fnames.resize(0);
    file_list_ext(fnames8, utf32to8(path), utf32to8(ext), keep_ext);
    for (Long i = 0; i < Size(fnames8); ++i)
        fnames.push_back(utf8to32(fnames8[i]));
}

inline void file_copy(Str32_I fname_out, Str32_I fname_in, Bool_I replace)
{
    file_copy(utf32to8(fname_out), utf32to8(fname_in), replace);
}

// check if is a chinese character
// does not include punctuations
// reference: https://stackoverflow.com/questions/1366068/whats-the-complete-range-for-chinese-characters-in-unicode
inline Bool is_chinese(Char32_I c)
{
    if ((c >= U'\u2E80' && c <= U'\u2FD5') ||
        (c >= U'\u3190' && c <= U'\u319f') ||
        (c >= U'\u3400' && c <= U'\u4DBF') ||
        (c >= U'\u4E00' && c <= U'\u9FCC') ||
        (c >= U'\uF900' && c <= U'\uFAAD'))
        return true;
    return false;
}

} // namespace slisc
