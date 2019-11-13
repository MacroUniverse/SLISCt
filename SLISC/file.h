#pragma once
#include "scalar_arith.h"
#include <sstream>
#include <fstream>
#include <codecvt>
#ifdef SLS_HAS_FILESYSTEM
#include <filesystem>
#endif

namespace slisc {

using std::stringstream;

inline void file_list(vector_O<Str> names, Str_I path);

#ifdef SLS_HAS_FILESYSTEM
// check if a file exist on Windws (case sensitive)
inline Bool file_exist_case(Str_I fname)
{
    Long ind = max((Long)fname.rfind('/'), (Long)fname.rfind('\\'));
    Str path, name;
    if (ind < 0) {
        path = "./";
        name = fname;
    }
    else {
        path = fname.substr(0, ind + 1);
        name = fname.substr(ind + 1);
    }
    vector<Str> names;
    file_list(names, path);
    if (is_in(name, names))
        return true;
    else
        return false;
}
#endif

// check if a file exist, default is case sensitive
inline Bool file_exist(Str_I fname, Bool_I case_sens = true) {
#ifndef SLS_HAS_FILESYSTEM
    ifstream f(fname);
    return f.good();
#else
    if (case_sens)
        return file_exist_case(fname);
    else {
        ifstream f(fname);
        return f.good();
    }
#endif
}

// read whole file to Str
inline void read_file(Str_O str, Str_I fname)
{
    if (!file_exist(fname))
        SLS_ERR("file \"" + fname + "\" does not exist!");
    ifstream fin(fname, std::ios::binary);
    stringstream ss;
    ss << fin.rdbuf();
    str = ss.str();
}

inline void file_rm(Str_I wildcard_name) {
    system(("rm " + wildcard_name).c_str());
}

// list all files in current directory
// only works for linux
#ifdef __GNUC__
inline void file_list(vector_O<Str> fnames, Str_I path)
{
    Str temp_fname, temp_fname_pref = "SLISC_temporary_file";

    // create unused temporary file name
    for (Long i = 0; i < 1000; ++i) {
        if (i == 999)
            SLS_ERR("too many temporary files!");
        temp_fname = temp_fname_pref + to_string(i);
        if (!file_exist(temp_fname)) break;
    }
    
    // save a list of all files (no folder) to temporary file
    Int ret = system(("ls -p " + path + " | grep -v / > " + temp_fname).c_str()); ++ret;

    // read the temporary file
    ifstream fin(temp_fname);
    for (Long i = 0; i < 10000; ++i) {
        Str name;
        std::getline(fin, name);
        if (fin.eof())
            break;
        fnames.push_back(name);
    }
    fin.close();

    // remove temporary file
    std::remove(temp_fname.c_str());
}
#else
#ifdef SLS_HAS_FILESYSTEM
// std::filesystem implementation of file_list()
// works in Visual Studio, not gcc 8
// directory example: "C:/Users/addis/Documents/GitHub/SLISC/"
inline void file_list(vector_O<Str> names, Str_I path)
{
    for (const auto & entry : std::filesystem::directory_iterator(path)) {
        std::stringstream ss;
        if (entry.is_directory())
            continue;
        auto temp = entry.path().filename();
        ss << temp;
        string str = ss.str();
        str = str.substr(1, str.size() - 2);
        // str = str.substr(path.size(), str.size()); // remove " " and path
        names.push_back(str);
    }
}
#else
inline void file_list(vector_O<Str> fnames, Str_I path)
{
    SLS_ERR("not implemented");
}
#endif
#endif

// choose files with a given extension from a list of files
inline void file_ext(vector_O<Str> fnames_ext, vector_I<Str> fnames, Str_I ext, Bool_I keep_ext = true)
{
    fnames_ext.resize(0);
    Long N_ext = ext.size();
    for (Long i = 0; i < Size(fnames); ++i) {
        const Str & str = fnames[i];
        // check position of '.'
        Long ind = fnames[i].size() - N_ext - 1;
        if (ind < 0 || str[ind] != '.')
            continue;
        // match extension
        if (str.rfind(ext) != str.size() - ext.size())
            continue;
        if (keep_ext)
            fnames_ext.push_back(str);
        else
            fnames_ext.push_back(str.substr(0, str.size() - N_ext - 1));
    }
}

// list all files in current directory, with a given extension
inline void file_list_ext(vector_O<Str> fnames, Str_I path, Str_I ext, Bool_I keep_ext = true)
{
    vector<Str> fnames0;
    file_list(fnames0, path);
    file_ext(fnames, fnames0, ext, keep_ext);
}

// copy a file (read then write)
inline void file_copy(Str_I fname_out, Str_I fname_in, Bool_I replace = false)
{
    if (!file_exist(fname_in))
        SLS_ERR("file not found!");
    if (file_exist(fname_out) && !replace) {
        while (true) {
            if (file_exist(fname_out)) {
                SLS_WARN("\n\nfile [" + fname_out + "] already exist! delete file to continue...\n"
                    "  (set argument `replace = false` to replace file by default)\n\n");
            }
            else {
                break;
            }
            pause(10);
        }
    }
    ifstream fin(fname_in, std::ios::binary);
    ofstream fout(fname_out, std::ios::binary);
    fout << fin.rdbuf();
    fin.close();
    fout.close();
}
} // namespace slisc
