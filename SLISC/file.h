#pragma once
#include "global.h"
#include <sstream>
#include <fstream>
#include <codecvt>
#include <filesystem>
#include "unicode.h"

namespace slisc {

inline void file_list(vector_O<Str> names, Str_I path);

#ifdef _MSC_VER
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
#ifndef _MSC_VER
	ifstream f(fname.c_str());
	return f.good();
#else
	if (case_sens)
		return file_exist_case(fname);
	else {
		ifstream f(fname.c_str());
		return f.good();
	}
#endif
}

// not case sensitive on Windows, see file_exist_case()
inline Bool file_exist(Str32_I fname) {
	return file_exist(utf32to8(fname));
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
		if (i == 999) error("too many temporary files!");
		temp_fname = temp_fname_pref + to_string(i);
		if (!file_exist(temp_fname)) break;
	}
	
	// save a list of all files (no folder) to temporary file
	system(("ls -p " + path + " | grep -v / > " + temp_fname).c_str());

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
#endif

// std::filesystem implementation of file_list()
// works in Visual Studio, not gcc 8
// directory example: "C:/Users/addis/Documents/GitHub/SLISC/"
#ifdef _MSC_VER
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

// list all files in current directory, with a given extension
inline void file_list_ext(vector_O<Str32> fnames, Str32_I path, Str32_I ext, Bool_I keep_ext = true)
{
	vector<Str> fnames8;
	fnames.resize(0);
	file_list_ext(fnames8, utf32to8(path), utf32to8(ext), keep_ext);
	for (Long i = 0; i < Size(fnames8); ++i)
		fnames.push_back(utf8to32(fnames8[i]));
}

} // namespace slisc
