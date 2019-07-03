#pragma once
#include "global.h"

namespace slisc {

// see if elm == vec[i], return i
// return -1 if not found
template <class T, class Tv>
inline Long search(const T &elm, const Tv &vec)
{
	for (Long i = 0; i < vec.size(); ++i) {
		if (elm == vec[i])
			return i;
	}
	return -1;
}

// search ind so that v[ind] == s
// same as lookupInt(), but Int operator-(v[i], s) must be implemented
// operator-(v[i], s) only need to return the correct sign
// return 0, output ind: if s is found
// return -1: if s is not found and s < v[0]
// return 1: if s is not found and v[end] < s
// return -2: if s is not found and v[0] < s < v[end], output ind so that v[ind] < s < v[ind+1]
// return -3: if s is not found because v.size() == 0
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
