#include "slisc.h"

namespace slisc {

// search ind so that v[ind] == s
// "<" and "==" operator must exist between v[] and s
// if s is found, return 0, output ind
// if s is not found and s < v[0], return -1
// if s is not found and v[end] < s, return 1
// if s is not found and v[0] < s < v[end], return -2 and output ind so that v[ind] < s < v[ind+1]
// if v.size() < 1, return -3

template <class T, class T1>
Int lookupInt(Long_O &ind, const T1 &v, const T &s)
{
	Long i, N = v.size(), ind1 = 0, ind2 = N - 1;
	if (N < 1) return -3;
	if (s < v[0]) return -1;
	if (v[ind2] < s) return 1;
	if (s == v[0]) {
		ind = 0; return 0;
	}
	if (s == v[ind2]) {
		ind = ind2; return 0;
	}
	for (i = 0; i < N; ++i) {
		ind = (ind1 + ind2) / 2;
		if (s == v[ind]) return 0;
		if (s < v[ind]) ind2 = ind;
		ind1 = ind;
		if (ind2 - ind1 == 1) {
			ind = ind1; return -2;
		}
	}
	return 0;
}
} // namespace slisc
