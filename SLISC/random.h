#pragma once
#include "slisc.h"
#include <ctime>

namespace slisc
{
	namespace internal
	{
		// Algorithm from Numerical Recipes 3ed
		// empty construcotr uses "std:time()" as seed, avoid calling it twice in the same second.
		// it is best if the whole program uses only one object
		struct Ran
		{
			Ullong u, v, w;
			Ran(Ullong j = std::time(nullptr)) : v(4101842887655102017LL), w(1) {
				u = j ^ v; int64();
				v = u; int64();
				w = v; int64();
			}
			Ullong int64() {
				u = u * 2862933555777941757LL + 7046029254386353087LL;
				v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
				w = 4294957665U * (w & 0xffffffff) + (w >> 32);
				Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
				return (x + v) ^ w;
			}
			Doub doub() { return 5.42101086242752217E-20 * int64(); }
			Uint int32() { return (Uint)int64(); }
		};
	}

	inline Doub rand()
	{
		static internal::Ran rand_gen;
		return rand_gen.doub();
	}

	template <typename T>
	inline void rand(NRbase<T> &v)
	{
		Long i, N = v.size();
		for (i = 0; i < N; ++i)
			v(i) = rand();
	}

	template <>
	inline void rand<Comp>(NRbase<Comp> &v)
	{
		Long i, N = 2 * v.size();
		Doub *p = (Doub*)v.ptr();
		for (i = 0; i < N; ++i)
			p[i] = rand();
	}
}