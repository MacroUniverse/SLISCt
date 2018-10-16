// generate random numbers

#pragma once
#include "slisc.h"
#include <ctime>

namespace slisc
{
	namespace internal
	{
		// Algorithm from Numerical Recipes 3ed
		// empty construcotr uses "std:time()" as seed, avoid calling it twice in the same second.
		// it is best if the whole program uses only one object, which is currently in rand() function.
		class Ran
		{
		private:
			Ullong u, v, w;
			Uint int32() { return (Uint)int64(); }
			Ullong int64() {
				u = u * 2862933555777941757LL + 7046029254386353087LL;
				v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
				w = 4294957665U * (w & 0xffffffff) + (w >> 32);
				Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
				return (x + v) ^ w;
			}
		public:
			Ran(Ullong j = std::time(nullptr)) :
				v(4101842887655102017LL), w(1) {
				u = j ^ v; int64();
				v = u; int64();
				w = v; int64();
			}
			Doub doub() { return 5.42101086242752217E-20 * int64(); }
		};
	}

	inline Doub rand()
	{
		static internal::Ran rand_gen;
		return rand_gen.doub();
	}

	inline Comp randComp()
	{ return Comp(rand(), rand()); }

	template <typename T>
	inline void rand(Base<T> &v)
	{
		Long i, N = v.size();
		for (i = 0; i < N; ++i)
			v(i) = rand();
	}

	// complex random number
	// uniform distribution
	template <>
	inline void rand<Comp>(Base<Comp> &v)
	{
		Long i, N = v.size();
		for (i = 0; i < N; ++i)
			v(i) = randComp();
	}
}
