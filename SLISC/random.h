// generate random numbers

#pragma once
#include "arithmatic.h"
#include <ctime>

#ifndef _RAND_SEED_
#define _RAND_SEED_ std::time(nullptr)
#endif

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
			Ran(Ullong j = _RAND_SEED_) :
				v(4101842887655102017LL), w(1) {
				u = j ^ v; int64();
				v = u; int64();
				w = v; int64();
			}
			Doub doub() { return 5.42101086242752217E-20 * int64(); }
		};
	}

	// generate random Doub in [0, 1]
	inline Doub rand()
	{
		static internal::Ran rand_gen;
		return rand_gen.doub();
	}

	// generate random Int in {0,1,2,...,N-1}
	inline Int randInt(Int N)
	{ return Int(round(N*rand() - 0.5)); }

	// generate a random permutation of {0,1,2,...,N-1}
	inline void randPerm(VecInt_O perm, Int_I N)
	{
		Int j, n, ind;
		VecInt pool;
		linspace(pool, 0, N - 1, N);
		perm.resize(N);
		for (n = N; n > 0; --n) {
			ind = randInt(n);
			perm[n - 1] = pool(ind);
			for (j = ind; j < n - 1; ++j)
				pool(j) = pool(j + 1);
		}
	}

	inline Comp randComp()
	{ return Comp(rand(), rand()); }

	template <typename T>
	inline void rand(Vbase<T> &v)
	{
		Long i, N = v.size();
		for (i = 0; i < N; ++i)
			v(i) = rand();
	}

	// complex random number
	// uniform distribution
	template <>
	inline void rand<Comp>(Vbase<Comp> &v)
	{
		Long i, N = v.size();
		for (i = 0; i < N; ++i)
			v(i) = randComp();
	}
}
