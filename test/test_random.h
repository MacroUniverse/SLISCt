#include "../SLISC/random.h"
#include "../SLISC/arithmatic.h"

void test_random()
{
	using slisc::internal::Ran;
	using slisc::rand;
	// same seed
	Ran rangen1(1234), rangen2(1234);
	Long i, N = 10;
	for (i = 0; i < 5; ++i)
		if (rangen1.doub() != rangen2.doub()) error("failed!");
	VecDoub vx(N), vx1(N);
	rand(vx); rand(vx1);
	if (vx == vx1) error("failed!");
	if (max(vx) > 1.) error("failed!");
	if (max(vx) < 0.) error("failed!");
	VecComp vc(N), vc1(N);
	rand(vc); rand(vc1);
	if (vc == vc1) error("failed!");
}
