#include "../SLISC/random.h"
#include "../SLISC/arithmatic.h"

void test_rand()
{
	// same seed
	Ran rangen1(1234), rangen2(1234);
	for (Int i = 0; i < 5; ++i)
		if (rangen1.doub() != rangen2.doub()) error("failed!");
	Ran rangen3;
	if (rand_gen.doub() == rangen3.doub()) error("failed!");
	if (rand_gen.doub() == rangen3.doub()) error("failed!");

	Long i, N = 10;
	VecDoub vx(N), vx1(N);
	rand(vx); rand(vx1);
	if (vx == vx1) error("failed!");
	if (max(vx) > 1.) error("failed!");
	if (max(vx) < 0.) error("failed!");
	VecComp vc(N), vc1(N);
	rand(vc); rand(vc1);
	if (vc == vc1) error("failed!");
}
