#include "ae_invariants.h"


namespace ae {

	using namespace bfl;

	bool operator==(const affine_invariants &ai1, const affine_invariants &ai2) {

		return	(ai1.degree == ai2.degree) &&
				(ai1.walsh_dist == ai2.walsh_dist) &&
				(ai1.ac_dist == ai2.ac_dist);
	}

	bool operator!=(const affine_invariants &ai1, const affine_invariants &ai2) {

		return !(ai1 == ai2);
	}

	bool operator<(const affine_invariants &ai1, const affine_invariants &ai2)
	{
		if (ai1.degree < ai2.degree)
			return true;
		else if (ai1.degree == ai2.degree)
		{
			auto c1 = compare(ai1.ac_dist, ai2.ac_dist);

			if ((c1 < 0) ||
				(c1 == 0 && (compare(ai1.walsh_dist, ai2.walsh_dist) < 0)))
				return true;
		}

		return false;
	}

}