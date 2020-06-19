#include "matrix.h"


uint64_t make_row_reduced(dmatrix &m)
{
	const uint64_t N = m.rows();
	const uint64_t M = m.cols();

	uint64_t rank{ 0 };

	for (uint64_t c = 0; c < M; c++)
	{
		for (uint64_t r = rank; r < N; r++)
		{
			if (get_bit(m[r], c))
			{
				//ct::swap_if<uint64_t>(m[r], m[rank], r != rank);
				//if (r != rank)
					std::swap(m[r], m[rank]);

				//for (uint64_t j = 0; j < rank; j++)
				//	ct::xor_if(m[j], m[rank], get_bit(m[j], col));

				for (uint64_t j = r + 1; j < N; j++)
					//ct::xor_if(m[j], m[rank], get_bit(m[j], c));
					if (get_bit(m[j], c))
						m[j] ^= m[rank];

				rank++;

				break;
			}
		}
	}

	return rank;

}

uint64_t rank(const dmatrix &mat)
{
	auto rr(mat);

	auto rank = make_row_reduced(rr);

	return rank;

}
