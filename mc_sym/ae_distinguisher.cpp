#include "ae_distinguisher.h"

namespace ae {

	template<>
	std::unique_ptr<class_distinguisher_base<6>> class_distinguisher() {

		return std::unique_ptr<class_distinguisher_base<6>>(new class_distinguisher_advanced<6>("../data/n6_group_data.txt"));
	}
}
