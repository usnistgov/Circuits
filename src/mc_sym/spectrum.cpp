#include "spectrum.h"

namespace bfl {

	int compare(const spectrum_distribution &d1, const spectrum_distribution &d2)
	{
		if		(d1.size() < d2.size())	return -1;
		else if (d1.size() > d2.size())	return 1;

		for (auto it1 = d1.begin(), it2 = d2.begin(); it1 != d1.end(); ++it1, ++it2)
		{
			if		(it1->first < it2->first)	return -1;
			else if (it1->first > it2->first)	return 1;
			else
			{
				if		(it1->second < it2->second)	return -1;
				else if (it1->second > it2->second)	return 1;
			}
		}

		return 0;
	}

	bool operator==(const spectrum_distribution &d1, const spectrum_distribution &d2) {

		return (compare(d1, d2) == 0);
	}

	bool operator<(const spectrum_distribution &d1, const spectrum_distribution &d2) {

		return (compare(d1, d2) == -1);
	}

}
