#include "topology.h"
#include "logger.h"

TopologyVector load_topologies(const std::string &fileName)
{
	TopologyVector vec;

	std::ifstream ifile(fileName);

	if (ifile.is_open())
	{
		while (!ifile.eof())
		{
			std::string line;

			std::getline(ifile, line);

			//LOG_INFO << "read line : " << line << LOG_ENDL;

			if (!line.empty())
				vec.push_back(Topology(line));
		}
	}

	LOG_INFO << "loaded " << vec.size() << " topologies" << LOG_ENDL;

	return vec;
}

TopologyVector load_topologies(uint64_t K) 
{
	if (K >= 1 && K <= 6)
		return load_topologies("../data/topologies/topologies_k" + std::to_string(K) + ".txt");
	else
		return {};
}

