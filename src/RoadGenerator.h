#include <map>
#include <vector>

#include "Map.h"

class RoadGenerator {
	std::map<Map::Edge*, int> road;
	std::map<Map::Center*, std::vector<Map::Edge*> > roadConnections;

	std::vector<Map::Center*> centers;

	public:
		RoadGenerator(std::vector<Map::Center*> _centers): road(), roadConnections() {
			centers = _centers;
		}

		void generate();

		std::map<Map::Center*, std::vector<Map::Edge*> > getRoadConnections() {
			return roadConnections;
		}

		std::map<Map::Edge*, int> getRoads() {
			return road;
		}
};
