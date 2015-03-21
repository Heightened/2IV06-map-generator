#include "RoadGenerator.h"

#include <queue>
#include <algorithm>

void RoadGenerator::generate() {
	road.clear();
	roadConnections.clear();

	std::queue<Map::Center*> queue;
	float elevationThresholds[4] = {0, 0.05f, 0.37f, 0.64f};
	int newLevel;

	std::map<Map::Corner*, int> cornerContour;
	std::map<Map::Center*, int> centerContour;

	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		if((*it)->coast || (*it)->ocean) {
			centerContour.insert(std::pair<Map::Center*, int>(*it, 1));
			queue.push(*it);
		} else {
			centerContour.insert(std::pair<Map::Center*, int>(*it, 0));
		}
	}

	Map::Center *q;
	while (!queue.empty()) {
		q = queue.front();
		queue.pop();

		for (std::set<Map::Center*>::iterator it = q->neighbours.begin(); it != q->neighbours.end(); it++) {
			newLevel = centerContour[q];
			while ((*it)->elevation >  elevationThresholds[newLevel] && !(*it)->water) {
				newLevel += 1;
			}
			if (newLevel < (centerContour[(*it)] == 0 ? 999 : centerContour[(*it)])) {
				centerContour[(*it)] = newLevel;
				queue.push(*it);
			}
		}
	}

	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		for (std::set<Map::Corner*>::iterator cit = (*it)->corners.begin(); cit != (*it)->corners.end(); cit++) {
			cornerContour[(*cit)] = std::min(cornerContour[(*cit)] == 0 ? 999 : cornerContour[(*cit)],
											 centerContour[(*it)] == 0 ? 999 : centerContour[(*it)]);
		}
	}

	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		for (std::vector<Map::Edge*>::iterator eit = (*it)->borders.begin(); eit != (*it)->borders.end(); eit++) {
			if (cornerContour[(*eit)->v0] != cornerContour[(*eit)->v1]) {
				road[(*eit)] = std::min(cornerContour[(*eit)->v0], cornerContour[(*eit)->v1]);
				roadConnections[(*it)].push_back(*eit);
			}
		}
	}
}
