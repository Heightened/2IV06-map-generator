#pragma once

#include <vector>
#include <set>
#include <map>
#include <glm/glm.hpp>

#include "PointSelector.h"
#include "GraphVisualisation.h"

namespace Map {
	class Center;
	class Edge;
	class Corner;

	class Center {
		public:
			Center(): point(glm::vec2(0, 0)), neighbours(), borders(), corners() {
				index = -1;
			}

			Center(int i, glm::vec2 _point): point(_point), neighbours(), borders(), corners() {
				index = i;
			}

			Center(const Center &c): point(c.point), neighbours(c.neighbours), borders(c.borders), corners(c.corners) {
				index = c.index;

				water = c.water;
				ocean = c.ocean;
				coast = c.coast;
				border = c.border;
			}

			int index;

			glm::vec2 point;
			bool water;
			bool ocean;
			bool coast;
			bool border;

			std::vector<Center> neighbours;
			std::vector<Edge> borders;
			std::vector<Corner> corners;
	};

	class Corner {
		public:
			Corner(int i, glm::vec2 p): point(p), touches(), protrudes(), adjacent() {
				index = i;
			}

			int index;

			glm::vec2 point;
			bool ocean;
			bool water;
			bool coast;
			bool border;

			std::set<Center> touches;
			std::set<Edge> protrudes;
			std::set<Corner> adjacent;
	};

	class Edge {
		public:
			Edge(int i, Center _d0, Center _d1, Corner _v0, Corner _v1, glm::vec2 _midway): d0(_d0), d1(_d1), v0(_v0), v1(_v1), midway(_midway) {
				index = i;
			}

			int index;

			Center d0;
			Center d1;
			Corner v0;
			Corner v1;

			glm::vec2 midway;
	};
};

/**
 * Generates a new Island based on a set of parameters
 */
class Generator {
	// Width of the Island to be generated
	int width;
	// Height of the Island to be generated
	int height;
	// The number of samples to use to generate the terrain
	int sampleSize;

	PointSelector *shape();
	std::vector<glm::vec2> placePoints(PointSelector*);
	Graph * buildGraph(std::vector<glm::vec2>);
	void addFeatures();

	Map::Corner makeCorner(std::map<int, std::vector<Map::Corner> > &cornerMap, std::vector<Map::Corner> &corners, glm::vec2 p);

	public:
		/**
		 * @param in int mapWidth Width of the Island to be generated
		 * @param in int mapHeight Height of the Island to be generated
		 * @param in int sampleSize The number of samples to use to generate the terrain
		 */
		Generator(int mapWidth, int mapHeight, int sampleSize);
		// Starts the generation process
		Graph * start();
};
