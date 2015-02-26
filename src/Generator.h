#pragma once

#include <vector>
#include <glm/glm.hpp>

#include "PointSelector.h"
#include "GraphVisualisation.h"

namespace Map {
	class Center;
	class Edge;
	class Corner;

	class Center {
		public:
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

	class Edge {
		public:
			Edge(int i, Center _d0, Center _d1, Center _v0, Center _v1, glm::vec2 _midway): d0(_d0), d1(_d1), v0(_v0), v1(_v1), midway(_midway) {
				index = i;
			}

			int index;

			Center d0;
			Center d1;
			Center v0;
			Center v1;

			glm::vec2 midway;
	};

	class Corner {
		public:
			int index;

			glm::vec2 point;
			bool ocean;
			bool water;
			bool coast;
			bool border;

			std::vector<Center> touches;
			std::vector<Edge> protrudes;
			std::vector<Corner> adjacent;
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
