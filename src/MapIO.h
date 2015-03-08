#pragma once

#include <vector>
#include <glm/glm.hpp>

#include "Map.h"

namespace IO {
	struct Center {
		int index;
		glm::vec2 point;
		bool water;
		bool ocean;
		bool coast;
		bool border;
		float elevation;
	};

	struct Corner {
		int index;
		glm::vec2 point;
		bool ocean;
		bool water;
		bool coast;
		bool border;
		float elevation;
	};

	struct Edge {
		int index;
		int center_d0;
		int center_d1;
		int corner_v0;
		int corner_v1;
		glm::vec2 midway;
	};

	// The neighbours relation on Center
	struct Neighbour_rel {
		// center_out has center_in in its neighbour list
		int center_out;
		int center_in;
	};

	struct Border_rel {
		int center_out;
		int edge_in;
	};

	struct Corner_rel {
		int center_out;
		int corner_in;
	};

	struct Touches_rel {
		int corner_out;
		int center_in;
	};

	struct Protrudes_rel {
		int corner_out;
		int edge_in;
	};

	struct Adjacent_rel {
		int corner_out;
		int corner_in;
	};

	void exportMap(FILE *file, std::vector<Map::Center*> centers);
	std::vector<Map::Center*> importMap(FILE *file);
};
