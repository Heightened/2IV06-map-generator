#pragma once

#include <vector>
#include <set>
#include <glm/glm.hpp>

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

				water = false;
				ocean = false;
				coast = false;
				border = false;

				elevation = 0.0f;
			}

			inline bool operator< (const Center& other) const{
				return index < other.index;
			}

			int index;

			glm::vec2 point;
			bool water;
			bool ocean;
			bool coast;
			bool border;

			float elevation;

			std::set<Center*> neighbours;
			std::vector<Edge> borders;
			std::set<Corner*> corners;
	};

	class Corner {
		public:
			Corner(int i, glm::vec2 p): point(p), touches(), protrudes(), adjacent() {
				index = i;

				water = false;
				ocean = false;
				coast = false;
				border = false;

				elevation = 0.0f;
			}

			inline bool operator< (const Corner& other) const{
				return index < other.index;
			}

			int index;

			glm::vec2 point;
			bool ocean;
			bool water;
			bool coast;
			bool border;

			float elevation;

			std::set<Center*> touches;
			std::vector<Edge> protrudes;
			std::set<Corner*> adjacent;
	};

	class Edge {
		public:
			Edge(int i, Center *_d0, Center *_d1, Corner *_v0, Corner *_v1, glm::vec2 _midway): midway(_midway) {
				index = i;

				d0 = _d0;
				d1 = _d1;
				v0 = _v0;
				v1 = _v1;
			}

			inline bool operator< (const Edge& other) const{
				return index < other.index;
			}

			int index;

			Center *d0;
			Center *d1;
			Corner *v0;
			Corner *v1;

			glm::vec2 midway;
	};
};
