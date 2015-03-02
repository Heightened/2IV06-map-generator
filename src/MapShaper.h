#pragma once

#include <vector>
#include <glm/glm.hpp>

/**
 * Abstract class used to shape a map
 */
class MapShaper {
	public:
		MapShaper(){}

		/**
		 * Checks if a point is part of the land
		 *
		 * @param in glm::vec2 The point to check
		 * @return bool Whether or not the given point is part of the land
		 */
		virtual bool isLand(glm::vec2) = 0;
};
