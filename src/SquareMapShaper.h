#pragma once

#include "MapShaper.h"

/**
 * Shapes the map into a square form
 */
class SquareMapShaper: public MapShaper {
	public:
		SquareMapShaper();
		virtual bool isLand(glm::vec2);
};
