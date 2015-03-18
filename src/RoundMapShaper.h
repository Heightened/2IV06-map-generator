#pragma once

#include "MapShaper.h"

/**
 * Shapes the map into a round form
 */
class RoundMapShaper: public MapShaper {
	public:
		RoundMapShaper();
		virtual bool isLand(glm::vec2);
};
