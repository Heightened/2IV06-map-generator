#pragma once

#include "MapShaper.h"

/**
 * Shapes the map into a radial form
 */
class RadialMapShaper: public MapShaper {
	float ISLAND_FACTOR;

	int bumps;
	float startAngle;
	float dipAngle;
	float dipWidth;

	public:
		RadialMapShaper();
		virtual bool isLand(glm::vec2);
};
