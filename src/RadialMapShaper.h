#pragma once

#include "MapShaper.h"

/**
 * Shapes the map into a radial form
 */
class RadialMapShaper: public MapShaper {
	int bumps;
	float startAngle;
	float dipAngle;
	float dipWidth;

	public:
		RadialMapShaper();
		virtual bool isLand(glm::vec2);
};
