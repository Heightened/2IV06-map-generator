#pragma once

#include "MapShaper.h"

/**
 * Shapes the map into a round form
 */
class CrescentMapShaper: public MapShaper {
	glm::vec2 offset;
	float main_size;
	float negative_size;

	public:
		CrescentMapShaper();
		virtual bool isLand(glm::vec2);
};
