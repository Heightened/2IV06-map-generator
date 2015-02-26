#pragma once

#include "PointSelector.h"

/**
 * Selects points that will result in random polygons
 */
class RandomPointSelector: public PointSelector {
	public:
		RandomPointSelector(int width, int height): PointSelector(width, height) {}
		virtual std::vector<glm::vec2> select(int number);
};
