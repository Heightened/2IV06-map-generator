#pragma once

#include "PointSelector.h"

/**
 * Selects points that will result in square polygons
 */
class SquarePointSelector: public PointSelector {
	public:
		SquarePointSelector(int width, int height): PointSelector(width, height) {}
		virtual std::vector<glm::vec2> select(int number);
};
