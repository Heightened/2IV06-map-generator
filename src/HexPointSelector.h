#pragma once

#include "PointSelector.h"

/**
 * Selects points that will result in hexagonal polygons
 */
class HexPointSelector: public PointSelector {
	public:
		HexPointSelector(int width, int height): PointSelector(width, height) {}
		virtual std::vector<glm::vec2> select(int number);
};
