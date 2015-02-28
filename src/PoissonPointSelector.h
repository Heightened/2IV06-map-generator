#pragma once

#include "PointSelector.h"

/**
 * Selects points using Poisson dart throwing
 */
class PoissonPointSelector: public PointSelector {
	public:
		PoissonPointSelector(int width, int height): PointSelector(width, height) {}
		virtual std::vector<glm::vec2> select(int number);
};
