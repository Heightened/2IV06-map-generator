#include "SquareMapShaper.h"

SquareMapShaper::SquareMapShaper(): MapShaper() {}

bool SquareMapShaper::isLand(glm::vec2 p) {
	return true;
};
