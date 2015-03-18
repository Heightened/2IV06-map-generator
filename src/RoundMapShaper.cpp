#include "RoundMapShaper.h"

RoundMapShaper::RoundMapShaper(): MapShaper() {}

bool RoundMapShaper::isLand(glm::vec2 p) {
	return (p.x*p.x + p.y*p.y) < 1;
};
