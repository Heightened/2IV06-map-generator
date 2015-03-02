#include "BlobMapShaper.h"

#include <cstdlib>
#include <cmath>

BlobMapShaper::BlobMapShaper(): MapShaper() {}

bool BlobMapShaper::isLand(glm::vec2 p) {
	float eye1 = glm::length(glm::vec2(p.x-0.2, p.y/2+0.2)) < 0.05;
	float eye2 = glm::length(glm::vec2(p.x+0.2, p.y/2+0.2)) < 0.05;
	float body = glm::length(p) < 0.8 - 0.18 * sin(5 * atan2(p.y, p.x));
	return body && !eye1 && !eye2;
};
