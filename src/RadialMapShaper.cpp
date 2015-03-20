#include "RadialMapShaper.h"

#define PI 3.1415926535897932384626433832795
#define ISLAND_FACTOR 1.07f

#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <cstdio>

float floatBetween(float a, float b) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    float r = random * diff;
    return a + r;
}

RadialMapShaper::RadialMapShaper(): MapShaper() {
	srand(time(NULL));

	bumps = 1 + rand() % 5;
	startAngle = floatBetween(0, 2*PI);
	dipAngle = floatBetween(0, 2*PI);
	dipWidth = floatBetween(0.2f, 0.7f);
}

bool RadialMapShaper::isLand(glm::vec2 p) {
	float angle = atan2(p.y, p.x);
	float length = 0.5f * (std::max(std::abs(p.x), std::abs(p.y)) + glm::length(p));

	float r1 = 0.5f + (0.40f * sin(startAngle + bumps*angle + cos((bumps+3)*angle)));
	float r2 = 0.7f - (0.20f * sin(startAngle + bumps*angle - sin((bumps+2)*angle)));

	if (std::abs(angle - dipAngle) < dipWidth ||
		std::abs(angle - dipAngle + 2*PI) < dipWidth ||
		std::abs(angle - dipAngle - 2*PI) < dipWidth) {
		r1 = 0.2f;
		r2 = 0.2f;
	}

	return (length < r1 || (length > r1*ISLAND_FACTOR && length < r2));
};
