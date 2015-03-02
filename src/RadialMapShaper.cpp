#include "RadialMapShaper.h"

#define PI 3.1415926535897932384626433832795
#define ISLAND_FACTOR 1.07f

#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>

RadialMapShaper::RadialMapShaper(): MapShaper() {
	srand(time(NULL));

	bumps = rand() % 6;
	startAngle = static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/2*PI));
	dipAngle = static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/2*PI));
	dipWidth = 0.2f - static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(0.7f - 0.2f)));
}

bool RadialMapShaper::isLand(glm::vec2 p) {
	float angle = atan2(p.y, p.x);
	float length = 0.5f * std::max(abs(p.x), abs(p.y)) + p.length();

	float r1 = 0.5f + 0.40f * sin(startAngle + bumps*angle + cos((bumps+3)*angle));
	float r2 = 0.7f - 0.20f * sin(startAngle + bumps*angle + cos((bumps+2)*angle));

	if (abs(angle - dipAngle) < dipWidth ||
		abs(angle - dipAngle + 2*PI) < dipWidth ||
		abs(angle - dipAngle - 2*PI) < dipWidth) {
		r1 = 0.2f;
		r2 = 0.2f;
	}

	return (length < r1 || (length > r1*ISLAND_FACTOR && length < r2));
};
