#include "HexPointSelector.h"

#include <cmath>

std::vector<glm::vec2> HexPointSelector::select(int number) {
	std::vector<glm::vec2> ret;

	int n = sqrt(number);

	for (int i = 0; i < n; i++) {
		for	(int j = 0; j < n; j++) {
			ret.push_back(glm::vec2(
				(0.5f + i)/n * width,
				(0.25f + 0.5f * (i % 2) + j)/n * height
			));
		}
	}

	return ret;
}
