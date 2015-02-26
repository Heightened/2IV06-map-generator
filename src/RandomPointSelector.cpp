#include "RandomPointSelector.h"

#include <cmath>
#include <ctime>

std::vector<glm::vec2> RandomPointSelector::select(int number) {
	std::vector<glm::vec2> ret;

	srand(time(NULL));

	for (int i = 0; i < number; i++) {
		ret.push_back(glm::vec2(
			10 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(width - 10 - 10))),
			10 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(height - 10 - 10)))
		));
	}

	return ret;
}
