#include "HexPointSelector.h"

std::vector<glm::vec2> HexPointSelector::select(int number) {
	std::vector<glm::vec2> ret;

	for (int i = 0; i < width; i++) {
		for	(int j = 0; j < height; j++) {
			ret.push_back(glm::vec2(i, j));
		}
	}

	return ret;
}
