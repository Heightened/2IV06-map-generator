#include "Generator.h"

#include <vector>
#include <glm/vec2.hpp>

#include "PointSelector.h"
#include "HexPointSelector.h"

#include <cstdio>

Generator::Generator() {

};

void Generator::start() {
	PointSelector *psel = new HexPointSelector(10, 10);
	std::vector<glm::vec2> points = psel->select(200);

	for(std::vector<glm::vec2>::iterator it = points.begin(); it != points.end(); ++it) {
		printf("(%f, %f)\n", it->x, it->y);
	}
};
