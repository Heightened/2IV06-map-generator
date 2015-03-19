#include "CrescentMapShaper.h"

#include <ctime>
#include <cstdlib>
#include <algorithm>

CrescentMapShaper::CrescentMapShaper(): MapShaper(), offset() {
	srand(time(NULL));

	int first = rand() % 50 + 50;
	int second = rand() % 40 + 30;

	main_size = std::max(first, second)/100.0f;
	negative_size = std::min(first, second)/100.0f;

	offset.x = (rand() % 200)/100.0f - 1.0f;
	offset.y = (rand() % 200)/100.0f - 1.0f;
}

bool CrescentMapShaper::isLand(glm::vec2 p) {
	glm::vec2 q(p);

	q += offset;

	return (p.x*p.x + p.y*p.y) < main_size && //Main landmass circle
		!((q.x*q.x + q.y*q.y) < negative_size); //Negative circle
};
