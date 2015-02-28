#include "PoissonPointSelector.h"

#include <cstdio>

#include "vendor/PDSampling.h"

std::vector<glm::vec2> PoissonPointSelector::select(int number) {
	std::vector<glm::vec2> ret;

	double radius = 0.1f;

	PDSampler *sampler;

	while(ret.size() < number){
		sampler = new PureSampler(radius);
		sampler->complete();
		int N = (int) sampler->points.size();

		for (int i = 0; i < N; i++) {
			//Fit generated noise to the target domain
			glm::vec2 p(((sampler->points[i].x/2.0f) + 0.5f) * width, ((sampler->points[i].y/2.0f) + 0.5f) * width);
			ret.push_back(p);
		}
	}


	return ret;
}
