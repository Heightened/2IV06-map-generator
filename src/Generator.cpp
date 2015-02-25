#include "Generator.h"

#include <cstdio>
#include <cfloat>

#include "HexPointSelector.h"
#include "vendor/VoronoiDiagramGenerator.h"

Generator::Generator(int _width, int _height, int _sampleSize) {
	width = _width;
	height = _height;
	sampleSize = _sampleSize;
};

PointSelector *Generator::shape() {
	return new HexPointSelector(width, height);
};

std::vector<glm::vec2> Generator::placePoints(PointSelector *psel) {
	return psel->select(sampleSize);
};

void Generator::buildGraph(std::vector<glm::vec2> points) {
	//Voronoi
	int numPoints = points.size();

	float xValues [numPoints];
	float yValues [numPoints];

	float minX = FLT_MAX;
	float maxX = FLT_MIN;
	float minY = FLT_MAX;
	float maxY = FLT_MIN;

	for (int i = 0; i < numPoints; i++) {
		xValues[i] = points[i].x;
		yValues[i] = points[i].y;

		if(points[i].x < minX) {
			minX = points[i].x;
		}

		if(points[i].x > maxX) {
			maxX = points[i].x;
		}

		if(points[i].y < minY) {
			minY = points[i].y;
		}

		if(points[i].y > maxY) {
			maxY = points[i].y;
		}
	}

	VoronoiDiagramGenerator *voronoiGen = new VoronoiDiagramGenerator();
	voronoiGen->setGenerateDelaunay(true);

	//TODO: Should {min,max}{X,Y} be the map size?
	voronoiGen->generateVoronoi(xValues, yValues, numPoints, minX, maxX, minY, maxY, 0, true);

	//buildGraph
	std::vector<Map::Center> centers;
	std::vector<Map::Center> edges;

	for(std::vector<glm::vec2>::iterator it = points.begin(); it != points.end(); ++it) {
		Map::Center c(centers.size(), *it);
		centers.push_back(c);
	}

	//Voronoi Edge
	float v0x;
	float v0y;
	float v1x;
	float v1y;
	//Delaunay Line
	float d0x;
	float d0y;
	float d1x;
	float d1y;

	voronoiGen->resetIterator();
	voronoiGen->resetDelaunayEdgesIterator();

	while(voronoiGen->getNext(v0x, v0y, v1x, v1y) && voronoiGen->getNextDelaunay(d0x, d0y, d1x, d1y)) {
		//v: (v0x, v0y) -> (v1x, v1y)
		//d: (d0x, d0y) -> (d1x, d1y)

	}

	//improveCorners
};

void Generator::addFeatures() {

};

void Generator::start() {
	//Shaping
	PointSelector *psel = shape();

	//Placing points
	std::vector<glm::vec2> points = placePoints(psel);

	//Building graph
	buildGraph(points);
};
