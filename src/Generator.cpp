#include "Generator.h"

#include <cstdio>
#include <cfloat>

#include "HexPointSelector.h"
#include "RandomPointSelector.h"
#include "PoissonPointSelector.h"
#include "vendor/VoronoiDiagramGenerator.h"

#define GENERATOR_MIN_DISTANCE 1.0e-6

struct vec2comp {
	bool operator() (const glm::vec2& lhs, const glm::vec2& rhs) const{
		if (lhs.x == rhs.x) {
			return lhs.y > rhs.y;
		} else {
			return lhs.x > rhs.x;
		}
	}
};

Generator::Generator(int _width, int _height, int _sampleSize) {
	width = _width;
	height = _height;
	sampleSize = _sampleSize;

	pointType = POINTSELECTOR_RANDOM;
};

PointSelector *Generator::shape() {
	switch (pointType) {
	case POINTSELECTOR_HEX:
		return new HexPointSelector(width, height);
	case POINTSELECTOR_POISSON:
		return new PoissonPointSelector(width, height);
	case POINTSELECTOR_RANDOM:
	default:
		return new RandomPointSelector(width, height);
	}
};

std::vector<glm::vec2> Generator::placePoints(PointSelector *psel) {
	return psel->select(sampleSize);
};

Map::Corner Generator::makeCorner(std::map<int, std::vector<Map::Corner> > &cornerMap, std::vector<Map::Corner> &corners, glm::vec2 p) {
	//Check if a corner already exists
	int bucket;
	for (bucket = (int)p.x - 1; bucket <= (int)p.x + 1; bucket++) {
		if (cornerMap.find(bucket) == cornerMap.end()){
			cornerMap.insert(std::pair<int, std::vector<Map::Corner> > (bucket, std::vector<Map::Corner>()));
		}

		for(std::vector<Map::Corner>::iterator it = cornerMap.find(bucket)->second.begin(); it != cornerMap.find(bucket)->second.end(); ++it) {
			float dx = p.x - it->point.x;
			float dy = p.y - it->point.y;
			if (dx * dx + dy * dy < GENERATOR_MIN_DISTANCE) {
				return *it;
			}
		}
	}
	bucket = (int)p.x;

	//Create a new corner
	Map::Corner q(corners.size(), p);
	q.border = ((int)p.x == 0 || (int)p.x == width || (int)p.y == 0 || (int)p.y == height);
	corners.push_back(q);
	cornerMap.find(bucket)->second.push_back(q);

	return q;
};

void Generator::buildGraph(std::vector<glm::vec2> points) {
	//Voronoi
	int numPoints = points.size();

	float* xValues = new float[numPoints];
	float* yValues = new float[numPoints];

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
	std::map<glm::vec2, Map::Center, vec2comp> centers;
	std::vector<Map::Edge> edges;
	std::map<int, std::vector<Map::Corner> > cornerMap;
	std::vector<Map::Corner> corners;

	for(std::vector<glm::vec2>::iterator it = points.begin(); it != points.end(); ++it) {
		Map::Center c(centers.size(), *it);
		centers.insert(std::pair<glm::vec2, Map::Center>(*it, c));

		if (polygonGraph) {
			polygonGraph->AddNode(*it);
		}
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

		glm::vec2 midway((v0x + v1x) / 2, (v0y + v1y) / 2);

		Map::Edge e(edges.size(),
			   centers[glm::vec2(d0x, d0y)],
			   centers[glm::vec2(d1x, d1y)],
			   makeCorner(cornerMap, corners, glm::vec2(v0x, v0y)),
			   makeCorner(cornerMap, corners, glm::vec2(v1x, v1y)),
			   midway);

		//TODO: Connect to neighbours

		edges.push_back(e);

		if (polygonGraph) {
			glm::vec2 a(v0x, v0y), b(v1x, v1y);

			polygonGraph->AddEdge(a, b);

			polygonGraph->AddNode(a);
			polygonGraph->AddNode(b);
		}
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
