#include "Generator.h"

#include <cstdio>
#include <cfloat>
#include <queue>
#include <algorithm>

#include "HexPointSelector.h"
#include "SquarePointSelector.h"
#include "RandomPointSelector.h"
#include "PoissonPointSelector.h"
#include "SquareMapShaper.h"
#include "RoundMapShaper.h"
#include "RadialMapShaper.h"
#include "BlobMapShaper.h"
#include "CrescentMapShaper.h"
#include "vendor/VoronoiDiagramGenerator.h"

#define GENERATOR_MIN_DISTANCE 1.0e-6
#define LAKE_THRESHOLD 0.3f

struct vec2comp {
	bool operator() (const glm::vec2& lhs, const glm::vec2& rhs) const{
		if (lhs.x == rhs.x) {
			return lhs.y > rhs.y;
		} else {
			return lhs.x > rhs.x;
		}
	}
};

struct CornerElevationComp {
	bool operator() (Map::Corner*& lhs, Map::Corner*& rhs) {
		return lhs->elevation < rhs->elevation;
	}
};

struct CornerMoistureComp {
	bool operator() (Map::Corner*& lhs, Map::Corner*& rhs) {
		return lhs->moisture < rhs->moisture;
	}
};


Generator::Generator(int _width, int _height, int _sampleSize): centers(), edges(), corners() {
	width = _width;
	height = _height;
	sampleSize = _sampleSize;

	spring_count = width/6;

	pointType = POINTSELECTOR_RANDOM;
	shaperType = MAPSHAPER_RADIAL;
};

PointSelector *Generator::select() {
	switch (pointType) {
	case POINTSELECTOR_HEX:
		return new HexPointSelector(width, height);
	case POINTSELECTOR_SQUARE:
		return new SquarePointSelector(width, height);
	case POINTSELECTOR_POISSON:
		return new PoissonPointSelector(width, height);
	case POINTSELECTOR_RANDOM:
	default:
		return new RandomPointSelector(width, height);
	}
};

MapShaper *Generator::shape() {
	switch (shaperType) {
	case MAPSHAPER_ROUND:
		return new RoundMapShaper();
	case MAPSHAPER_SQUARE:
		return new SquareMapShaper();
	case MAPSHAPER_BLOB:
		return new BlobMapShaper();
	case MAPSHAPER_CRESCENT:
		return new CrescentMapShaper();
	case MAPSHAPER_RADIAL:
	default:
		return new RadialMapShaper();
	}
};

std::vector<glm::vec2> Generator::placePoints(PointSelector *psel) {
	return psel->select(sampleSize);
};

Map::Corner *Generator::makeCorner(std::map<int, std::vector<Map::Corner*> > &cornerMap, std::vector<Map::Corner*> &corners, glm::vec2 p) {
	//Check if a corner already exists
	int bucket;
	for (bucket = (int)p.x - 1; bucket <= (int)p.x + 1; bucket++) {
		if (cornerMap.find(bucket) == cornerMap.end()){
			cornerMap.insert(std::pair<int, std::vector<Map::Corner*> > (bucket, std::vector<Map::Corner*>()));
		}

		for(std::vector<Map::Corner*>::iterator it = cornerMap.find(bucket)->second.begin(); it != cornerMap.find(bucket)->second.end(); ++it) {
			float dx = p.x - (*it)->point.x;
			float dy = p.y - (*it)->point.y;
			if (dx * dx + dy * dy < GENERATOR_MIN_DISTANCE) {
				return *it;
			}
		}
	}
	bucket = (int)p.x;

	//Create a new corner
	Map::Corner *q = new Map::Corner(corners.size(), p);
	//You're a border if you're within 2% of the map sizes
	q->border = (p.x < width * 0.02 || p.x > width - (width * 0.02) || p.y < height * 0.02 || p.y > height - (height * 0.02));
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
	std::map<glm::vec2, Map::Center*, vec2comp> centerMap;
	std::map<int, std::vector<Map::Corner*> > cornerMap;

	for(std::vector<glm::vec2>::iterator it = points.begin(); it != points.end(); ++it) {
		Map::Center *c = new Map::Center(centerMap.size(), *it);
		centerMap.insert(std::pair<glm::vec2, Map::Center*>(*it, c));
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

		//Edges point to corners, edges point to centers
		Map::Edge* edge = new Map::Edge(edges.size(),
			   centerMap[glm::vec2(d0x, d0y)],
			   centerMap[glm::vec2(d1x, d1y)],
			   makeCorner(cornerMap, corners, glm::vec2(v0x, v0y)),
			   makeCorner(cornerMap, corners, glm::vec2(v1x, v1y)),
			   midway);

		//Centers point to edges, corners point to edges
		edge->d0->borders.push_back(edge);
		edge->d1->borders.push_back(edge);
		edge->v0->protrudes.push_back(edge);
		edge->v1->protrudes.push_back(edge);

		//Centers point to centers
		edge->d0->neighbours.insert(edge->d1);
		edge->d1->neighbours.insert(edge->d0);

		//Corners point to corners
		edge->v0->adjacent.insert(edge->v1);
		edge->v1->adjacent.insert(edge->v0);

		//Centers point to corners
		edge->d0->corners.insert(edge->v0);
		edge->d0->corners.insert(edge->v1);
		edge->d1->corners.insert(edge->v0);
		edge->d1->corners.insert(edge->v1);

		//Corners point to centers
		edge->v0->touches.insert(edge->d0);
		edge->v0->touches.insert(edge->d1);
		edge->v1->touches.insert(edge->d0);
		edge->v1->touches.insert(edge->d1);

		edges.push_back(edge);
	}

	// Dump keys from centerMap
	for (std::map<glm::vec2, Map::Center*, vec2comp>::iterator it = centerMap.begin(); it != centerMap.end(); it++) {
		centers.push_back(it->second);
	}

	if (polygonGraph) {
		// Add centers to graph
		for (std::map<glm::vec2, Map::Center*, vec2comp>::iterator it = centerMap.begin(); it != centerMap.end(); it++) {
			polygonGraph->AddNode(it->second->point);
		}

		// Add edges to graph
		for (std::vector<Map::Edge*>::iterator it = edges.begin(); it != edges.end(); it++) {
			polygonGraph->AddEdge((*it)->v0->point, (*it)->v1->point);
		}
	}

	//improveCorners
	//TODO
};

void Generator::addFeatures(MapShaper* mshape) {
	assignElevations(mshape);
	assignMoisture();
	assignBiomes();
};

void Generator::assignElevations(MapShaper* mshape) {
	assignElevationsCorner(mshape);
	assignElevationsCoastAndLand();
	assignElevationsRedistribute();
	assignElevationsPolygons();

	if (heightGraph) {
		// Add centers to graph
		for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
			heightGraph->AddNode(glm::vec3((*it)->point.x, (*it)->point.y, (*it)->elevation));
		}

		// Add edges to graph
		for (std::vector<Map::Edge*>::iterator it = edges.begin(); it != edges.end(); it++) {
			heightGraph->AddEdge(glm::vec3((*it)->v0->point.x, (*it)->v0->point.y, (*it)->v0->elevation), glm::vec3((*it)->v1->point.x, (*it)->v1->point.y, (*it)->v1->elevation));
		}
	}
};

void Generator::assignElevationsCorner(MapShaper* mshape) {
	std::queue<Map::Corner*> queue;

	for (std::vector<Map::Corner*>::iterator it = corners.begin(); it != corners.end(); it++) {
		//Island shape

		//Transform point to (-1, +1) range
		glm::vec2 p((*it)->point);
		p.x = (p.x/width)*2.0f - 1.0f;
		p.y = (p.y/height)*2.0f - 1.0f;

		(*it)->water = !mshape->isLand(p);

		if ((*it)->border) {
			(*it)->elevation = 0.0f;
			queue.push(*it);
		} else {
			(*it)->elevation = FLT_MAX;
		}
	}

	Map::Corner *q;
	while (!queue.empty()) {
		q = queue.front();
		queue.pop();

		for (std::set<Map::Corner*>::iterator it = q->adjacent.begin(); it != q->adjacent.end(); it++) {
			float newElevation = q->elevation + 0.01f;

			if (!(q->water || (*it)->water)) {
				//TODO: Don't always increase elevation inlands
				newElevation += 1;
				//TODO: Extra randomness needed for hex
			}

			if (newElevation < (*it)->elevation) {
				(*it)->elevation = newElevation;
				queue.push(*it);
			}
		}
	}
}

void Generator::assignElevationsCoastAndLand() {
	std::queue<Map::Center*> queue;

	// Compute polygon attributes 'ocean' and 'water' based on the
	// corner attributes. Count the water corners per
	// polygon. Oceans are all polygons connected to the edge of the
	// map. In the first pass, mark the edges of the map as ocean;
	// in the second pass, mark any water-containing polygon
	// connected an ocean as ocean.
	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		int numWater = 0;

		for (std::set<Map::Corner*>::iterator cit = (*it)->corners.begin(); cit != (*it)->corners.end(); cit++) {
			if ((*cit)->border) {
				(*it)->border = true;
				(*it)->ocean = true;
				(*it)->water = true;
				queue.push(*it);
			}
			if ((*cit)->water) {
				numWater += 1;
			}
		}

		(*it)->water = (*it)->ocean || numWater >= (*it)->corners.size() * LAKE_THRESHOLD;
	}

	Map::Center *c;
	while (!queue.empty()) {
		c = queue.front();
		queue.pop();

		for (std::set<Map::Center*>::iterator it = c->neighbours.begin(); it != c->neighbours.end(); it++) {
			if ((*it)->water && !(*it)->ocean) {
				(*it)->ocean = true;
				queue.push(*it);
			}
		}
	}

	// Set the polygon attribute 'coast' based on its neighbors. If
	// it has at least one ocean and at least one land neighbor,
	// then this is a coastal polygon.
	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		int numOcean = 0;
		int numLand = 0;

		for (std::set<Map::Center*>::iterator nit = (*it)->neighbours.begin(); nit != (*it)->neighbours.end(); nit++) {
			numOcean += (*nit)->ocean;
			numLand += !(*nit)->water;
		}

		(*it)->coast = (numOcean > 0) && (numLand > 0);
	}

	// Set the corner attributes based on the computed polygon
	// attributes. If all polygons connected to this corner are
	// ocean, then it's ocean; if all are land, then it's land;
	// otherwise it's coast.
	for (std::vector<Map::Corner*>::iterator it = corners.begin(); it != corners.end(); it++) {
		int numOcean = 0;
		int numLand = 0;
		for (std::set<Map::Center*>::iterator cit = (*it)->touches.begin(); cit != (*it)->touches.end(); cit++) {
			numOcean += (*cit)->ocean;
			numLand += !(*cit)->water;
		}
		(*it)->ocean = (numOcean == (*it)->touches.size());
		(*it)->coast = (numOcean > 0) && (numLand > 0);
		(*it)->water = (*it)->border || ((numLand != (*it)->touches.size()) && !(*it)->coast);

		// Set the elevation of ocean and coast tiles to 0
		if ((*it)->ocean || (*it)->coast) {
			(*it)->elevation = 0.0f;
		}
	}
}

void Generator::assignElevationsRedistribute() {
	std::vector<Map::Corner*> landCorners;

	for (std::vector<Map::Corner*>::iterator it = corners.begin(); it != corners.end(); it++) {
		if (!((*it)->ocean || (*it)->coast)) {
			landCorners.push_back(*it);
		}
	}
	std::sort(landCorners.begin(), landCorners.end(), CornerElevationComp());

	float SCALE_FACTOR = 1.1f;
	int i = 0;
	float x;
	float y;
	for (std::vector<Map::Corner*>::iterator it = landCorners.begin(); it != landCorners.end(); it++) {
		y = (float)i / (landCorners.size() - 1);

		x = sqrt(SCALE_FACTOR) - sqrt(SCALE_FACTOR * (1-y));

		if (x > 1.0f) {
			x = 1.0f;
		}

		(*it)->elevation = x;

		i++;
	}
}

void Generator::assignElevationsPolygons() {
	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		float sum = 0.0f;
		for (std::set<Map::Corner*>::iterator cit = (*it)->corners.begin(); cit != (*it)->corners.end(); cit++) {
			sum += (*cit)->elevation;
		}
		(*it)->elevation = sum/(*it)->corners.size();
	}
}

void Generator::calculateDownslopes() {
	Map::Corner *r;
	for (std::vector<Map::Corner*>::iterator it = corners.begin(); it != corners.end(); it++) {
		r = *it;
		for (std::set<Map::Corner*>::iterator ait = (*it)->adjacent.begin(); ait != (*it)->adjacent.end(); ait++) {
			if ((*ait)->elevation <= r->elevation) {
				r = *ait;
			}
		}
		(*it)->downslope = r;
	}
}

void Generator::calculateWatersheds() {
	for (std::vector<Map::Corner*>::iterator it = corners.begin(); it != corners.end(); it++) {
		(*it)->watershed = *it;
		if (!(*it)->ocean && !(*it)->coast) {
			(*it)->watershed = (*it)->downslope;
		}
	}

	for (int i = 0; i < 100; i++) {
		bool changed = false;
		Map::Corner *r;

		for (std::vector<Map::Corner*>::iterator it = corners.begin(); it != corners.end(); it++) {
			if (!((*it)->ocean || (*it)->coast || (*it)->watershed->coast)) {
				r = (*it)->downslope->watershed;
				if (!r->ocean) {
					(*it)->watershed = r;
					changed = true;
				}
			}
		}

		if (!changed) {
			break;
		}
	}

	for (std::vector<Map::Corner*>::iterator it = corners.begin(); it != corners.end(); it++) {
		(*it)->watershed->watershed_size += 1;
	}
}

static Map::Edge* lookupEdgeFromCorner(Map::Corner* q, Map::Corner* s) {
	for (std::vector<Map::Edge*>::iterator it = q->protrudes.begin(); it != q->protrudes.end(); it++) {
		if ((*it)->v0 == s || (*it)->v1 == s) {
			return *it;
		}
	}
	return new Map::Edge(-1, NULL, NULL, NULL, NULL, glm::vec2(-1, -1));
}

void Generator::createRivers() {
	for (int i = 0; i < spring_count; i++) {
		std::vector<Map::Corner*>::iterator it = corners.begin();
		std::advance(it, std::rand() % corners.size());
		Map::Corner *q = *it;

		if (q->ocean || q->elevation < 0.3 || q->elevation > 0.9) {
			continue;
		}

		while (!q->coast) {
			if (q == q->downslope) {
				break;
			}
			Map::Edge* e = lookupEdgeFromCorner(q, q->downslope);
			e->river += 1;
			q->river += 1;
			q->downslope->river += 1;
			q = q->downslope;
		}
	}
}

void Generator::assignCornerMoisture() {
	std::queue<Map::Corner*> queue;
	float newMoisture;

	for (std::vector<Map::Corner*>::iterator it = corners.begin(); it != corners.end(); it++) {
		if (((*it)->water || (*it)->river > 0) && !(*it)->ocean) {
			(*it)->moisture = ((*it)->river > 0) ? std::min(3.0f, (0.2f * (*it)->river)) : 1.0f;
			queue.push(*it);
		} else {
			(*it)->moisture = 0.0f;
		}
	}

	Map::Corner *c;
	while (!queue.empty()) {
		c = queue.front();
		queue.pop();

		for (std::set<Map::Corner*>::iterator it = c->adjacent.begin(); it != c->adjacent.end(); it++) {
			newMoisture = c->moisture * 0.9;
			if (newMoisture > (*it)->moisture) {
				(*it)->moisture = newMoisture;
				queue.push(*it);
			}
		}
	}

	for (std::vector<Map::Corner*>::iterator it = corners.begin(); it != corners.end(); it++) {
		if ((*it)->ocean || (*it)->coast) {
			(*it)->moisture = 1.0f;
		}
	}
}

void Generator::assignMoistureRedistribute() {
	std::vector<Map::Corner*> landCorners;

	for (std::vector<Map::Corner*>::iterator it = corners.begin(); it != corners.end(); it++) {
		if (!((*it)->ocean || (*it)->coast)) {
			landCorners.push_back(*it);
		}
	}
	std::sort(landCorners.begin(), landCorners.end(), CornerMoistureComp());

	int i = 0;
	for (std::vector<Map::Corner*>::iterator it = landCorners.begin(); it != landCorners.end(); it++) {
		(*it)->moisture = (float) i / (landCorners.size() - 1);
		i++;
	}
}

void Generator::assignPolygonMoisture() {
	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		float sum = 0.0f;
		for (std::set<Map::Corner*>::iterator cit = (*it)->corners.begin(); cit != (*it)->corners.end(); cit++) {
			sum += (*cit)->moisture;
		}
		(*it)->moisture = sum/(*it)->corners.size();
	}
}


void Generator::assignMoisture() {
	calculateDownslopes();
	calculateWatersheds();
	createRivers();
	assignCornerMoisture();
	assignMoistureRedistribute();
	assignPolygonMoisture();
};

//TODO: Abstract BiomePicker with different implementations?
static Map::Biome getBiome(Map::Center p) {
	if (p.ocean) {
		return Map::OCEAN;
	} else if (p.water) {
		if (p.elevation < 0.1) {
			return Map::MARSH;
		}
		if (p.elevation > 0.8) {
			return Map::ICE;
		}
		return Map::LAKE;
	} else if (p.coast) {
		return Map::BEACH;
	} else if (p.elevation > 0.8) {
		if (p.moisture > 0.5) {
			return Map::SNOW;
		} else if (p.moisture > 0.33) {
			return Map::TUNDRA;
		} else if (p.moisture > 0.16) {
			return Map::BARE;
		} else {
			return Map::SCORCHED;
		}
	} else if (p.elevation > 0.6) {
		if (p.moisture > 0.66) {
			return Map::TAIGA;
		} else if (p.moisture > 0.33) {
			return Map::SHRUBLAND;
		} else {
			return Map::TEMPERATE_DESERT;
		}
	} else if (p.elevation > 0.3) {
		if (p.moisture > 0.83) {
			return Map::TEMPERATE_RAIN_FOREST;
		} else if (p.moisture > 0.5) {
			return Map::TEMPERATE_DECIDUOUS_FOREST;
		} else if (p.moisture > 0.16) {
			return Map::GRASSLAND;
		} else {
			return Map::TEMPERATE_DESERT;
		}
	} else {
		if (p.moisture > 0.66) {
			return Map::TROPICAL_RAIN_FOREST;
		} else if (p.moisture > 0.33) {
			return Map::TROPICAL_SEASONAL_FOREST;
		} else if (p.moisture > 0.16) {
			return Map::GRASSLAND;
		} else {
			return Map::SUBTROPICAL_DESERT;
		}
	}
}

void Generator::assignBiomes() {
	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		(*it)->biome = getBiome(**it);
	}
};

void Generator::start() {
	//Shaping
	MapShaper *mshape = shape();
	PointSelector *psel = select();

	//Placing points
	std::vector<glm::vec2> points = placePoints(psel);

	//Building graph
	buildGraph(points);

	//Add features
	addFeatures(mshape);
};
