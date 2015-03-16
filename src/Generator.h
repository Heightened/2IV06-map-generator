#pragma once

#include <vector>
#include <set>
#include <map>
#include <glm/glm.hpp>

#include "Map.h"
#include "PointSelector.h"
#include "MapShaper.h"
#include "GraphVisualisation.h"

/**
 * Point selector type
 */
enum PointSelectorType {
	POINTSELECTOR_RANDOM,
	POINTSELECTOR_HEX,
	POINTSELECTOR_POISSON,
	POINTSELECTOR_SQUARE
};

/**
 * Map shaper type
 */
enum MapShaperType {
	MAPSHAPER_RADIAL,
	MAPSHAPER_SQUARE,
	MAPSHAPER_BLOB
};

/**
 * Generates a new Island based on a set of parameters
 */
class Generator {
	// Width of the Island to be generated
	int width;
	// Height of the Island to be generated
	int height;
	// The number of samples to use to generate the terrain
	int sampleSize;

	// The point selector to use
	PointSelectorType pointType;
	// The map shaper to use
	MapShaperType shaperType;

	// A Graph representing the current polygons in the map
	Graph *polygonGraph;
	// A Graph representing the current polygons in the map, with height
	Graph *heightGraph;

	// Data structures representing the map
	std::vector<Map::Center*> centers;
	std::vector<Map::Edge*> edges;
	std::vector<Map::Corner*> corners;

	MapShaper *shape();
	PointSelector *select();
	std::vector<glm::vec2> placePoints(PointSelector*);
	void buildGraph(std::vector<glm::vec2>);
	void addFeatures(MapShaper*);
	void assignElevations(MapShaper*);
	void assignElevationsCorner(MapShaper*);
	void assignElevationsCoastAndLand();
	void assignElevationsRedistribute();
	void assignElevationsPolygons();
	void assignMoisture();
	void calculateDownslopes();
	void calculateWatersheds();
	void createRivers();
	void assignCornerMoisture();
	void assignMoistureRedistribute();
	void assignPolygonMoisture();
	void assignBiomes();

	Map::Corner *makeCorner(std::map<int, std::vector<Map::Corner*> > &cornerMap, std::vector<Map::Corner*> &corners, glm::vec2 p);

	public:
		/**
		 * @param in int mapWidth Width of the Island to be generated
		 * @param in int mapHeight Height of the Island to be generated
		 * @param in int sampleSize The number of samples to use to generate the terrain
		 */
		Generator(int mapWidth, int mapHeight, int sampleSize);

		// Resets the generator
		void reset() {
			centers.clear();
			edges.clear();
			corners.clear();
		};

		// Starts the generation process
		void start();

		/**
		 * Sets the point selector type to use
		 * @param in PointSelectorType type The new type
		 */
		void setPointSelectorType(PointSelectorType t) {
			pointType = t;
		}

		/**
		 * Sets the map shaper type to use
		 * @param in MapShaperType type The new type
		 */
		void setMapShaperType(MapShaperType t) {
			shaperType = t;
		}

		/**
		 * Returns a Graph representing the current polygons in the map
		 * @return Graph* The current polygons in the map
		 */
		Graph *getPolygonGraph() {
			return polygonGraph;
		}

		/**
		 * Sets the Graph representing the current polygons in the map
		 * @param in Graph* g The Graph instance to use
		 */
		void setPolygonGraph(Graph *g) {
			polygonGraph = g;
		}

		/**
		 * Returns a Graph representing the current polygons with height in the map
		 * @return Graph* The current polygons in the map with height
		 */
		Graph *getHeightGraph() {
			return heightGraph;
		}

		/**
		 * Sets the Graph representing the current polygons with height in the map
		 * @param in Graph* g The Graph instance to use
		 */
		void setHeightGraph(Graph *g) {
			heightGraph = g;
		}

		std::vector<Map::Center*> getCenters() {
			return centers;
		}

		float getWidth() {
			return (float) width;
		}

		float getHeight() {
			return (float) height;
		}
};
