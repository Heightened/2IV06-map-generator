#pragma once

#ifndef PREDEFINEDOBJECTS
	#include "Objects.h"
#endif

#include <vector>

class Graph {
	int edgeCount;
	std::vector<float*> edges;

	int nodeCount;
	std::vector<float*> nodes;
public:
	Graph();
	//Tip: use glm vector based calls to avoid confusion
	void AddEdge(float x0, float x1, float y0, float y1, float z0 = 0, float z1 = 0);
	void AddEdge(glm::vec3 a, glm::vec3 b);
	void AddEdge(glm::vec2 a, glm::vec2 b);
	void AddNode(float x, float y, float z = 0);
	void AddNode(glm::vec3 node);
	void AddNode(glm::vec2 node);

	void RemoveEdge(int i);
	void RemoveNode(int i);

	//returns the number of edges that were added
	int getEdgeCount();
	//returns the 6 coordinates corresponding to the x0, x1, y0, y1, z0, z1 with which edge i was constructed
	float* getEdge(int i);

	//returns the number of nodes that were added
	int getNodeCount();
	//returns the 3 coordinates corresponding to the x, y and optional z with which node i was constructed
	float* getNode(int i);
};

//Divides graph into smaller sections to be visualized as multiple GraphVertices
std::vector<Graph*> DivideGraph(Graph* g);

class GraphVertices : public Attribute {
public:
	//constructs a graph using spheres for nodes and beams for edges
	GraphVertices(Graph* g, int edgeCount, int nodeCount, float scale, float thickness);
};
