#include "GraphVisualisation.h"

#include <cmath>

#include <wx/log.h> 

Graph::Graph() {
	edgeCount = 0;
	nodeCount = 0;
}

void Graph::AddEdge(float x0, float x1, float y0, float y1, float z0, float z1) {
	float* edge = new float[6];
	edge[0] = x0; edge[1] = x1; edge[2] = y0; edge[3] = y1; edge[4] = z0; edge[5] = z1;
	edges.push_back(edge);
	edgeCount++;
	wxLogError(wxT("edge %i added!"), edgeCount);
}

void Graph::AddEdge(glm::vec3 a, glm::vec3 b) {
	AddEdge(a.x, b.x, a.y, b.y, a.z, b.z);
}

void Graph::AddEdge(glm::vec2 a, glm::vec2 b) {
	AddEdge(a.x, b.x, a.y, b.y, 0.0f, 0.0f);
}

void Graph::AddNode(float x, float y, float z) {
	float* node = new float[3];
	node[0] = x; node[1] = y; node[2] = z;
	nodes.push_back(node);
	nodeCount++;
	wxLogError(wxT("node %i added!"), nodeCount);
}

void Graph::AddNode(glm::vec3 node) {
	AddNode(node.x, node.y, node.z);
}

void Graph::AddNode(glm::vec2 node) {
	AddNode(node.x, node.y, 0.0f);
}

int Graph::getEdgeCount() {
	wxLogError(wxT("edgeCount: %i"), edgeCount);
	return edgeCount;
}

float* Graph::getEdge(int i) {
	return edges[i];
}

int Graph::getNodeCount() {
	wxLogError(wxT("nodeCount: %i"), nodeCount);
	return nodeCount;
}

float* Graph::getNode(int i) {
	return nodes[i];
}

GraphVertices::GraphVertices(Graph* g, int edgeCount, int nodeCount, float scale, float thickness) : Attribute(edgeCount * 6 * 6 * 3 + nodeCount * 20 * 3 * 3 * 4) {
	Attribute graphVertices = Attribute(0, 0);
	
	for (int i = 0; i < edgeCount; i++) {
		float* edge = g->getEdge(i);

		glm::vec3 from(edge[0],edge[2],edge[4]);
		glm::vec3 to(edge[1],edge[3],edge[5]);

		Attribute edgeVertices = BoxVertices(glm::distance(from, to), thickness, thickness);

		glm::vec3 initial(1,0,0);
		glm::vec3 desired(to - from);

		edgeVertices.translate(glm::vec3(0,0,-thickness*0.5f));

		glm::vec3 cross(glm::cross(initial, desired));
		if (glm::length(cross) != 0) {
			glm::vec3 axis(glm::normalize(cross));
			//glm::vec3 axis(0,1,0);
			float angle = acos(glm::dot(initial, desired) / (glm::length(initial) * glm::length(desired)));
			//float angle = 0;
			edgeVertices.rotateRad(axis, angle);
		}

		glm::vec3 position(from + desired * 0.5f);
		edgeVertices.translate(position);
		graphVertices = graphVertices + edgeVertices;
	}
	
	for (int i = 0; i < nodeCount; i++) {
		float* node = g->getNode(i);

		glm::vec3 position(node[0], node[1], node[2]);

		Attribute nodeVertices = SphereVertices(thickness, 1);
		nodeVertices.translate(position);
		graphVertices = graphVertices + nodeVertices;
	}
	
	graphVertices.scale(scale);

	values = new GLfloat[length];
	for (int i = 0; i < length; i++) {
		values[i] = graphVertices[i];
	}
}
