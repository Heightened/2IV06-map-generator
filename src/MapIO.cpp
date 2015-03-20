#include <cstdio>
#include <ctime>
#include <algorithm>
#include <map>

#include "MapIO.h"

void IO::exportMap(FILE *file, std::vector<Map::Center*> centers, Map::Info info) {
	std::set<Map::Corner*> corners;
	std::set<Map::Edge*> edges;

	// Relationships
	std::vector<IO::Neighbour_rel> neighbours;
	std::vector<IO::Border_rel> borders;
	std::vector<IO::Corner_rel> corner_rel;

	std::vector<IO::Touches_rel> touches;
	std::vector<IO::Protrudes_rel> protrudes;
	std::vector<IO::Adjacent_rel> adjacent;

	// Transform Map::Center, Map::Corner, Map::Edge to IO structs, save their relations
	int centercount = centers.size();
	IO::Center* io_centers = new IO::Center[centercount];
	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		io_centers[(*it)->index].index = (*it)->index;
		io_centers[(*it)->index].point = (*it)->point;
		io_centers[(*it)->index].water = (*it)->water;
		io_centers[(*it)->index].ocean = (*it)->ocean;
		io_centers[(*it)->index].coast = (*it)->coast;
		io_centers[(*it)->index].border = (*it)->border;
		io_centers[(*it)->index].elevation = (*it)->elevation;
		io_centers[(*it)->index].moisture = (*it)->moisture;
		io_centers[(*it)->index].biome = (*it)->biome;

		for (std::set<Map::Center*>::iterator cit = (*it)->neighbours.begin(); cit != (*it)->neighbours.end(); cit++) {
			IO::Neighbour_rel n;
			n.center_out = (*it)->index;
			n.center_in = (*cit)->index;
			neighbours.push_back(n);
		}

		for (std::set<Map::Corner*>::iterator cit = (*it)->corners.begin(); cit != (*it)->corners.end(); cit++) {
			corners.insert(*cit);

			IO::Corner_rel c;
			c.center_out = (*it)->index;
			c.corner_in = (*cit)->index;
			corner_rel.push_back(c);
		}

		for (std::vector<Map::Edge*>::iterator eit = (*it)->borders.begin(); eit != (*it)->borders.end(); eit++) {
			edges.insert(*eit);

			IO::Border_rel b;
			b.center_out = (*it)->index;
			b.edge_in = (*eit)->index;
			borders.push_back(b);
		}
	}

	int cornercount = corners.size();
	IO::Corner* io_corners = new IO::Corner[cornercount];
	for (std::set<Map::Corner*>::iterator it = corners.begin(); it != corners.end(); it++) {
		io_corners[(*it)->index].index = (*it)->index;
		io_corners[(*it)->index].point = (*it)->point;
		io_corners[(*it)->index].water = (*it)->water;
		io_corners[(*it)->index].ocean = (*it)->ocean;
		io_corners[(*it)->index].coast = (*it)->coast;
		io_corners[(*it)->index].border = (*it)->border;
		io_corners[(*it)->index].elevation = (*it)->elevation;
		io_corners[(*it)->index].moisture = (*it)->moisture;
		io_corners[(*it)->index].watershed_size = (*it)->watershed_size;
		io_corners[(*it)->index].corner_downslope = (*it)->downslope->index;
		io_corners[(*it)->index].corner_watershed = (*it)->watershed->index;

		for (std::set<Map::Center*>::iterator cit = (*it)->touches.begin(); cit != (*it)->touches.end(); cit++) {
			IO::Touches_rel t;
			t.corner_out = (*it)->index;
			t.center_in = (*cit)->index;
			touches.push_back(t);
		}

		for (std::vector<Map::Edge*>::iterator eit = (*it)->protrudes.begin(); eit != (*it)->protrudes.end(); eit++) {
			IO::Protrudes_rel p;
			p.corner_out = (*it)->index;
			p.edge_in = (*eit)->index;
			protrudes.push_back(p);
		}

		for (std::set<Map::Corner*>::iterator cit = (*it)->adjacent.begin(); cit != (*it)->adjacent.end(); cit++) {
			IO::Adjacent_rel a;
			a.corner_out = (*it)->index;
			a.corner_in = (*cit)->index;
			adjacent.push_back(a);
		}
	}

	int edgecount = edges.size();
	IO::Edge* io_edges = new IO::Edge[edgecount];
	for (std::set<Map::Edge*>::iterator it = edges.begin(); it != edges.end(); it++) {
		io_edges[(*it)->index].index = (*it)->index;
		io_edges[(*it)->index].center_d0 = (*it)->d0->index;
		io_edges[(*it)->index].center_d1 = (*it)->d1->index;
		io_edges[(*it)->index].corner_v0 = (*it)->v0->index;
		io_edges[(*it)->index].corner_v1 = (*it)->v1->index;
		io_edges[(*it)->index].midway = (*it)->midway;
		io_edges[(*it)->index].river = (*it)->river;
	}

	// Make arrays from relations
	int neighbours_count = neighbours.size();
	IO::Neighbour_rel* neighbours_a = (IO::Neighbour_rel*) malloc(neighbours_count*sizeof(IO::Neighbour_rel));
	std::copy(neighbours.begin(), neighbours.end(), neighbours_a);

	int borders_count = borders.size();
	IO::Border_rel* borders_a = (IO::Border_rel*) malloc(borders_count*sizeof(IO::Border_rel));
	std::copy(borders.begin(), borders.end(), borders_a);

	int corner_rel_count = corner_rel.size();
	IO::Corner_rel* corner_rel_a = (IO::Corner_rel*) malloc(corner_rel_count*sizeof(IO::Corner_rel));
	std::copy(corner_rel.begin(), corner_rel.end(), corner_rel_a);

	int touches_count = touches.size();
	IO::Touches_rel* touches_a = (IO::Touches_rel*) malloc(touches_count*sizeof(IO::Touches_rel));
	std::copy(touches.begin(), touches.end(), touches_a);

	int protrudes_count = protrudes.size();
	IO::Protrudes_rel* protrudes_a = (IO::Protrudes_rel*) malloc(protrudes_count*sizeof(IO::Protrudes_rel));
	std::copy(protrudes.begin(), protrudes.end(), protrudes_a);

	int adjacent_count = adjacent.size();
	IO::Adjacent_rel* adjacent_a = (IO::Adjacent_rel*) malloc(adjacent_count*sizeof(IO::Adjacent_rel));
	std::copy(adjacent.begin(), adjacent.end(), adjacent_a);

	// Save to disk
	// Map info
	fwrite(&info, sizeof(Map::Info), 1, file);
	// Center, Corner, Edge
	fwrite(&centercount, sizeof(int), 1, file);
	fwrite(io_centers, sizeof(IO::Center), centercount, file);
	fwrite(&cornercount, sizeof(int), 1, file);
	fwrite(io_corners, sizeof(IO::Corner), cornercount, file);
	fwrite(&edgecount, sizeof(int), 1, file);
	fwrite(io_edges, sizeof(IO::Edge), edgecount, file);
	// Relationships
	fwrite(&neighbours_count, sizeof(int), 1, file);
	fwrite(neighbours_a, sizeof(IO::Neighbour_rel), neighbours_count, file);
	fwrite(&borders_count, sizeof(int), 1, file);
	fwrite(borders_a, sizeof(IO::Border_rel), borders_count, file);
	fwrite(&corner_rel_count, sizeof(int), 1, file);
	fwrite(corner_rel_a, sizeof(IO::Corner_rel), corner_rel_count, file);
	fwrite(&touches_count, sizeof(int), 1, file);
	fwrite(touches_a, sizeof(IO::Touches_rel), touches_count, file);
	fwrite(&protrudes_count, sizeof(int), 1, file);
	fwrite(protrudes_a, sizeof(IO::Protrudes_rel), protrudes_count, file);
	fwrite(&adjacent_count, sizeof(int), 1, file);
	fwrite(adjacent_a, sizeof(IO::Adjacent_rel), adjacent_count, file);
};

Map::Center *centerFromId(std::vector<Map::Center*> centers, int id) {
	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		if ((*it)->index == id) {
			return (*it);
		}
	}
	return NULL;
};

Map::Corner *cornerFromId(std::vector<Map::Corner*> corners, int id) {
	for (std::vector<Map::Corner*>::iterator it = corners.begin(); it != corners.end(); it++) {
		if ((*it)->index == id) {
			return (*it);
		}
	}
	return NULL;
};

Map::Edge *edgeFromId(std::vector<Map::Edge*> edge, int id) {
	for (std::vector<Map::Edge*>::iterator it = edge.begin(); it != edge.end(); it++) {
		if ((*it)->index == id) {
			return (*it);
		}
	}
};

std::vector<Map::Center*> IO::importMap(FILE *file, Map::Info* info) {
	int centercount;
	int cornercount;
	int edgecount;

	int neighbours_count;
	int borders_count;
	int corner_rel_count;
	int touches_count;
	int protrudes_count;
	int adjacent_count;

	fread(info, sizeof(Map::Info), 1, file);

	fread(&centercount, sizeof(int), 1, file);
	IO::Center* io_centers = (IO::Center*) malloc(centercount*sizeof(IO::Center));
	fread(io_centers, sizeof(IO::Center), centercount, file);

	fread(&cornercount, sizeof(int), 1, file);
	IO::Corner* io_corners = (IO::Corner*) malloc(cornercount*sizeof(IO::Corner));
	fread(io_corners, sizeof(IO::Corner), cornercount, file);

	fread(&edgecount, sizeof(int), 1, file);
	IO::Edge* io_edges = (IO::Edge*) malloc(edgecount*sizeof(IO::Edge));
	fread(io_edges, sizeof(IO::Edge), edgecount, file);


	fread(&neighbours_count, sizeof(int), 1, file);
	IO::Neighbour_rel* neighbours_a = (IO::Neighbour_rel*) malloc(neighbours_count*sizeof(IO::Neighbour_rel));
	fread(neighbours_a, sizeof(IO::Neighbour_rel), neighbours_count, file);

	fread(&borders_count, sizeof(int), 1, file);
	IO::Border_rel* borders_a = (IO::Border_rel*) malloc(borders_count*sizeof(IO::Border_rel));
	fread(borders_a, sizeof(IO::Border_rel), borders_count, file);

	fread(&corner_rel_count, sizeof(int), 1, file);
	IO::Corner_rel* corner_rel_a = (IO::Corner_rel*) malloc(corner_rel_count*sizeof(IO::Corner_rel));
	fread(corner_rel_a, sizeof(IO::Corner_rel), corner_rel_count, file);

	fread(&touches_count, sizeof(int), 1, file);
	IO::Touches_rel* touches_a = (IO::Touches_rel*) malloc(touches_count*sizeof(IO::Touches_rel));
	fread(touches_a, sizeof(IO::Touches_rel), touches_count, file);

	fread(&protrudes_count, sizeof(int), 1, file);
	IO::Protrudes_rel* protrudes_a = (IO::Protrudes_rel*) malloc(protrudes_count*sizeof(IO::Protrudes_rel));
	fread(protrudes_a, sizeof(IO::Protrudes_rel), protrudes_count, file);

	fread(&adjacent_count, sizeof(int), 1, file);
	IO::Adjacent_rel* adjacent_a = (IO::Adjacent_rel*) malloc(adjacent_count*sizeof(IO::Adjacent_rel));
	fread(adjacent_a, sizeof(IO::Adjacent_rel), adjacent_count, file);

	std::map<int, Map::Center*> centers;
	std::map<int, Map::Corner*> corners;
	std::map<int, Map::Edge*> edges;

	//Rebuild Map data structures

	for (int i = 0; i < centercount; i++) {
		Map::Center *c = new Map::Center(io_centers[i].index, io_centers[i].point);
		c->water = io_centers[i].water;
		c->ocean = io_centers[i].ocean;
		c->coast = io_centers[i].coast;
		c->border = io_centers[i].border;
		c->elevation = io_centers[i].elevation;
		c->moisture = io_centers[i].moisture;
		c->biome = io_centers[i].biome;

		centers.insert(std::pair<int, Map::Center*>(io_centers[i].index, c));
	}

	for (int i = 0; i < cornercount; i++) {
		Map::Corner *c = new Map::Corner(io_corners[i].index, io_corners[i].point);
		c->water = io_corners[i].water;
		c->ocean = io_corners[i].ocean;
		c->coast = io_corners[i].coast;
		c->border = io_corners[i].border;
		c->elevation = io_corners[i].elevation;
		c->moisture = io_corners[i].moisture;
		c->watershed_size = io_corners[i].watershed_size;

		corners.insert(std::pair<int, Map::Corner*>(io_corners[i].index, c));
	}

	for (int i = 0; i < edgecount; i++) {
		Map::Edge *e = new Map::Edge(
			io_edges[i].index,
			centers[io_edges[i].center_d0],
			centers[io_edges[i].center_d1],
			corners[io_edges[i].corner_v0],
			corners[io_edges[i].corner_v1],
			io_edges[i].midway
		);
		e->river = io_edges[i].river;
		edges.insert(std::pair<int, Map::Edge*>(io_edges[i].index, e));
	}

	//Rebuild one to one relations
	for (int i = 0; i < cornercount; i++) {
		corners[io_corners[i].index]->downslope = corners[io_corners[i].corner_downslope];
		corners[io_corners[i].index]->watershed = corners[corners, io_corners[i].corner_watershed];
	}

	//Rebuild relations
	//The structure here is always the same.
	//We use the *FromId functions to get a pointer to the out node,
	//Then we add the in node to the relation structure on the out node

	for (int i = 0; i < neighbours_count; i++) {
		Map::Center* out = centers[neighbours_a[i].center_out];
		if (out != NULL) out->neighbours.insert(centers[neighbours_a[i].center_in]);
	}

	for (int i = 0; i < borders_count; i++) {
		Map::Center* out = centers[borders_a[i].center_out];
		if (out != NULL) out->borders.push_back(edges[borders_a[i].edge_in]);
	}

	for (int i = 0; i < corner_rel_count; i++) {
		Map::Center* out = centers[corner_rel_a[i].center_out];
		if (out != NULL) out->corners.insert(corners[corner_rel_a[i].corner_in]);
	}

	for (int i = 0; i < touches_count; i++) {
		Map::Corner* out = corners[touches_a[i].corner_out];
		if (out != NULL) out->touches.insert(centers[touches_a[i].center_in]);
	}

	for (int i = 0; i < protrudes_count; i++) {
		Map::Corner* out = corners[protrudes_a[i].corner_out];
		if (out != NULL) out->protrudes.push_back(edges[protrudes_a[i].edge_in]);
	}

	for (int i = 0; i < adjacent_count; i++) {
		Map::Corner* out = corners[adjacent_a[i].corner_out];
		if (out != NULL) out->adjacent.insert(corners[adjacent_a[i].corner_in]);
	}

	std::vector<Map::Center*> ret;

	for (std::map<int, Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		ret.push_back(&(*it->second));
	}

	return ret;
};

