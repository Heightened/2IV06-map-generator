#include <cstdio>
#include <algorithm>

#include "MapIO.h"

void IO::exportMap(FILE *file, std::vector<Map::Center*> centers) {
	std::set<Map::Corner*> corners;
	std::set<Map::Edge> edges;

	// Relationships
	std::vector<IO::Neighbour_rel> neighbours;
	std::vector<IO::Border_rel> borders;
	std::vector<IO::Corner_rel> corner_rel;

	std::vector<IO::Touches_rel> touches;
	std::vector<IO::Protrudes_rel> protrudes;
	std::vector<IO::Adjacent_rel> adjacent;

	// Transform Map::Center, Map::Corner, Map::Edge to IO structs, save their relations
	int centercount = centers.size();
	IO::Center io_centers[centercount];
	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		io_centers[(*it)->index].index = (*it)->index;
		io_centers[(*it)->index].point = (*it)->point;
		io_centers[(*it)->index].water = (*it)->water;
		io_centers[(*it)->index].ocean = (*it)->ocean;
		io_centers[(*it)->index].coast = (*it)->coast;
		io_centers[(*it)->index].border = (*it)->border;
		io_centers[(*it)->index].elevation = (*it)->elevation;

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

		for (std::vector<Map::Edge>::iterator eit = (*it)->borders.begin(); eit != (*it)->borders.end(); eit++) {
			edges.insert(*eit);

			IO::Border_rel b;
			b.center_out = (*it)->index;
			b.edge_in = eit->index;
			borders.push_back(b);
		}
	}

	int cornercount = corners.size();
	IO::Corner io_corners[cornercount];
	for (std::set<Map::Corner*>::iterator it = corners.begin(); it != corners.end(); it++) {
		io_centers[(*it)->index].index = (*it)->index;
		io_centers[(*it)->index].point = (*it)->point;
		io_centers[(*it)->index].water = (*it)->water;
		io_centers[(*it)->index].ocean = (*it)->ocean;
		io_centers[(*it)->index].coast = (*it)->coast;
		io_centers[(*it)->index].border = (*it)->border;
		io_centers[(*it)->index].elevation = (*it)->elevation;

		for (std::set<Map::Center*>::iterator cit = (*it)->touches.begin(); cit != (*it)->touches.end(); cit++) {
			IO::Touches_rel t;
			t.corner_out = (*it)->index;
			t.center_in = (*cit)->index;
			touches.push_back(t);
		}

		for (std::vector<Map::Edge>::iterator eit = (*it)->protrudes.begin(); eit != (*it)->protrudes.end(); eit++) {
			IO::Protrudes_rel p;
			p.corner_out = (*it)->index;
			p.edge_in = eit->index;
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
	IO::Edge io_edges[edgecount];
	for (std::set<Map::Edge>::iterator it = edges.begin(); it != edges.end(); it++) {
		io_edges[it->index].index = it->index;
		io_edges[it->index].center_d0 = it->d0->index;
		io_edges[it->index].center_d1 = it->d1->index;
		io_edges[it->index].corner_v0 = it->v0->index;
		io_edges[it->index].corner_v1 = it->v1->index;
		io_edges[it->index].midway = it->midway;
	}

	// Make arrays from relations
	int neighbours_count = neighbours.size();
	IO::Neighbour_rel neighbours_a[neighbours_count];
	std::copy(neighbours.begin(), neighbours.end(), neighbours_a);

	int borders_count = borders.size();
	IO::Border_rel borders_a[borders_count];
	std::copy(borders.begin(), borders.end(), borders_a);

	int corner_rel_count = corners.size();
	IO::Corner_rel corner_rel_a[corner_rel_count];
	std::copy(corner_rel.begin(), corner_rel.end(), corner_rel_a);

	int touches_count = touches.size();
	IO::Touches_rel touches_a[touches_count];
	std::copy(touches.begin(), touches.end(), touches_a);

	int protrudes_count = protrudes.size();
	IO::Protrudes_rel protrudes_a[protrudes_count];
	std::copy(protrudes.begin(), protrudes.end(), protrudes_a);

	int adjacent_count = adjacent.size();
	IO::Adjacent_rel adjacent_a[adjacent_count];
	std::copy(adjacent.begin(), adjacent.end(), adjacent_a);

	// Save to disk
	// Center, Corner, Edge
	fwrite(&centercount, sizeof(int), 1, file);
	fwrite(&io_centers, sizeof(IO::Center), centercount, file);
	fwrite(&cornercount, sizeof(int), 1, file);
	fwrite(&io_corners, sizeof(IO::Corner), cornercount, file);
	fwrite(&edgecount, sizeof(int), 1, file);
	fwrite(&io_edges, sizeof(IO::Edge), edgecount, file);
	// Relationships
	fwrite(&neighbours_count, sizeof(int), 1, file);
	fwrite(&neighbours_a, sizeof(IO::Neighbour_rel), neighbours_count, file);
	fwrite(&borders_count, sizeof(int), 1, file);
	fwrite(&borders_a, sizeof(IO::Border_rel), borders_count, file);
	fwrite(&corner_rel_count, sizeof(int), 1, file);
	fwrite(&corner_rel_a, sizeof(IO::Corner_rel), corner_rel_count, file);
	fwrite(&touches_count, sizeof(int), 1, file);
	fwrite(&touches_a, sizeof(IO::Touches_rel), touches_count, file);
	fwrite(&protrudes_count, sizeof(int), 1, file);
	fwrite(&protrudes_a, sizeof(IO::Protrudes_rel), protrudes_count, file);
	fwrite(&adjacent_count, sizeof(int), 1, file);
	fwrite(&adjacent_a, sizeof(IO::Adjacent_rel), adjacent_count, file);
};

std::vector<Map::Center*> IO::importMap(FILE *file) {

};
