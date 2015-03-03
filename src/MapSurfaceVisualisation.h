#include "Generator.h"

class MapSurfaceCellVertices : public Attribute {
public:
	MapSurfaceCellVertices(Map::Center* center, int size);
};

class MapSurface {
	int cellcount;
	ColoredObject** cells;
public:
	MapSurface(std::vector<Map::Center*> centers);
	void draw();
};
