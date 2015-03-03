#include <set>

#include "MapSurfaceVisualisation.h"

#include "wx/log.h"

typedef std::vector<Map::Edge, std::allocator<Map::Edge>> BORDERLIST;
typedef std::set<Map::Corner*, std::less<Map::Corner*>, std::allocator<Map::Corner*>> CORNERSET;

MapSurfaceCellVertices::MapSurfaceCellVertices(Map::Center* center, int size) : Attribute(size * 3 * 3) {
	Attribute vertices = Attribute(0, 0);
	
	glm::vec3 c(center->point, center->elevation);
	wxLogError(wxT("elevation: %f, size: %i"), center->elevation, size);

	for(BORDERLIST::iterator it = center->borders.begin(); it != center->borders.end(); ++it) {
		glm::vec3 a((*it).v0->point, (*it).v0->elevation);
		glm::vec3 b((*it).v1->point, (*it).v1->elevation);
		
		//The correctly ordered face will have a normal with positive z
		if (glm::cross(a-c, b-c).z < 0) {
			glm::vec3 temp = a;
			a = b;
			b = temp;
		}

		GLfloat triangle[] = {
			a.x, a.y, a.z,
			b.x, b.y, b.z,
			c.x, c.y, c.z
		};

		vertices = vertices + Attribute(3 * 3, &triangle[0]);
	}

	vertices.scale(1/15.0f);

/*
	CORNERSET::iterator itb = --center->corners.end();
	for (CORNERSET::iterator ita = center->corners.begin(); ita != center->corners.end(); ++ita) {
		wxLogError(wxT("add a triangle at (%f, %f, %f)"), c.x, c.y, c.z);
		
		glm::vec3 b((*ita)->point, (*ita)->elevation);
		glm::vec3 a((*itb)->point, (*itb)->elevation);
		
		GLfloat triangle[] = {
			a.x, a.y, a.z,
			b.x, b.y, b.z,
			c.x, c.y, c.z
		};

		vertices = vertices + Attribute(3 * 3, &triangle[0]);

		itb = ita;
	}
*/

	values = new GLfloat[length];
	for (int i = 0; i < length; i++) {
		values[i] = vertices[i];
	}
}

MapSurface::MapSurface(std::vector<Map::Center*> centers) {
	cellcount = centers.size();

	cells = new ColoredObject*[cellcount];

	int i = 0;
	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); ++it) {
		MapSurfaceCellVertices vertices((*it), (*it)->borders.size());
		Normals normals(vertices.size()/3, vertices, false);

		glm::vec3 c(0.2f, 0.8f, 0.2f);
		if ((*it)->ocean) {
			c = glm::vec3(0.1f, 0.1f, 0.3f);
		} else if ((*it)->water) {
			c = glm::vec3(0.2f, 0.2f, 0.6f);
		} else if ((*it)->coast) {
			c = glm::vec3(0.6f, 0.8f, 0.2f);
		}

		SolidColor color(vertices.size()/3, c.r, c.g, c.b);
		const GLfloat* attributes[] = {vertices, normals, color};
		cells[i] = new ColoredObject(vertices.size()/3, &attributes[0]);
		
		i++;
	}
}

void MapSurface::draw() {
	for (int i = 0; i < cellcount; i++) {
		cells[i]->draw();
	}
}
