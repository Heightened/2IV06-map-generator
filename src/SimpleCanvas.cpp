#include "SimpleCanvas.h"

#include <algorithm>
#include <wx/dcclient.h>

#include "Generator.h"

wxBEGIN_EVENT_TABLE(SimpleCanvas, wxGLCanvas)
	EVT_PAINT(SimpleCanvas::Paint)
wxEND_EVENT_TABLE()

#define NORMALIZE(point, width, height) (point.x/width)*2 - 1.0f, (point.y/height)*2 - 1.0f
#define rgb(r,g,b) r/255.0f, g/255.0f, b/255.0f

class CounterClockwiseCompare {
	glm::vec2 center;
	public:
		CounterClockwiseCompare(glm::vec2 c): center(c) {}
		bool operator()(glm::vec2 const& a, glm::vec2 const& b) {
			if (a.x - center.x >= 0 && b.x - center.x < 0) {
				return false;
			}
			if (a.x - center.x < 0 && b.x - center.x >= 0) {
				return true;
			}
			if (a.x - center.x == 0 && b.x - center.x == 0) {
				if (a.y - center.y >= 0 || b.y - center.y >= 0) {
					return a.y < b.y;
				}
				return b.y < a.y;
			}

			// compute the cross product of vectors (center -> a) x (center -> b)
			int det = (a.x - center.x) * (b.y - center.y) - (b.x - center.x) * (a.y - center.y);
			if (det < 0) {
				return false;
			}
			if (det > 0) {
				return true;
			}

			// points a and b are on the same line from the center
			// check which point is closer to the center
			int d1 = (a.x - center.x) * (a.x - center.x) + (a.y - center.y) * (a.y - center.y);
			int d2 = (b.x - center.x) * (b.x - center.x) + (b.y - center.y) * (b.y - center.y);
			return d1 < d2;
		}
};

static void setBiome(Map::Biome b, float alpha) {
	switch(b) {
	case Map::LAKE:
		glColor4f(rgb(91, 132, 173), alpha);
		break;
	case Map::BEACH:
		glColor4f(rgb(172,159,139), alpha);
		break;
	case Map::ICE:
		glColor4f(rgb(222, 230, 245), alpha);
		break;
	case Map::MARSH:
		glColor4f(rgb(64, 108, 86), alpha);
		break;
	case Map::SNOW:
		glColor4f(rgb(248, 248, 248), alpha);
		break;
	case Map::TUNDRA:
		glColor4f(rgb(221, 221, 187), alpha);
		break;
	case Map::BARE:
		glColor4f(rgb(187, 187, 187), alpha);
		break;
	case Map::SCORCHED:
		glColor4f(rgb(153, 153, 153), alpha);
		break;
	case Map::TAIGA:
		glColor4f(rgb(204, 212, 187), alpha);
		break;
	case Map::SHRUBLAND:
		glColor4f(rgb(196, 204, 187), alpha);
		break;
	case Map::TEMPERATE_DESERT:
		glColor4f(rgb(228, 232, 202), alpha);
		break;
	case Map::TEMPERATE_RAIN_FOREST:
		glColor4f(rgb(164, 196, 168), alpha);
		break;
	case Map::TEMPERATE_DECIDUOUS_FOREST:
		glColor4f(rgb(180, 201, 169), alpha);
		break;
	case Map::GRASSLAND:
		glColor4f(rgb(196, 212, 170), alpha);
		break;
	case Map::TROPICAL_RAIN_FOREST:
		glColor4f(rgb(156, 187, 169), alpha);
		break;
	case Map::TROPICAL_SEASONAL_FOREST:
		glColor4f(rgb(169, 204, 164), alpha);
		break;
	case Map::SUBTROPICAL_DESERT:
		glColor4f(rgb(233, 221, 199), alpha);
		break;
	case Map::OCEAN:
	default:
		glColor4f(rgb(54, 54, 97), alpha);
	}
}

void SimpleCanvas::GenerateGeometry() {
	centers = gen->getCenters();
}

void SimpleCanvas::GenerateGeometry(std::vector<Map::Center*> _centers) {
	centers = _centers;
}

void SimpleCanvas::Paint(wxPaintEvent& WXUNUSED(event)) {
	if (!init) {
		Initialize();
	}

	wxPaintDC(this);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);
	GLint xsize = (GLint)GetSize().x;
	GLint ysize = (GLint)GetSize().y;
	glViewport(0, 0, xsize, ysize);

	float width = gen->getWidth();
	float height = gen->getHeight();

	float minElevation = FLT_MAX;
	float maxElevation = FLT_MIN;

	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		if ((*it)->elevation < minElevation) {
			minElevation = (*it)->elevation;
		}
		if ((*it)->elevation > maxElevation) {
			maxElevation = (*it)->elevation;
		}
	}

	glColor3f(0.0, 0.0, 0.0);

	float alpha;
	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		alpha = 1.0f - ((*it)->elevation - minElevation)/(maxElevation-minElevation)*0.7f;

		// Alpha values for height don't really work well with biomes, so disable them for now
		setBiome((*it)->biome, 1.0f);
		std::vector<glm::vec2> cornerPoints;
		for (std::set<Map::Corner*>::iterator eit = (*it)->corners.begin(); eit != (*it)->corners.end(); eit++) {
			cornerPoints.push_back((*eit)->point);
		}
		std::sort(cornerPoints.begin(), cornerPoints.end(), CounterClockwiseCompare((*it)->point));

		glBegin(GL_POLYGON);
		for (std::vector<glm::vec2>::iterator eit = cornerPoints.begin(); eit != cornerPoints.end(); eit++) {
			glVertex2f(NORMALIZE((*eit), width, height));
		}
		glEnd();

		glColor3f(34/255.0f, 85/255.0f, 136/255.0f);
		for (std::vector<Map::Edge*>::iterator eit = (*it)->borders.begin(); eit != (*it)->borders.end(); eit++) {
			if ((*eit)->river > 0) {
				glLineWidth(sqrt((*eit)->river));
				glBegin(GL_LINES);
				glVertex2f(NORMALIZE((*eit)->v0->point, width, height));
				glVertex2f(NORMALIZE((*eit)->v1->point, width, height));
				glEnd();
			}
		}
		glLineWidth(1);
	}

	glFlush();
	SwapBuffers();
}

void SimpleCanvas::Initialize() {
	SetCurrent();
	glewExperimental = true;
	bool glew = glewInit();
	if (glew != GLEW_OK) {
		exit(EXIT_FAILURE);
	}

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glEnable(GL_CULL_FACE);
	glClearColor(0,0,0,1.0f);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glHint(GL_LINE_SMOOTH, GL_NICEST);
	glEnable(GL_LINE_SMOOTH);

	GenerateGeometry();

	glFlush();

	init = true;
}
