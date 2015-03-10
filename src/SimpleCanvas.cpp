#include "SimpleCanvas.h"

#include <algorithm>
#include <wx/dcclient.h>

#include "Generator.h"

wxBEGIN_EVENT_TABLE(SimpleCanvas, wxGLCanvas)
	EVT_PAINT(SimpleCanvas::Paint)
wxEND_EVENT_TABLE()

#define NORMALIZE(point, width, height) (point.x/width)*2 - 1.0f, (point.y/height)*2 - 1.0f

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

		glColor4f(179/255.0f, 166/255.0f, 146/255.0f, alpha);
		if ((*it)->water) {
			glColor4f(91/255.0f, 132/255.0f, 173/255.0f, alpha);
		}
		if ((*it)->ocean) {
			glColor4f(54/255.0f, 54/255.0f, 97/255.0f, alpha);
		}
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

		glColor3f(0.0, 0.0, 0.0);
		glBegin(GL_LINE_LOOP);
		for (std::vector<glm::vec2>::iterator eit = cornerPoints.begin(); eit != cornerPoints.end(); eit++) {
			glVertex2f(NORMALIZE((*eit), width, height));
		}
		glEnd();
	}

	glColor3f(1.0, 0.0, 0.0);

	glPointSize(2.0f);
	glBegin(GL_POINTS);
	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		glVertex2f(NORMALIZE((*it)->point, width, height));
	}
	glEnd();

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

	GenerateGeometry();

	glFlush();

	init = true;
}
