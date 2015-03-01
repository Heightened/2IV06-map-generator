#include "SimpleCanvas.h"

#include <wx/dcclient.h>

#include "Generator.h"

wxBEGIN_EVENT_TABLE(SimpleCanvas, wxGLCanvas)
	EVT_PAINT(SimpleCanvas::Paint)
wxEND_EVENT_TABLE()

#define NORMALIZE(point, width, height) (point.x/width)*2 - 1.0f, (point.y/height)*2 - 1.0f

void SimpleCanvas::GenerateGeometry() {
	centers = gen->getCenters();
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

	glColor3f(0.0, 0.0, 0.0);

	glBegin(GL_LINES);
	for (std::vector<Map::Center*>::iterator it = centers.begin(); it != centers.end(); it++) {
		for (std::vector<Map::Edge>::iterator eit = (*it)->borders.begin(); eit != (*it)->borders.end(); eit++) {
			glVertex2f(NORMALIZE(eit->v0->point, width, height));
			glVertex2f(NORMALIZE(eit->v1->point, width, height));
		}
	}
	glEnd();

	glBegin(GL_POINTS);
		glPointSize(0.002f);
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


	GenerateGeometry();

	glFlush();

	init = true;
}
