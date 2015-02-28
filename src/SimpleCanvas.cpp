#include "SimpleCanvas.h"

#include <wx/dcclient.h>

#include "Generator.h"

wxBEGIN_EVENT_TABLE(SimpleCanvas, wxGLCanvas)
	EVT_PAINT(SimpleCanvas::Paint)
wxEND_EVENT_TABLE()

void SimpleCanvas::GenerateGeometry() {
	g = new Graph();
	printf("Starting gen\n");
	Generator *gen = new Generator(600, 600, 2000);
	gen->setPolygonGraph(g);

	gen->start();
	printf("Finished gen\n");
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

	int edges = g->getEdgeCount();
	float *edge = new float[4];

	glBegin(GL_LINES);
		glColor3f(0.0, 0.0, 0.0);

		for (int i = 0; i < edges; i++) {
			edge = g->getEdge(i);

			glVertex2f((edge[0]/600.0f)*2 - 1.0f, (edge[2]/600.0f)*2 - 1.0f);
			glVertex2f((edge[1]/600.0f)*2 - 1.0f, (edge[3]/600.0f)*2 - 1.0f);
		}
	glEnd();


	int nodes = g->getNodeCount();
	float *node = new float[2];

	glBegin(GL_POINTS);
		glPointSize(0.002f);
		for (int i = 0; i< nodes; i++) {
			node = g->getNode(i);
			glVertex2f((node[0]/600.0f)*2 - 1.0f, (node[1]/600.0f)*2 - 1.0f);
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
