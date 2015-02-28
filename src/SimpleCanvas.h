#pragma once

#include <GL/glew.h>

#include <wx/glcanvas.h>


#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <string>
#include <fstream>
#include <vector>
#include <iostream>

#include "Canvas.h"
#include "GraphVisualisation.h"

class SimpleCanvas : public Canvas {
	Graph *g;
	wxDECLARE_EVENT_TABLE();
public:
	SimpleCanvas(wxWindow* parent, wxSize size, Generator *gen): Canvas(parent, size, gen){};
	virtual void GenerateGeometry();
	virtual void Paint(wxPaintEvent& WXUNUSED(event));
	virtual void Initialize();
};
