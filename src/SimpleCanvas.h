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
#include "Generator.h"
//TODO: Seperate data structure definition from generator

class SimpleCanvas : public Canvas {
	std::vector<Map::Center> centers;
	wxDECLARE_EVENT_TABLE();
public:
	SimpleCanvas(wxWindow* parent, wxSize size, Generator *gen): Canvas(parent, size, gen){};
	virtual void GenerateGeometry();
	virtual void Paint(wxPaintEvent& WXUNUSED(event));
	virtual void Initialize();
};
