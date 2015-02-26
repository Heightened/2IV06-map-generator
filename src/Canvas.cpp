#include "Canvas.h"
#include "Generator.h"

#include <wx/dcclient.h>

#include "GraphVisualisation.h"

GLuint initShaders(const char * VertexShaderFile, const char * FragmentShaderFile) {
	const GLubyte * a = glGetString(GL_SHADING_LANGUAGE_VERSION_ARB);
	//wxLogError(wxT("%i\n"), a);

	//Line string to buffer code of shaders
	std::string Line = "";

	// Read the Vertex Shader code from the file
	std::string VShaderSource;
	std::ifstream VShaderStream(VertexShaderFile);
	while(getline(VShaderStream, Line)) {
		VShaderSource += "\n" + Line;
	}
	VShaderStream.close();

	//wxLogError(wxT("%s\n"),VShaderSource);

	// Read the Fragment Shader code from the file
	std::string FShaderSource;
	std::ifstream FShaderStream(FragmentShaderFile, std::ios::in);
	while(getline(FShaderStream, Line)) {
		FShaderSource += "\n" + Line;
	}
	FShaderStream.close();

	//wxLogError(wxT("%s\n"),FShaderSource);

	// Create the vertex shader
	GLuint VShader = glCreateShader(GL_VERTEX_SHADER);
	char const * VShaderSourcePtr = VShaderSource.c_str();
	glShaderSource(VShader, 1, &VShaderSourcePtr, NULL);
	glCompileShader(VShader);

	GLint Result = GL_FALSE;
	int InfoLogLength;

	// Check Vertex Shader
	glGetShaderiv(VShader, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(VShader, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 0 ){
		std::vector<char> VShaderErrorMessage(InfoLogLength+1);
		glGetShaderInfoLog(VShader, InfoLogLength, NULL, &VShaderErrorMessage[0]);
		wxLogError(wxT("%s\n"), &VShaderErrorMessage[0]);
	}

	// Create fragment Shader
	GLuint FShader = glCreateShader(GL_FRAGMENT_SHADER);
	char const * FShaderSourcePtr = FShaderSource.c_str();
	glShaderSource(FShader, 1, &FShaderSourcePtr, NULL);
	glCompileShader(FShader);

	// Check Fragment Shader
	glGetShaderiv(FShader, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(FShader, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 0 ){
		std::vector<char> FragmentShaderErrorMessage(InfoLogLength+1);
		glGetShaderInfoLog(FShader, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
		wxLogError(wxT("%s\n"), &FragmentShaderErrorMessage[0]);
	}

	// Link the program
	GLuint Program = glCreateProgram();
	glAttachShader(Program, VShader);
	glAttachShader(Program, FShader);
	glLinkProgram(Program);

	glDeleteShader(VShader);
	glDeleteShader(FShader);

	// Check the program
	glGetProgramiv(Program, GL_LINK_STATUS, &Result);
	glGetProgramiv(Program, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 0 ){
		std::vector<char> ProgramErrorMessage(InfoLogLength+1);
		glGetProgramInfoLog(Program, InfoLogLength, NULL, &ProgramErrorMessage[0]);
		wxLogError(wxT("%s\n"), &ProgramErrorMessage[0]);
	}

	return Program;
}

GLuint initShaders() {
	return initShaders("PreviewVertexShader.glsl", "PreviewFragmentShader.glsl");
}

//Class to maintain matrices used in shaders.

void UniformMatrices::setFov(float _fov) {
	projection = glm::perspective(_fov, aspect, 0.1f, 100.0f);
}

void UniformMatrices::setWindowSize(int width, int height) {
	aspect = width / (float)height;
	projection = glm::perspective(fov, aspect, 0.1f, 100.0f);
}

void UniformMatrices::setPosition(glm::vec3 _position) {
	position = _position;
	view = glm::lookAt(position, focus, glm::vec3(0,0,1));
}

void UniformMatrices::movePosition(glm::vec3 direction) {
	position += direction;
	view = glm::lookAt(position, focus, glm::vec3(0,0,1));
}

void UniformMatrices::movePositionRelative(glm::vec3 direction) {
	//TODO relative motion in it's own space.
	position += direction;
	view = glm::lookAt(position, focus, glm::vec3(0,0,1));
}

glm::vec3 UniformMatrices::getPosition() {
	return position;
}

void UniformMatrices::setFocus(glm::vec3 _focus) {
	focus = _focus;
	view = glm::lookAt(position, focus, glm::vec3(0,0,1));
}

void UniformMatrices::moveFocus(glm::vec3 direction) {
	focus += direction;
	view = glm::lookAt(position, focus, glm::vec3(0,0,1));
}

void UniformMatrices::aim(glm::vec3 direction) {
	focus = position + direction;
	view = glm::lookAt(position, focus, glm::vec3(0,0,1));
}

glm::vec3 UniformMatrices::getFocus() {
	return focus;
}

void UniformMatrices::move(glm::vec3 direction) {
	movePosition(direction);
	moveFocus(direction);
}

void UniformMatrices::moveRelative(glm::vec3 direction) {
	//TODO relative motion in it's own space.
	movePosition(direction);
	moveFocus(direction);
}

void UniformMatrices::send(glm::mat4 _model) {
	model = _model;

	combined = projection * view * _model;

	glUniformMatrix4fv(combinedMatrixId, 1, GL_FALSE, &combined[0][0]);
	glUniformMatrix4fv(modelId, 1, GL_FALSE, &model[0][0]);
	glUniformMatrix4fv(viewId, 1, GL_FALSE, &view[0][0]);
}

UniformMatrices::UniformMatrices(GLuint shaders, glm::mat4 projection, glm::vec3 _position, glm::vec3 _focus, glm::mat4 model, const char * combinedName, const char * viewName, const char * modelName) {
	if (glIsProgram(shaders)) wxLogError(wxT("Shaders are valid\n"));

	combinedMatrixId = glGetUniformLocation(shaders, combinedName);
	viewId = glGetUniformLocation(shaders, viewName);
	modelId = glGetUniformLocation(shaders, modelName);

	wxLogError(wxT("combinedMatrixId: %i \n"), combinedMatrixId);
	wxLogError(wxT("viewId: %i \n"), viewId);
	wxLogError(wxT("modelId: %i \n"), modelId);

	this->projection = projection;
	this->view = view;

	position = _position;
	focus = _focus;
	this->view = glm::lookAt(position, focus, glm::vec3(0,0,1));

	send(model);
}

//wxGLCanvas subclass

static void CheckGLError() {
    GLenum errLast = GL_NO_ERROR;

    for ( ;; )
    {
        GLenum err = glGetError();
        if ( err == GL_NO_ERROR )
            return;

        // normally the error is reset by the call to glGetError() but if
        // glGetError() itself returns an error, we risk looping forever here
        // so check that we get a different error than the last time
        if ( err == errLast )
        {
            wxLogError(wxT("OpenGL error state couldn't be reset."));
            return;
        }

        errLast = err;

        wxLogError(wxT("OpenGL error %d"), err);
    }
}

wxBEGIN_EVENT_TABLE(Canvas, wxGLCanvas)
	EVT_PAINT(Canvas::Paint)
wxEND_EVENT_TABLE()

Canvas::Canvas(wxWindow* parent, wxSize size) : wxGLCanvas(parent, wxID_ANY, NULL, wxDefaultPosition, size, 0, "Preview", wxNullPalette) {
	this->size = size;
}

void Canvas::GenerateGeometry() {
	Generator *gen = new Generator(600, 600, 2000);

	Graph *g = gen->start();

	int edges = g->getEdgeCount();
	int nodes = g->getNodeCount();
	GraphVertices vertices(g, edges, nodes, 3.0f, 0.2f);
	//SphereVertices vertices(1.0f, 3); // 0,1,2 or 3 are the only values for subdivisions, higher values will cause exceptions
	//BoxVertices vertices(2.0f, 2.0f, 2.0f);
	Normals normals(vertices.size()/3, vertices, true);
	SolidColor color(vertices.size()/3, 0.8f, 0.5f, 0.0f);
	const GLfloat* attributes[] = {vertices, normals, color};
	object = new ColoredObject(vertices.size()/3, &attributes[0]);
}

void Canvas::Paint(wxPaintEvent& WXUNUSED(event)) {
	wxPaintDC dc(this);
	
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	int width = size.GetWidth();
	int height = size.GetHeight();
	glViewport(0, 0, width, height);
	glCullFace(GL_BACK);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //glUseProgram(shaders);
	
	//Draw geometry
	glm::mat4 Matrix(1,0,0,0,
						0,1,0,0,
						0,0,1,0,
						0,0,0,1);
	for(int i = 0; i < 1; i++) {
		viewer->send(Matrix);
		object->draw();
		Matrix = glm::translate(Matrix, glm::vec3(3,0,0));
	}
	
	viewer->setWindowSize(width, height);
	viewer->setFov(60.0f);

	SwapBuffers();
}

void Canvas::Initialize(wxGLContext* context) {
	glewExperimental = true;
	SetCurrent(*context);
	bool glew = glewInit();
	if (glew != GLEW_OK) {
		exit(EXIT_FAILURE);
	}
	
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glEnable(GL_CULL_FACE);
	glClearColor(0,0,0,1.0f);
	
	GLuint shaders = initShaders();
	glUseProgram(shaders);
	viewer = new UniformMatrices(shaders);
	
	GenerateGeometry();

	glFlush();

	CheckGLError();
}
