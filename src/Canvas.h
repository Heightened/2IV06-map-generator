#include <GL/glew.h>

#include <wx/glcanvas.h>


#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <string>
#include <fstream>
#include <vector>
#include <iostream>

#include "Objects.h"

GLuint initShaders(const char * VertexShaderFile, const char * FragmentShaderFile);

GLuint initShaders();

//Uniform matrices
class UniformMatrices {
	//Perspective matrix should remain unchanged unless fov is adjusted
	glm::mat4 projection;
	//View matrix should be updated whenever the camera is moved
	glm::mat4 view;
	//The modelview matrix is adjusted to apply transformations to existing vertices 
	glm::mat4 model;

	//This matrix combines the perspective, model and view matrices to acquire the on-screen position of a vertex.
	//It is only updated when the matrices are sent to the shaders.
	glm::mat4 combined;

	//Uniform identifiers for matrices to be used in shaders.
	GLuint combinedMatrixId, viewId, modelId;

	//Reusable parameters
	float fov, aspect;
	glm::vec3 position, focus;

public:
	void setFov(float fov);
	void setWindowSize(int width, int height);

	void setPosition(glm::vec3 position);
	void movePosition(glm::vec3 direction);
	void movePositionRelative(glm::vec3 direction);
	glm::vec3 getPosition();

	void setFocus(glm::vec3 focus);
	void moveFocus(glm::vec3 direction);
	void aim(glm::vec3 direction);
	glm::vec3 getFocus();

	void move(glm::vec3 direction);
	void moveRelative(glm::vec3 direction);

	//Send all uniform matrices to the shader to be used.
	void send(glm::mat4 model);

	//Constructs a UniformMatrices instance to maintain matrices to be used in shaders
	//UniformMatrices(GLuint shaders, glm::mat4 perspective, glm::mat4 view, glm::mat4 model);
	UniformMatrices(GLuint shaders,
					glm::mat4 projection = glm::perspective(60.0f, 2.0f, 0.1f, 1000.0f),
					glm::vec3 position = glm::vec3(12,12,12), glm::vec3 focus = glm::vec3(0,0,0),
					glm::mat4 model = glm::mat4(1.0f), 
					const char * combinedName = "Combined", 
					const char * viewName = "View", 
					const char * modelName = "Model");
};

//Canvas

class Canvas : public wxGLCanvas {
	wxSize size;
	UniformMatrices* viewer;
	ColoredObject* object;
	wxDECLARE_EVENT_TABLE();
public:
	Canvas(wxWindow* parent, wxSize size);
	void GenerateGeometry();
	void Paint(wxPaintEvent& WXUNUSED(event));
	void Initialize(wxGLContext* context);
};
