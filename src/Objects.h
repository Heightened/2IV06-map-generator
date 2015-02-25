#include <glew.h>
#include <glm/glm.hpp>

//Buffer operations
GLuint generateBuffer(GLenum target, int size, const GLvoid * vectors, GLenum usage);
void loadBuffer(GLuint buffer, GLuint attribute, int vectorsize, GLenum type, bool normalize);

//Object classes
class Object {
	const GLfloat * vertices;
	const GLfloat * normals;

protected:
	const int vertexCount;
	GLuint verticesId;
	GLuint normalsId;

public:
	virtual void draw() {};
	
	Object(int vertexCount, const GLfloat ** attributes);
};

class ColoredObject : public Object {
	const GLfloat * colors;

protected:
	GLuint colorsId;

public:
	virtual void draw();

	ColoredObject(int vertexCount, const GLfloat ** attributes);
};

class TexturedObject : public Object {
	const GLfloat * uvcoords;

protected:
	GLuint uvcoordsId;

public:
	virtual void draw();

	TexturedObject(int vertexCount, const GLfloat ** attributes);
};

class Attribute {
protected:
	int length;
	GLfloat * values;

	//Constructor that allows an empty Attribute object to be constructed by child classes
	Attribute(int l);
public:
	int size();
	operator GLfloat * ();
	//Translations are only applicable on attributes consiting of vec3s
	void translate(glm::vec3 v);

	Attribute(int l, GLfloat * v);

	friend Attribute operator+ (const Attribute &first, const Attribute &second);
};

class QuadVertices : public Attribute {
public:
	QuadVertices(glm::vec3 p, glm::vec3 x, glm::vec3 y);
};

class BoxVertices : public Attribute {
public:
	BoxVertices(float x, float y, float z);
};

class SphereVertices : public Attribute {
public:
	SphereVertices(float r, int subdivisions);
};

class Normals : public Attribute {
public:
	Normals(int vertexCount, GLfloat * vertices, bool smooth);
};

class SolidColor : public Attribute {
public:
	SolidColor(int vertexCount, float r, float g, float b);
};
