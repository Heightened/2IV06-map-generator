#pragma once

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#define PREDEFINEDOBJECTS 

//Buffer operations
GLuint generateBuffer(GLenum target, size_t size, const GLvoid * vectors, GLenum usage);
void loadBuffer(GLuint buffer, GLuint attribute, size_t vectorsize, GLenum type, bool normalize);

//Object classes
class Object {
	const GLfloat * vertices;
	const GLfloat * normals;

protected:
	const size_t vertexCount;
	GLuint verticesId;
	GLuint normalsId;

public:
	virtual void draw() {};
	
	Object(size_t vertexCount, const GLfloat ** attributes);
	~Object();
};

class ColoredObject : public Object {
	const GLfloat * colors;

protected:
	GLuint colorsId;

public:
	virtual void draw();

	ColoredObject(size_t vertexCount, const GLfloat ** attributes);
	~ColoredObject();
};

class TexturedObject : public Object {
	const GLfloat * uvcoords;

protected:
	GLuint uvcoordsId;

public:
	virtual void draw();

	TexturedObject(size_t vertexCount, const GLfloat ** attributes);
	~TexturedObject();
};

class Attribute {
protected:
	size_t length;
	GLfloat * values;

	//Constructor that allows an empty Attribute object to be constructed by child classes
	Attribute(size_t l);
public:
	size_t size();
	operator GLfloat * ();
	//Translations are only applicable on attributes consiting of vec3s
	void translate(glm::vec3 v);
	void scale(float scalar);
	//Rotations are only applicable on attributes consiting of vec3s
	void rotateRad(glm::vec3 axis, float angle);
	void rotate(glm::vec3 axis, float angle);

	Attribute(size_t l, GLfloat * v);

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
	Normals(size_t vertexCount, GLfloat * vertices, bool smooth);
};

class SolidColor : public Attribute {
public:
	SolidColor(size_t vertexCount, float r, float g, float b);
};
