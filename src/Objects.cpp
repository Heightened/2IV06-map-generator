#define PI 3.1415926535897932384626433832795

#include "Objects.h"

#include <cmath>
#include <cstdio>

#include <wx/log.h> 

//Methods to set up buffers to contain for example vertex locations or similar attributes

GLuint generateBuffer(GLenum target, size_t size, const GLvoid * vectors, GLenum usage) {
	GLuint id = -1;
	glGenBuffers(1, &id);
	glBindBuffer(target, id);
	glBufferData(target, size, vectors, usage);
	return id;
}

void loadBuffer(GLuint buffer, GLuint attribute, size_t vectorsize, GLenum type, bool normalize) {
	glEnableVertexAttribArray(attribute);
	glBindBuffer(GL_ARRAY_BUFFER, buffer);
	glVertexAttribPointer(attribute, vectorsize, type, normalize, 0, 0); 
}

//Objects using buffers

Object::Object(size_t _vertexCount, const GLfloat ** attributes) : vertexCount(_vertexCount), vertices(attributes[0]), normals(attributes[1]) {
	verticesId = generateBuffer(GL_ARRAY_BUFFER, vertexCount*3*4, vertices, GL_STATIC_DRAW);
	normalsId = generateBuffer(GL_ARRAY_BUFFER, vertexCount*3*4, normals, GL_STATIC_DRAW);
}

Object::~Object() {
	delete[] vertices;
	delete[] normals;

	const GLuint buffers[] = {verticesId, normalsId};
	glDeleteBuffers(2, buffers);
}

void ColoredObject::draw() {
	loadBuffer(verticesId, 0, 3, GL_FLOAT, false);
	loadBuffer(normalsId, 1, 3, GL_FLOAT, false);
	loadBuffer(colorsId, 2, 3, GL_FLOAT, false);

	glDrawArrays(GL_TRIANGLES, 0, vertexCount); // 3 indices per vertex

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(2);
}

ColoredObject::ColoredObject(size_t _vertexCount, const GLfloat ** attributes) : Object(_vertexCount, attributes), colors(attributes[2]) {
	colorsId = generateBuffer(GL_ARRAY_BUFFER, vertexCount*3*4, colors, GL_STATIC_DRAW);
	printf("\nvertexCount: %i\n", vertexCount);
}

ColoredObject::~ColoredObject() {
	delete[] colors;

	const GLuint buffers[] = {colorsId};
	glDeleteBuffers(1, buffers);
}

void TexturedObject::draw() {
	loadBuffer(verticesId, 0, 3, GL_FLOAT, false);
	loadBuffer(normalsId, 1, 3, GL_FLOAT, false);
	loadBuffer(uvcoordsId, 3, 2, GL_FLOAT, false);

	glDrawArrays(GL_TRIANGLES, 0, vertexCount); // 3 indices per vertex

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(3);
}

TexturedObject::TexturedObject(size_t _vertexCount, const GLfloat ** attributes) : Object(_vertexCount, attributes), uvcoords(attributes[2]) {
	uvcoordsId = generateBuffer(GL_ARRAY_BUFFER, vertexCount*3*4, uvcoords, GL_STATIC_DRAW);
}

TexturedObject::~TexturedObject() {
	delete[] uvcoords;

	const GLuint buffers[] = {uvcoordsId};
	glDeleteBuffers(1, buffers);
}

//Attribute class

Attribute::Attribute(size_t l) : length(l), values(new GLfloat[l]) {
}

size_t Attribute::size() {
	return length;
}

Attribute::operator GLfloat * () { 
	return values; 
}

void Attribute::translate(glm::vec3 v) {
	size_t n = size();
	for (size_t i = 0; i < n; i += 3) {
		values[i] += v[0];
		values[i+1] += v[1];
		values[i+2] += v[2];
	}
}

void Attribute::scale(float scalar) {
	size_t n = size();
	for (size_t i = 0; i < n; i++) {
		values[i] *= scalar;
	}
}

void Attribute::rotateRad(glm::vec3 axis, float angle) {
	//convert to degrees
	angle = (float)(angle / PI) * 180;

	rotate(axis, angle);
}

void Attribute::rotate(glm::vec3 axis, float angle) {
	size_t n = size();
	
	for (size_t i = 0; i < n; i += 3) {
		glm::vec3 value(values[i], values[i+1], values[i+2]);
		glm::mat4x4 rotation(1.0f);
		rotation = glm::rotate(rotation, angle, axis);
		value = glm::vec3(rotation * glm::vec4(value, 1.0f));
		values[i] = value.x;
		values[i+1] = value.y;
		values[i+2] = value.z;
	}
}

Attribute::Attribute(size_t l, GLfloat * v) : length(l), values(new GLfloat[l]) {
	for (int i = 0; i < l; i++) {
		values[i] = v[i];
	}
}

Attribute operator+ (const Attribute &first, const Attribute &second) {
	size_t n = first.length + second.length;
	GLfloat * result = new GLfloat[n];
	for (size_t i = 0; i < first.length; i++) {
		result[i] = first.values[i];
	}
	for (size_t i = 0; i < second.length; i++) {
		result[i+first.length] = second.values[i];
	}
	//printf("\n%i\n", n);
	return Attribute(n, result);
}

//End Attribute class

//Vertex position Attributes

QuadVertices::QuadVertices(glm::vec3 p, glm::vec3 x, glm::vec3 y) : Attribute(6*3) {
	glm::vec3 v = x + y;
	//This GLfloat only exists within the scope of this constructor and therefor a pointer to it will not have the desired result.
	GLfloat quadVertices[] = {
		p.x+0.0f, p.y+0.0f, p.z+0.0f,
		p.x+v.x, p.y+v.y, p.z+v.z,
		p.x+y.x, p.y+y.y, p.z+y.z,

		p.x+0.0f, p.y+0.0f, p.z+0.0f,
		p.x+x.x, p.y+x.y, p.z+x.z,
		p.x+v.x, p.y+v.y, p.z+v.z
	};
	values = new GLfloat[length];
	for (int i = 0; i < length; i++) {
		values[i] = quadVertices[i];
	}
}

BoxVertices::BoxVertices(float x, float y, float z) : Attribute(6*6*3) {
	glm::vec3 xv(x,0,0), yv(0,y,0), zv(0,0,z),
			o(-x/2,-y/2,0), a(o+xv), b(a+yv), c(o+yv), d(o+zv);

	Attribute boxVertices = //Upstanding faces
							QuadVertices(o, xv, zv) +
							QuadVertices(a, yv, zv) +
							QuadVertices(b, -xv, zv) +
							QuadVertices(c, -yv, zv) +
							//Top and bottom faces
							QuadVertices(o, yv, xv) +
							QuadVertices(d, xv, yv);
	
	values = new GLfloat[length];
	for (int i = 0; i < length; i++) {
		values[i] = boxVertices[i];
	}
}

SphereVertices::SphereVertices(float r, int subdivisions) : Attribute(20 * 3 * 3 * (int)glm::pow(4.0f, (float)subdivisions)) {
	//printf("size: %i", size());

	float z0 = -r, z1 = -r/2, z2 = r/2, z3 = r;
	
	glm::vec2 dg[10]; //Decagon, containing the vertices of two pentagons
	float theta = 0;
	for (int i = 0; i < 10; i++, theta += PI/5.0f) {
		dg[i] = (r * glm::sin(glm::acos(z2/r))) * glm::vec2(glm::cos(theta), glm::sin(theta));
	}

	Attribute vertices(0, 0);

	//Top and Bottom interlock in the middle
	for (int i = 0; i < 10; i++) {
		if (i % 2 == 0) {
			GLfloat triangles[] = {
				//Middle section triangle
				dg[i].x, dg[i].y, z1,
				dg[(i+2)%10].x, dg[(i+2)%10].y, z1,
				dg[(i+1)%10].x, dg[(i+1)%10].y, z2,
				//Bottom triangle
				dg[(i+2)%10].x, dg[(i+2)%10].y, z1,
				dg[i].x, dg[i].y, z1,
				0, 0, z0
			};
			vertices = vertices + Attribute(6 * 3, &triangles[0]);
		} else {
			GLfloat triangles[] = {
				//Middle section triangle
				dg[(i+2)%10].x, dg[(i+2)%10].y, z2,
				dg[i].x, dg[i].y, z2,
				dg[(i+1)%10].x, dg[(i+1)%10].y, z1,
				//Top triangle
				dg[i].x, dg[i].y, z2,
				dg[(i+2)%10].x, dg[(i+2)%10].y, z2,
				0, 0, z3
			};
			vertices = vertices + Attribute(6 * 3, &triangles[0]);
		}
	}

	if (subdivisions > 0) {
		//Subdivide
		for (int n = subdivisions; n > 0; n--) {

			Attribute newvertices(0, 0);
			for (size_t i = 0; i < length; i += 9) {
				//Take triangle ABC with A = Vi -> Vi+2, B = Vi+3 -> Vi+5 and C = Vi+6 -> Vi+8
				glm::vec3 
					A(vertices[i], vertices[i+1], vertices[i+2]), 
					B(vertices[i+3], vertices[i+4], vertices[i+5]), 
					C(vertices[i+6], vertices[i+7], vertices[i+8]),

					AB((A.x+B.x)/2, (A.y+B.y)/2, (A.z+B.z)/2),
					BC((B.x+C.x)/2, (B.y+C.y)/2, (B.z+C.z)/2),
					CA((C.x+A.x)/2, (C.y+A.y)/2, (C.z+A.z)/2);

				AB = r * AB / glm::length(AB);
				BC = r * BC / glm::length(BC);
				CA = r * CA / glm::length(CA);

				GLfloat triangles[] = {
					A.x, A.y, A.z,
					AB.x, AB.y, AB.z,
					CA.x, CA.y, CA.z,

					B.x, B.y, B.z,
					BC.x, BC.y, BC.z,
					AB.x, AB.y, AB.z,

					C.x, C.y, C.z,
					CA.x, CA.y, CA.z,
					BC.x, BC.y, BC.z,

					AB.x, AB.y, AB.z,
					BC.x, BC.y, BC.z,
					CA.x, CA.y, CA.z
				};
				newvertices = newvertices + Attribute(12 * 3, &triangles[0]);
			}
			vertices = newvertices;
		}
	}

	//Set values
	values = new GLfloat[length];
	for (size_t i = 0; i < length; i++) {
		values[i] = vertices[i];
	}
}

//Vertex Attributes

Normals::Normals(size_t vertexCount, GLfloat * vertices, bool smooth) : Attribute(vertexCount*3) {
	//Process vertices to find normals for each triangle
	for (size_t i = 0; i < length; i += 9) {
		glm::vec3 a(vertices[i], vertices[i+1], vertices[i+2]);
		glm::vec3 b(vertices[i+3], vertices[i+4], vertices[i+5]);
		glm::vec3 c(vertices[i+6], vertices[i+7], vertices[i+8]);
		
		glm::vec3 normal = glm::cross((b - a), (c - a));
		//printf("%i < %i: %f, %f, %f\n",i, length, normal.x, normal.y, normal.z);

		for (int j = 0; j < 3; j++) {
			values[i+j*3] = normal.x;
			values[i+j*3+1] = normal.y;
			values[i+j*3+2] = normal.z;
		}
	}

	//Normalize
	for (size_t i = 0; i < length; i += 3) {
		glm::vec3 normal = glm::vec3(values[i], values[i+1], values[i+2]);
		float magnitude = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);

		values[i] = normal.x / magnitude;
		values[i+1] = normal.y / magnitude;
		values[i+2] = normal.z / magnitude;
	}

	//Smooth by correcting for shared vertices if specified
	if (smooth) {
		//int * commons = new int[length/3];
		GLfloat * sum = new GLfloat[length];

		for (size_t i = 0; i < length; i += 3) {
			//commons[i/3] = 0;
			sum[i] = 0;
			sum[i+1] = 0;
			sum[i+2] = 0;
			for (size_t j = 0; j < length; j += 3) {
				if (vertices[i] == vertices[j] && vertices[i+1] == vertices[j+1] && vertices[i+2] == vertices[j+2]) {
					//two vertices are equal
					//The normal should only be taken in to account if an identical one hasn't allready been added
					//For now we'll approximate by normalizing the vector components individually
					//commons[i/3]++;
					sum[i] += values[j];
					sum[i+1] += values[j+1];
					sum[i+2] += values[j+2];
				}
			}
		}

		values = sum;
		/*for (int i = 0; i < length; i += 3) {
			//int commoncount = commons[i/3];
			values[i] = (sum[i] != 0) ? sum[i] / glm::abs(sum[i]) : 0;// / commoncount;
			values[i+1] = (sum[i+1] != 0) ? sum[i+1] / glm::abs(sum[i+1]) : 0;// / commoncount;
			values[i+2] = (sum[i+2] != 0) ? sum[i+2] / glm::abs(sum[i+2]) : 0;// / commoncount;
		}*/
	}

	//Normalize
	for (size_t i = 0; i < length; i += 3) {
		glm::vec3 normal = glm::vec3(values[i], values[i+1], values[i+2]);
		float magnitude = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
		values[i] = normal.x / magnitude;
		values[i+1] = normal.y / magnitude;
		values[i+2] = normal.z / magnitude;
	}
}

SolidColor::SolidColor(size_t vertexCount, float r, float g, float b) : Attribute(vertexCount*3) {
	//values = new GLfloat[length];
	for (size_t i = 0; i < length; i += 3) {
		//printf("%i: %f, %f, %f\n",i/3, r, g, b);
		values[i] = r;
		values[i+1] = g;
		values[i+2] = b;
	}
}
