#include <glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <string>
#include <fstream>
#include <vector>
#include <iostream>


GLuint initShaders(const char * VertexShaderFile, const char * FragmentShaderFile);

GLuint initShaders();

bool initCanvas();
