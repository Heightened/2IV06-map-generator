#include "Canvas.h"

GLuint initShaders(const char * VertexShaderFile, const char * FragmentShaderFile) {
	const GLubyte * a = glGetString(GL_SHADING_LANGUAGE_VERSION_ARB);
	std::cout << a << std::endl;

	//Line string to buffer code of shaders
	std::string Line = "";

	// Read the Vertex Shader code from the file
	std::string VShaderSource;
	std::ifstream VShaderStream(VertexShaderFile, std::ios::in);
	while(getline(VShaderStream, Line)) {
		VShaderSource += "\n" + Line;
	}
	VShaderStream.close();

	// Read the Fragment Shader code from the file
	std::string FShaderSource;
	std::ifstream FShaderStream(FragmentShaderFile, std::ios::in);
	while(getline(FShaderStream, Line)) {
		FShaderSource += "\n" + Line;
	}
	FShaderStream.close();

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
		printf("%s\n", &VShaderErrorMessage[0]);
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
		printf("%s\n", &FragmentShaderErrorMessage[0]);
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
		printf("%s\n", &ProgramErrorMessage[0]);
	}

	return Program;
}

GLuint initShaders() {
	initShaders("shaders/PreviewVertexShader.glsl", "shaders/PreviewFragmentShader.glsl");
}

bool initCanvas() {
	glewExperimental = true;
	bool glew = glewInit();
	if (glew != GLEW_OK) {
		exit(EXIT_FAILURE);
	}

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glEnable(GL_CULL_FACE);
	glClearColor(0,0,0,1.0f);
}
