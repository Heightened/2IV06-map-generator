#version 420 core

layout(location = 0) in vec3 vertex;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec3 color;

out vec3 fPosition;
out vec3 fNormal;
out vec3 fColor;

out vec3 eyeVector;
out vec3 lightVector;
out vec3 normalVector;

//View uniforms
uniform mat4 Combined;
uniform mat4 Model;
uniform mat4 View;

//Lighting uniforms
//uniform vec3 light;

void main() {
	vec3 light = vec3(50,50,50);

	//Position of the vertex on the screen
	gl_Position = Combined * vec4(vertex,1);
	
	//Position of the vertex in the world
	fPosition = (Model * vec4(vertex,1)).xyz;

	//Pass on unused inputs
	fNormal = normal;
	fColor = color;
	
	//Direction of reflecting light from vertex to the camera eye
	eyeVector = vec3(0,0,0) - (View * Model * vec4(vertex,1)).xyz;

	//Vector from the vertex to the lightsource
	lightVector = (View * vec4(light,1)).xyz + eyeVector;
	
	//Normal of the the vertex in camera space
	normalVector = (View * Model * vec4(normal,0)).xyz;
}
