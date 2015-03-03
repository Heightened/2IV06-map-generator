#version 420 core

in vec3 fPosition;
in vec3 fNormal;
in vec3 fColor;

in vec3 eyeVector;
in vec3 lightVector;
in vec3 normalVector;

out vec3 color;

void main() {
	vec3 light = vec3(50,50,50);

	color = fColor;

	// Light emission properties
	vec3 LightColor = vec3(1,1,1);
	float LightPower = 1000.0f;
	
	// Material properties
	vec3 MaterialDiffuseColor =	fColor;
	vec3 MaterialAmbientColor = 0.1f * MaterialDiffuseColor;
	vec3 MaterialSpecularColor = vec3(0.3,0.3,0.3);

	// Light angles
	float distance = length( light - fPosition );
	vec3 n = normalize(normalVector);
	vec3 l = normalize(lightVector);
	float cosTheta = clamp( dot( n,l ), 0,1 );
	vec3 E = normalize(eyeVector);
	vec3 R = reflect(-l,n);
	float cosAlpha = clamp( dot( E,R ), 0,1 );

	// Combine to calculate color
	color = /*vec3(0.2f,0.2f,0.2f) + */LightPower / (distance*distance) * (MaterialDiffuseColor * cosTheta + MaterialSpecularColor * cosAlpha);
}
