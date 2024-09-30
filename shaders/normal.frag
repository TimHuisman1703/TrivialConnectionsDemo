#version 450

in vec3 worldPos;
in vec3 worldNormal;

out vec4 outColor;

void main() {
    outColor = vec4(0.6 + 0.15 * worldNormal, 1.0);
}