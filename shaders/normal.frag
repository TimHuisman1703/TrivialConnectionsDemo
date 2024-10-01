#version 450

in vec3 fragNormal;

out vec4 outColor;

void main() {
    outColor = vec4(0.6 + 0.15 * fragNormal, 1.0);
}