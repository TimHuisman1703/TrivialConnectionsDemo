#version 450

uniform vec4 albedo;

out vec4 outColor;

void main() {
    outColor = albedo;
}