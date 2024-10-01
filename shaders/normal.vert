#version 450

uniform mat4 mvp;

in vec3 pos;
in vec3 normal;

out vec3 fragNormal;

void main() {
    gl_Position = mvp * vec4(pos - normal * 0.002, 1.0);

    fragNormal = normal;
}