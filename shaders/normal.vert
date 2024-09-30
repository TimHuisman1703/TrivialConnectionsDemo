#version 450

uniform mat4 mvp;

in vec3 pos;
in vec3 normal;

out vec3 worldPos;
out vec3 worldNormal;

void main() {
    gl_Position = mvp * vec4(pos, 1.0);

    worldPos = pos;
    worldNormal = normal;
}