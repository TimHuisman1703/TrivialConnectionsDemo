#version 450

struct Vertex {
	vec3 position;
	int k;
};

uniform mat4 mvp;
uniform int selected_vertex_idx;
layout(std140) uniform vertex_buffer
{
    int vertex_amount;
    Vertex vertices[1000];
};

in vec3 pos;

out int k;

#define RADIUS 0.02
#define SELECTED_RADIUS 0.035

void main() {
    Vertex v = vertices[gl_InstanceID];
    vec3 worldPos = (pos * (selected_vertex_idx == gl_InstanceID ? SELECTED_RADIUS : RADIUS)) + v.position;
    gl_Position = mvp * vec4(worldPos, 1);

    k = v.k;
}
