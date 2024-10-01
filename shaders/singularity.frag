#version 450

uniform bool draw_zero_k;
uniform bool transparent;

flat in int k;

out vec4 outColor;

void main() {
	float a = pow(0.6, abs(k));
	float b = pow(0.9, abs(k));

	float alpha = transparent ? 0.35 : 1.0;

	if (k != 0 || (draw_zero_k && !transparent)) {
		if (k > 0) {
			outColor = vec4(1.0, b, a, alpha);
		} else {
			outColor = vec4(a, b, 1.0, alpha);
		}
	} else {
		discard;
	}
}