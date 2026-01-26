#pragma once
#include <string>

// Maintain shaders as strings

std::string frag = R"(
  #version 400

  in Fragment {
    vec4 color;
    vec2 mapping;
  }
  fragment;

  layout(location = 0) out vec4 fragmentColor;

  void main() {
    float r = dot(fragment.mapping, fragment.mapping);
    if (r > 1) discard;
    fragmentColor = vec4(fragment.color.rgb, 1 - r * r); // smooth alpha falloff
  }
)";

std::string vert = R"(
  #version 400

  layout(location = 0) in vec3 vertexPosition;
  layout(location = 1) in vec4 vertexColor;
  layout(location = 2) in vec2 vertexSize;
  // vertexSize is 2D texture cordinate, but we only use the x

  uniform mat4 al_ModelViewMatrix;
  uniform mat4 al_ProjectionMatrix;

  out Vertex {
    vec4 color;
    float size;
  } vertex;

  void main() {
    gl_Position = al_ModelViewMatrix * vec4(vertexPosition, 1.0);
    vertex.color = vertexColor;
    vertex.size = vertexSize.x;
  }
)";

std::string geo = R"(
  #version 400

  // take in a point and output a triangle strip with 4 vertices (aka a "quad")
  //
  layout(points) in;
  layout(triangle_strip, max_vertices = 4) out;
  // XXX ~ what about winding? should use triangles?

  uniform mat4 al_ProjectionMatrix;
  uniform float pointSize;

  in Vertex {
    vec4 color;
    float size;
  }
  vertex[];

  out Fragment {
    vec4 color;
    vec2 mapping;
  }
  fragment;

  void main() {
    mat4 m = al_ProjectionMatrix;   // rename to make lines shorter
    vec4 v = gl_in[0].gl_Position;  // al_ModelViewMatrix * gl_Position

    float r = pointSize;
    r *= vertex[0].size;

    gl_Position = m * (v + vec4(-r, -r, 0.0, 0.0));
    fragment.color = vertex[0].color;
    fragment.mapping = vec2(-1.0, -1.0);
    EmitVertex();

    gl_Position = m * (v + vec4(r, -r, 0.0, 0.0));
    fragment.color = vertex[0].color;
    fragment.mapping = vec2(1.0, -1.0);
    EmitVertex();

    gl_Position = m * (v + vec4(-r, r, 0.0, 0.0));
    fragment.color = vertex[0].color;
    fragment.mapping = vec2(-1.0, 1.0);
    EmitVertex();

    gl_Position = m * (v + vec4(r, r, 0.0, 0.0));
    fragment.color = vertex[0].color;
    fragment.mapping = vec2(1.0, 1.0);
    EmitVertex();

    EndPrimitive();
  }
)";
