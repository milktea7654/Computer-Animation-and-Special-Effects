#pragma once

#include <string>

static std::string point_vs_src = R"(
#version 330
in vec2 attrPosition;
in vec3 attrColor;

uniform vec2 domainSize;
uniform float pointSize;
uniform float drawDisk;

out vec3 fragColor;
out float fragDrawDisk;

void main() 
{
    gl_Position =
	    vec4(attrPosition * vec2(2.0 / domainSize.x, 2.0 / domainSize.y) - vec2(1.0, 1.0), 0.0, 1.0);

    gl_PointSize = pointSize;
    fragColor = attrColor;
    fragDrawDisk = drawDisk;
}
)";

static std::string point_fs_src = R"(
#version 330
precision mediump float;
in vec3 fragColor;

uniform float drawDisk;
uniform float pointSize;

out vec4 color;
void main() 
{
	if (drawDisk == 1.0) 
    {
        vec2 coord = gl_PointCoord - vec2(0.5);
        float distSquared = dot(coord, coord);
        
        if (distSquared > 0.25) 
            discard;
	}
	color = vec4(fragColor, 1.0);
}
)";

static std::string mono_color_vs_src = R"(
#version 330

uniform vec2 attrPosition;
uniform vec2 domainSize;
uniform float pointSize;

void main() 
{
    gl_Position =
	    vec4(attrPosition * vec2(2.0 / domainSize.x, 2.0 / domainSize.y) - vec2(1.0, 1.0), 0.0, 1.0);
    gl_PointSize = pointSize;
}
)";

static std::string mono_color_fs_src = R"(
#version 330
precision mediump float;

out vec4 color;
void main() 
{
    vec2 coord = gl_PointCoord - vec2(0.5);
    float distSquared = dot(coord, coord);
        
    if (distSquared > 0.25) 
        discard;
	color = vec4(1.0, 1.0, 0.0, 1.0);
}
)";