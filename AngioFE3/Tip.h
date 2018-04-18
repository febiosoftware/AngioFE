#pragma once
#include <FECore/vec3d.h>
#include "AngioElement.h"

class Segment;

class Tip
{
public:
	vec3d local_pos;
	AngioElement * angio_element;
	double time;
	vec3d growth_velocity;
	int face;//-1 for inside the element
	Segment * parent=nullptr;//may be nullptr for no parent segment
	int initial_fragment_id;
	bool is_branch = false;
	vec3d direction;//only used if the tip is a branch
	vec3d GetDirection(FEMesh* mesh) const;
	vec3d GetPosition(FEMesh * mesh) const;
};

