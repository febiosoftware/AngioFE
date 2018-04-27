#pragma once
#include <FECore/vec3d.h>
#include "AngioElement.h"

class Segment;

class Tip
{
public:
	Tip(){}
	Tip(Tip * other, FEMesh * mesh);//used to copy another tip
	vec3d local_pos;
	AngioElement * angio_element= nullptr;
	double time=0.0;
	vec3d growth_velocity;
	//int face;//-1 for inside the element
	AngioElement* face= nullptr;//the element where growth will take place in on the next growth step
	Segment * parent=nullptr;//may be nullptr for no parent segment
	int initial_fragment_id = -1;
	bool is_branch = false;
	bool use_direction = false;
	vec3d direction;//only used if the tip is a branch
	vec3d GetDirection(FEMesh* mesh) const;
	vec3d GetPosition(FEMesh * mesh) const;
};

