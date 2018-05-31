#pragma once
#include <FECore/vec3d.h>
#include "AngioElement.h"

class Segment;
class FEAngio;

class Tip
{
public:
	Tip(){}
	Tip(Tip * other, FEMesh * mesh);//used to copy another tip
	//sets the local positions clamps the values to -1 to 1
	void SetLocalPosition(vec3d pos);
	vec3d GetLocalPosition() const;
	AngioElement * angio_element= nullptr;//the element that contains the local_pos coordinates
	double time=0.0;
	double growth_velocity;
	//int face;//-1 for inside the element
	AngioElement* face= nullptr;//the element where growth originated from
	Segment * parent=nullptr;//may be nullptr for no parent segment
	int initial_fragment_id = -1;
	bool is_branch = false;
	bool use_direction = false;
	vec3d direction;//only used if the tip is a branch
	vec3d GetDirection(FEMesh* mesh) const;
	vec3d GetPosition(FEMesh * mesh) const;
	void PrintTipInfo(FEMesh *mesh) const;
	void PrintTipInfo(FEMesh *mesh, std::string title) const;
private:
	vec3d local_pos;
};

