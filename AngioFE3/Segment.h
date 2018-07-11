#pragma once
#include <FECore/vec3d.h>
#include <FECore/FEDomain.h>
#include <FECore/FESolidDomain.h>
#include "Tip.h"


//-----------------------------------------------------------------------------
// Microvessels are represented by a collection of line segments. 
// Growth is represented by the addition of new segments onto the 
// active tips of exisiting segments. Within the simulation, these 
// line segments are found in the Segment class.
class Segment  
{
public:


    Segment(){}
	vec3d Direction(FEMesh * mesh) const;
	vec3d NatcAtTime(double time)const;
	int GetInitialFragmentID()const;
	double Length(FEMesh * mesh)const;
	double RefLength(FEMesh * mesh) const;//calculate length in reference frame
	Tip * front = nullptr;
	Tip * back = nullptr;
	Segment * parent = nullptr;
	bool is_branch_base = false;
	bool processed = false;//used by the brancher to determine if branches can grow from this segment
private:
	
};
