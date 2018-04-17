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

private:
	Tip * front, * back;
	Segment * parent;
	bool is_branch_base = false;
};
