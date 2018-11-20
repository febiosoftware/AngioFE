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

	//! Constructor
    Segment(){}
	//! Returns the direction of the segment
	vec3d Direction(FEMesh * mesh) const;
	//! Return the natural coordinates that the segment would be at at a given time
	vec3d NatcAtTime(double time)const;
	//! returns the initial fragment id
	int GetInitialFragmentID()const;
	//! returns the length of the segment
	double Length(FEMesh * mesh)const;
	//! returns the length of the segment at a given time
	double LengthAtTime(FEMesh * mesh, double time)const;
	//! returns the length of a segment in the reference configuration of the mesh
	double RefLength(FEMesh * mesh) const;//calculate length in reference frame
	//! front tip
	Tip * front = nullptr;
	//! back tip
	Tip * back = nullptr;
	//! the segment that this segment is grown from
	Segment * parent = nullptr;
	//! is the first segment in the branch
	bool is_branch_base = false;
	//! Used by the brancher to determine if branches can grow from this segment
	bool processed = false;
private:
	
};
