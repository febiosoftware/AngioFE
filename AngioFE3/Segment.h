#pragma once
#include <FECore/vec3d.h>
#include <FECore/FEDomain.h>
#include <FECore/FESolidDomain.h>
#include "Tip.h"


//! two tips that form a line segment(in the natural coordinate element space), vascular networks are composed of these
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
	//! Length since last branch
	double LengthSinceBranch = 0;
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
	double time_emerged = 0;
private:
	
};
