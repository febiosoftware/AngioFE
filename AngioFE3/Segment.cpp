#include "StdAfx.h"
#include "Segment.h"
#include <FECore/vec3d.h>

vec3d Segment::Direction(FEMesh * mesh) const
{
	vec3d p0 = front->GetPosition(mesh);
	vec3d p1 = back->GetPosition(mesh);
	vec3d r = p0 - p1;
	r.unit();
	return r;
}
