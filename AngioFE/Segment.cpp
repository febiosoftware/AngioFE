#include "StdAfx.h"
#include "Segment.h"
#include <FECore/vec3d.h>
#include "angio3d.h"
#include <iostream>

vec3d Segment::Direction(FEMesh * mesh) const
{
	vec3d p0 = front->GetPosition(mesh);
	vec3d p1 = back->GetPosition(mesh);
	vec3d r = p0 - p1;
	if (r.norm() < 0.001) {
		r = front->direction;
	}
	r.unit();
	return r;
}

vec3d Segment::NatcAtTime(double time)const
{
	assert(time <= front->time && time >= back->time);
	vec3d front_local_pos = front->GetLocalPosition();
	vec3d back_local_pos = back->GetLocalPosition();
	return mix(front_local_pos, back_local_pos,(time -back->time)/(front->time - back->time) );
}

int Segment::GetInitialFragmentID()const
{
	return front->initial_fragment_id;
}

double Segment::Length(FEMesh * mesh)const
{
	return (front->GetPosition(mesh) - back->GetPosition(mesh)).norm();
}

double Segment::LengthAtTime(FEMesh * mesh, double time)const
{
	double contrib = 0;
	if(back->time > time)
	{
		contrib = 0;
	}
	else if((back->time <= time) && (front->time <= time))
	{
		contrib = 1;
	}
	else
	{
		contrib = (time - back->time) / (front->time - back->time);
	}
	assert(contrib >= -0.001 && contrib <= 1.001);

	return (front->GetPosition(mesh) - back->GetPosition(mesh)).norm() * contrib;
}

double Segment::RefLength(FEMesh * mesh)const
{
	return (front->GetRefPosition(mesh) - back->GetRefPosition(mesh)).norm();
}
