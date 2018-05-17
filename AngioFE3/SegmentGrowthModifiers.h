#pragma once
#include <FECore/FEMaterial.h>
//the component of vessel growth that is dependent on position within the mesh
class PositionDependentDirection : public FEMaterial
{
	virtual vec3d ApplyModifiers(vec3d prev)=0;
};

//The component of the mixture that represents the direction of the conttribution from the previous segment
class PreviousSegmentContribution : public FEMaterial
{
	virtual vec3d ApplyModifiers(vec3d prev)=0;
};

//determines the mixture between the position dependent direction and previous segment contribution 
class ContributionMix : public FEMaterial
{
	virtual double ApplyModifiers(double prev)=0;
};
