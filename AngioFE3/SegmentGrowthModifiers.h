#pragma once
#include <FECore/FEMaterial.h>
//the component of vessel growth that is dependent on position within the mesh
class PositionDependentDirection : public FEMaterial
{
public:
	explicit PositionDependentDirection(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~PositionDependentDirection() {}
	virtual vec3d ApplyModifiers(vec3d prev)=0;
};

//The component of the mixture that represents the direction of the conttribution from the previous segment
class PreviousSegmentContribution : public FEMaterial
{
public:
	explicit PreviousSegmentContribution(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~PreviousSegmentContribution() {}
	virtual vec3d ApplyModifiers(vec3d prev)=0;
};

//determines the mixture between the position dependent direction and previous segment contribution 
class ContributionMix : public FEMaterial
{
public:
	explicit ContributionMix(FEModel* pfem) : FEMaterial(pfem){}
	virtual ~ContributionMix(){}
	virtual double ApplyModifiers(double prev)=0;
};

//get the length over one unit of time at the given position
class SegmentGrowthLength :public FEMaterial
{
public:
	explicit SegmentGrowthLength(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~SegmentGrowthLength() {}
	virtual double ApplyModifiers(double prev) = 0;
};


class FiberPDD :public PositionDependentDirection 
{
public:
	explicit FiberPDD(FEModel* pfem) : PositionDependentDirection(pfem) {}
	virtual ~FiberPDD() {}
	vec3d ApplyModifiers(vec3d prev) override;
};

class ConcentrationPDD :public PositionDependentDirection
{
public:
	explicit ConcentrationPDD(FEModel* pfem) : PositionDependentDirection(pfem) {}
	virtual ~ConcentrationPDD() {}
	vec3d ApplyModifiers(vec3d prev) override;
};

class DensityGradientPDD :public PositionDependentDirection
{
public:
	explicit DensityGradientPDD(FEModel* pfem) : PositionDependentDirection(pfem) {}
	virtual ~DensityGradientPDD() {}
	vec3d ApplyModifiers(vec3d prev) override;
};

class PreviousSegmentPSC :public PreviousSegmentContribution
{
public:
	explicit PreviousSegmentPSC(FEModel* pfem) : PreviousSegmentContribution(pfem) {}
	virtual ~PreviousSegmentPSC() {}
	vec3d ApplyModifiers(vec3d prev) override;
};