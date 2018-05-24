#pragma once
#include <FECore/FEMaterial.h>
#include <FEBioMech/FESPRProjection.h>
class Tip;
//the component of vessel growth that is dependent on position within the mesh
class PositionDependentDirection : public FEMaterial
{
public:
	explicit PositionDependentDirection(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~PositionDependentDirection() {}
	virtual vec3d ApplyModifiers(vec3d prev, Tip* tip, FEMesh* mesh)=0;
	virtual void Update(FEMesh * mesh){} //may be used to get values from loadcurves that modify the behavior as a whole
};

//The component of the mixture that represents the direction of the conttribution from the previous segment
class PreviousSegmentContribution : public FEMaterial
{
public:
	explicit PreviousSegmentContribution(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~PreviousSegmentContribution() {}
	virtual vec3d ApplyModifiers(vec3d prev, Tip* tip, FEMesh* mesh)=0;
	virtual void Update(FEMesh * mesh) {}
};

//determines the mixture between the position dependent direction and previous segment contribution 
class ContributionMix : public FEMaterial
{
public:
	explicit ContributionMix(FEModel* pfem) : FEMaterial(pfem){}
	virtual ~ContributionMix(){}
	virtual double ApplyModifiers(double prev, Tip* tip, FEMesh* mesh)=0;
	virtual void Update(FEMesh * mesh) {}
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
	vec3d ApplyModifiers(vec3d prev, Tip* tip, FEMesh* mesh) override;
	void Update(FEMesh * mesh) override;
private:
	double contribution = 1.0;
	FESPRProjection projection;
};

class ConcentrationPDD :public PositionDependentDirection
{
public:
	explicit ConcentrationPDD(FEModel* pfem) : PositionDependentDirection(pfem) {}
	virtual ~ConcentrationPDD() {}
	vec3d ApplyModifiers(vec3d prev, Tip* tip, FEMesh* mesh) override;
};

class DensityGradientPDD :public PositionDependentDirection
{
public:
	explicit DensityGradientPDD(FEModel* pfem) : PositionDependentDirection(pfem) {}
	virtual ~DensityGradientPDD() {}
	vec3d ApplyModifiers(vec3d prev, Tip* tip, FEMesh* mesh) override;
};

class PreviousSegmentPSC :public PreviousSegmentContribution
{
	//overrides the previous segment contribution with the previous segment direction
public:
	explicit PreviousSegmentPSC(FEModel* pfem) : PreviousSegmentContribution(pfem) {}
	virtual ~PreviousSegmentPSC() {}
	vec3d ApplyModifiers(vec3d prev, Tip* tip, FEMesh* mesh) override;
};