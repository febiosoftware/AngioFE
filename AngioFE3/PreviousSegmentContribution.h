#pragma once
#include <FECore/FEMaterial.h>
#include <FEBioMech/FESPRProjection.h>
#include "AngioElement.h"
#include "VariableInterpolation.h"
class FEAngio;

class Tip;

//The component of the mixture that represents the direction of the conttribution from the previous segment
class PreviousSegmentContribution : public FEMaterial
{
public:
	explicit PreviousSegmentContribution(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~PreviousSegmentContribution() {}
	virtual vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh) = 0;
	virtual void Update(FEMesh * mesh) {}
};

//The component of the mixture that represents the direction of the conttribution from the previous segment
class PreviousSegmentContributionManager : public FEMaterial
{
public:
	explicit PreviousSegmentContributionManager(FEModel* pfem) : FEMaterial(pfem) { AddProperty(&psc_modifiers, "psc_modifier"); psc_modifiers.m_brequired = false; }
	virtual ~PreviousSegmentContributionManager() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh);
	void Update(FEMesh * mesh) {}
private:
	FEVecPropertyT<PreviousSegmentContribution> psc_modifiers;
};

class PreviousSegmentPSC : public PreviousSegmentContribution
{
	//overrides the previous segment contribution with the previous segment direction
public:
	explicit PreviousSegmentPSC(FEModel* pfem) : PreviousSegmentContribution(pfem) {}
	virtual ~PreviousSegmentPSC() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh) override;
};