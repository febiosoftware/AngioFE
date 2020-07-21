#pragma once
#include <FECore/FEMaterial.h>
#include <FECore/FESPRProjection.h>
#include "AngioElement.h"
#include "VariableInterpolation.h"
class FEAngio;

class Tip;

//! The component of the mixture that represents the direction of the conttribution from the previous segment
class PreviousSegmentContribution : public FEMaterial
{
public:
	//! constructor
	explicit PreviousSegmentContribution(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~PreviousSegmentContribution() {}
	//! returns the direction contribution from the previous segment
	virtual vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh) = 0;
	//! may update this based on the timestep
	virtual void Update(FEMesh * mesh) {}
};

//! The component of the mixture that represents the direction of the conttribution from the previous segment
class PreviousSegmentContributionManager : public FEMaterial
{
public:
	//! constructor
	explicit PreviousSegmentContributionManager(FEModel* pfem) : FEMaterial(pfem) { AddClassProperty(this, &psc_modifiers, "psc_modifiers", FEProperty::Optional); }
	virtual ~PreviousSegmentContributionManager() {}
	//! returns the combination contribution of all PSC modifiers
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh);
	//! may update this based on the timestep
	void Update(FEMesh * mesh) {}
private:
private:
	std::vector<PreviousSegmentContribution*>	psc_modifiers;	//!< pointers to elastic materials
	//PreviousSegmentContribution* psc_modifiers;
};

//! a contribution that is entirely the previous segment direction
class PreviousSegmentPSC : public PreviousSegmentContribution
{
	//overrides the previous segment contribution with the previous segment direction
public:
	//! constructor
	explicit PreviousSegmentPSC(FEModel* pfem) : PreviousSegmentContribution(pfem) {}
	virtual ~PreviousSegmentPSC() {}
	//returns the direction of the previous segment
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh) override;
};