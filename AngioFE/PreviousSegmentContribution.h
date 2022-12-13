#pragma once
#include <FECore/FEMaterial.h>
#include <FECore/FESPRProjection.h>
#include "AngioElement.h"
#include "VariableInterpolation.h"
class FEAngio;

class Tip;

//! The component of the mixture that represents the direction of the conttribution from the previous segment
class PreviousSegmentContribution : public FEMaterialProperty
{
	FECORE_BASE_CLASS(PreviousSegmentContribution)
public:
	//! constructor
	explicit PreviousSegmentContribution(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~PreviousSegmentContribution() {}
	//! returns the direction contribution from the previous segment
	virtual vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh) = 0;
	//! may update this based on the timestep
	virtual void Update(FEMesh * mesh) {}
};

//! The component of the mixture that represents the direction of the conttribution from the previous segment
class PreviousSegmentContributionManager : public FEMaterialProperty
{
	FECORE_BASE_CLASS(PreviousSegmentContributionManager)
public:
	//! constructor
	explicit PreviousSegmentContributionManager(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~PreviousSegmentContributionManager() {}
	//! returns the combination contribution of all PSC modifiers
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh);
	//! may update this based on the timestep
	void Update(FEMesh * mesh) {}
protected:
	DECLARE_FECORE_CLASS()
private:
	std::vector<PreviousSegmentContribution*>	psc_modifiers;	//!< pointers to elastic materials
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


//! The component of the mixture that represents the direction of the conttribution from the previous segment
class ProtoPreviousSegmentContribution : public FEMaterialProperty
{
	FECORE_BASE_CLASS(ProtoPreviousSegmentContribution)
public:
	//! constructor
	explicit ProtoPreviousSegmentContribution(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~ProtoPreviousSegmentContribution() {}
	//! returns the direction contribution from the previous segment
	virtual vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh) = 0;
	//! may update this based on the timestep
	virtual void Update(FEMesh * mesh) {}
};

//! The component of the mixture that represents the direction of the conttribution from the previous segment
class ProtoPreviousSegmentContributionManager : public FEMaterialProperty
{
	FECORE_BASE_CLASS(ProtoPreviousSegmentContributionManager)
public:
	//! constructor
	explicit ProtoPreviousSegmentContributionManager(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~ProtoPreviousSegmentContributionManager() {}
	//! returns the combination contribution of all PSC modifiers
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh);
	//! may update this based on the timestep
	void Update(FEMesh * mesh) {}
protected:
	DECLARE_FECORE_CLASS()
private:
	std::vector<ProtoPreviousSegmentContribution*>	proto_psc_modifiers;	//!< pointers to elastic materials
};

//! a contribution that is entirely the previous segment direction
class ProtoPreviousSegmentPSC : public ProtoPreviousSegmentContribution
{
	//overrides the previous segment contribution with the previous segment direction
public:
	//! constructor
	explicit ProtoPreviousSegmentPSC(FEModel* pfem) : ProtoPreviousSegmentContribution(pfem) {}
	virtual ~ProtoPreviousSegmentPSC() {}
	//returns the direction of the previous segment
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh) override;
};