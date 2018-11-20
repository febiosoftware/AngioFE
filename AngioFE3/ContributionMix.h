#pragma once
#include <FECore/FEMaterial.h>
#include <FEBioMech/FESPRProjection.h>
#include "AngioElement.h"
#include "VariableInterpolation.h"
class FEAngio;

class Tip;

//determines the mixture between the position dependent direction and previous segment contribution 
class ContributionMix : public FEMaterial
{
public:
	explicit ContributionMix(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~ContributionMix() {}
	virtual double ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh) = 0;
	virtual void Update(FEMesh * mesh) {}
};

class ContributionMixManager : public FEMaterial
{
public:
	explicit ContributionMixManager(FEModel* pfem) : FEMaterial(pfem) { AddProperty(&cm_modifiers, "psc_modifier"); cm_modifiers.m_brequired = false; }
	virtual ~ContributionMixManager() {}
	double ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh);
	void Update(FEMesh * mesh) {}
private:
	FEVecPropertyT<ContributionMix> cm_modifiers;
};

class PSCPDDContributionMix : public ContributionMix
{
public:
	explicit PSCPDDContributionMix(FEModel* pfem) : ContributionMix(pfem) {}
	virtual ~PSCPDDContributionMix() {}
	double ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh) override;
	void Update(FEMesh * mesh) override;
protected:
	//! parameter list
	DECLARE_PARAMETER_LIST();
private:
	double psc_weight = 0.5;
};