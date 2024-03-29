#pragma once
#include <FECore/FEMaterial.h>
#include <FECore/FESPRProjection.h>
#include "AngioElement.h"
#include "VariableInterpolation.h"
class FEAngio;

class Tip;

//! determines the mixture between the position dependent direction and previous segment contribution 
class ContributionMix : public FEMaterialProperty
{
	FECORE_BASE_CLASS(ContributionMix)
public:
	//! constructor
	explicit ContributionMix(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~ContributionMix() {}
	//! return the contribution mix at a given location
	virtual double ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh) = 0;
	//! updates the contribution mix to a given time
	virtual void Update(FEMesh * mesh) {}
};

//! combines the contribution mixes
class ContributionMixManager : public FEMaterialProperty
{
	FECORE_BASE_CLASS(ContributionMixManager)
public:
	//! constructor
	explicit ContributionMixManager(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~ContributionMixManager() {}
	//! return the contribution mix at a given location
	double ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh);
	//! updates the contribution mix manager to a given time
	void Update(FEMesh * mesh) {}
protected:
	DECLARE_FECORE_CLASS()
private:
	std::vector<ContributionMix*>	cm_modifiers;	//!< pointers to elastic materials

};

//! set the contribution mix to a value or a load curve
class PSCPDDContributionMix : public ContributionMix
{
public:
	//! constructor
	explicit PSCPDDContributionMix(FEModel* pfem) : ContributionMix(pfem) {}
	virtual ~PSCPDDContributionMix() {}
	//! return the contribution mix at a given location
	double ApplyModifiers(double dt, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh) override;
	//! updates the contribution mix to a given time
	void Update(FEMesh * mesh) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	double psc_weight = 1.0;
};

//! set the contribution mix to a value or a load curve
class DensFAContributionMix : public ContributionMix
{
public:
	//! constructor
	explicit DensFAContributionMix(FEModel* pfem) : ContributionMix(pfem) {}
	virtual ~DensFAContributionMix() {}
	//! return the contribution mix at a given location
	double ApplyModifiers(double dt, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh) override;
	//! updates the contribution mix to a given time
	void Update(FEMesh* mesh) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	double a0 = 0.1;
	double a_min = 0.05;
	double b = -33.28;
	double c = 0.3;
};

class ProtoContributionMix : public FEMaterialProperty
{
	FECORE_BASE_CLASS(ProtoContributionMix)
public:
	//! constructor
	explicit ProtoContributionMix(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~ProtoContributionMix() {}
	//! return the contribution mix at a given location
	virtual double ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh) = 0;
	//! updates the contribution mix to a given time
	virtual void Update(FEMesh * mesh) {}
};

//! combines the contribution mixes
class ProtoContributionMixManager : public FEMaterialProperty
{
	FECORE_BASE_CLASS(ProtoContributionMixManager)
public:
	//! constructor
	explicit ProtoContributionMixManager(FEModel* pfem) : FEMaterialProperty(pfem) { }
	virtual ~ProtoContributionMixManager() {}
	//! return the contribution mix at a given location
	double ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh);
	//! updates the contribution mix manager to a given time
	void Update(FEMesh * mesh) {}
protected:
	DECLARE_FECORE_CLASS()
private:
	std::vector<ProtoContributionMix*>	proto_cm_modifiers;	//!< pointers to elastic materials
};

//! set the contribution mix to a value or a load curve
class ProtoPSCPDDContributionMix : public ProtoContributionMix
{
public:
	//! constructor
	explicit ProtoPSCPDDContributionMix(FEModel* pfem) : ProtoContributionMix(pfem) {}
	virtual ~ProtoPSCPDDContributionMix() {}
	//! return the contribution mix at a given location
	double ApplyModifiers(double dt, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh) override;
	//! updates the contribution mix to a given time
	void Update(FEMesh * mesh) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	double proto_psc_weight = 1.0;
};