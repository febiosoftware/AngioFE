#pragma once
#include <FECore/FEMaterial.h>
#include <FECore/FESPRProjection.h>
#include "AngioElement.h"
#include "VariableInterpolation.h"
class FEAngio;

class Tip;

//! determines the mixture between the position dependent direction and previous segment contribution 
class ContributionMix : public FEMaterial
{
public:
	//! constructor
	explicit ContributionMix(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~ContributionMix() {}
	//! return the contribution mix at a given location
	virtual double ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh) = 0;
	//! updates the contribution mix to a given time
	virtual void Update(FEMesh * mesh) {}
};

//! combines the contribution mixes
class ContributionMixManager : public FEMaterial
{
public:
	//! constructor
	explicit ContributionMixManager(FEModel* pfem) : FEMaterial(pfem) { AddClassProperty(this, &cm_modifiers, "psc_modifier", FEProperty::Optional); }
	virtual ~ContributionMixManager() {}
	//! return the contribution mix at a given location
	double ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh);
	//! updates the contribution mix manager to a given time
	void Update(FEMesh * mesh) {}
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
	//! parameter list
	DECLARE_FECORE_CLASS();
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
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	double a0 = 0.1;
	double a = 0.26;
	double b = -62.23;
	double c = 0.3402;
};

class ProtoContributionMix : public FEMaterial
{
public:
	//! constructor
	explicit ProtoContributionMix(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~ProtoContributionMix() {}
	//! return the contribution mix at a given location
	virtual double ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh) = 0;
	//! updates the contribution mix to a given time
	virtual void Update(FEMesh * mesh) {}
};

//! combines the contribution mixes
class ProtoContributionMixManager : public FEMaterial
{
public:
	//! constructor
	explicit ProtoContributionMixManager(FEModel* pfem) : FEMaterial(pfem) { AddClassProperty(this, &proto_cm_modifiers, "proto_psc_modifier", FEProperty::Optional); }
	virtual ~ProtoContributionMixManager() {}
	//! return the contribution mix at a given location
	double ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh);
	//! updates the contribution mix manager to a given time
	void Update(FEMesh * mesh) {}
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
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	double proto_psc_weight = 1.0;
};