#pragma once
#include <FECore/FEMaterial.h>
#include "AngioElement.h"

class FEAngio;
class InitialModifier : public FEMaterial
{
public:
	//! constructor
	explicit InitialModifier(FEModel * pfem) :FEMaterial(pfem){}
	//! applies an initial modifier
	virtual void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)=0;
};

class InitialModifierManager: public FEMaterial
{
public:
	//! constructor
	explicit InitialModifierManager(FEModel * pfem) :FEMaterial(pfem) { AddProperty(&initial_modifiers, "initial_modifier"); initial_modifiers.m_brequired = false; }
	//! apply all initial modifiers
	void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio);
private:
	FEVecPropertyT<InitialModifier> initial_modifiers;
};

class FiberRandomizer : public InitialModifier
{
public:
	//! constructor
	explicit FiberRandomizer(FEModel * pfem) : InitialModifier(pfem) {}
	//! apply random fibers to all elements
	void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio) override;
};

class DensityInitializer : public InitialModifier
{
public:
	//! constructor
	explicit DensityInitializer(FEModel * pfem) : InitialModifier(pfem) {}
	//! initialize all densities
	void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio) override;
protected:
	//! parameter list
	DECLARE_PARAMETER_LIST();
private:
	double initial_density = 3.0;
};