#pragma once
#include <FECore/FEMaterial.h>
#include "AngioElement.h"

class FEAngio;
class InitialModifier : public FEMaterial
{
public:
	explicit InitialModifier(FEModel * pfem) :FEMaterial(pfem){}
	virtual void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)=0;
};

class InitialModifierManager: public FEMaterial
{
public:
	explicit InitialModifierManager(FEModel * pfem) :FEMaterial(pfem) { AddProperty(&initial_modifiers, "initial_modifier"); }
	void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio);
private:
	FEVecPropertyT<InitialModifier> initial_modifiers;
};

class FiberRandomizer : public InitialModifier
{
public:
	explicit FiberRandomizer(FEModel * pfem) : InitialModifier(pfem) {}
	void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio) override;
};

class DensityInitializer : public InitialModifier
{
public:
	explicit DensityInitializer(FEModel * pfem) : InitialModifier(pfem) {}
	void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio) override;
protected:
	DECLARE_PARAMETER_LIST();
private:
	double initial_density = 3.0;
};