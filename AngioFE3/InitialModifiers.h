#pragma once
#include <FECore/FEMaterial.h>
#include "AngioElement.h"

class FEAngio;
//! Anything that should be applied to all elements within an angio material, applied before nodedatainterpolation values
class InitialModifier : public FEMaterial
{
public:
	//! constructor
	explicit InitialModifier(FEModel * pfem) :FEMaterial(pfem){}
	//! applies an initial modifier
	virtual void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)=0;
};

//! applies all initial modifiers
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

//! sets fibers to a random orientation
class FiberRandomizer : public InitialModifier
{
public:
	//! constructor
	explicit FiberRandomizer(FEModel * pfem) : InitialModifier(pfem) {}
	//! apply random fibers to all elements
	void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio) override;
};

//! sets fibers to an EFD distribution
class DiscreteFiberEFDRandomizer : public InitialModifier
{
public:
	//! constructor
	explicit DiscreteFiberEFDRandomizer(FEModel * pfem) : InitialModifier(pfem) {}
	//! apply random fibers to all elements
	void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio) override;
protected:
	//! parameter list
	DECLARE_PARAMETER_LIST();
private:
	vec3d efd_spa = vec3d(1, 1, 1);
	vec3d efd_axes_a = vec3d(1, 0, 0);
	vec3d efd_axes_b = vec3d(0, 1, 0);
	vec3d efd_axes_c = vec3d(0, 0, 1);
};

class EFDFiberInitializer : public InitialModifier 
{
public:
	//! constructor
	explicit EFDFiberInitializer(FEModel * pfem) : InitialModifier(pfem) {}
	//! apply Isotropic EFD to each mdoel
	void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio) override;
protected:
	//! parameter list
	DECLARE_PARAMETER_LIST();
private:
	vec3d initial_spa = vec3d(1,1,1);
	vec3d initial_axes_a = vec3d(1, 0, 0);
	vec3d initial_axes_b = vec3d(0, 1, 0);
	vec3d initial_axes_c = vec3d(0, 0, 1);
};

//! sets the ecm density within a material
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

class RepulseInitializer : public InitialModifier
{
public:
	//! constructor
	explicit RepulseInitializer(FEModel * pfem) : InitialModifier(pfem) {}
	//! initialize all densities
	void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio) override;
protected:
	//! parameter list
	DECLARE_PARAMETER_LIST();
private:
	double initial_repulse_value = 0.0;
};