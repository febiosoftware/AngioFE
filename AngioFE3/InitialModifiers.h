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
	explicit InitialModifierManager(FEModel * pfem) :FEMaterial(pfem) { AddClassProperty(this, &initial_modifiers, "initial_modifier", FEProperty::Optional); }
	//! apply all initial modifiers
	void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio);
private:
	std::vector<InitialModifier*>	initial_modifiers;	//!< pointers to initial modifiers
	//InitialModifier* initial_modifiers;
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
	DECLARE_FECORE_CLASS();
private:
	vec3d efd_spd = vec3d(1, 1, 1);
	vec3d efd_axes_a = vec3d(1, 0, 0);
	vec3d efd_axes_b = vec3d(0, 1, 0);
	vec3d efd_axes_c = vec3d(0, 0, 1);
};

class DiscreteFiberEFDMatRandomizer : public InitialModifier
{
public:
	//! constructor
	explicit DiscreteFiberEFDMatRandomizer(FEModel * pfem) : InitialModifier(pfem) {}
	//! apply Isotropic EFD to each model
	void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio) override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	mat3ds m_SPD;
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
	DECLARE_FECORE_CLASS();
private:
	mat3ds m_SPDa;
};

class EFDMatPointFiberInitializer : public InitialModifier
{
public:
	//! constructor
	explicit EFDMatPointFiberInitializer(FEModel * pfem) : InitialModifier(pfem) {}
	//! apply Isotropic EFD to each mdoel
	void ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio) override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	mat3ds m_SPD;
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
	DECLARE_FECORE_CLASS();
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
	DECLARE_FECORE_CLASS();
private:
	double initial_repulse_value = 0.0;
};