#pragma once
#include <FECore/FEMaterial.h>
#include <FECore/FESPRProjection.h>
#include "AngioElement.h"
#include "VariableInterpolation.h"
class FEAngio;
class Tip;


//! Get the segment velocity at the given position
class SegmentGrowthVelocity : public FEMaterial
{
public:
	//! constructor
	explicit SegmentGrowthVelocity(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~SegmentGrowthVelocity() {}
	//! returns the velocity at the position after it has been modified by this modifier
	virtual double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_elem, FEMesh* mesh) = 0;
	virtual void UpdateScale() = 0;
};

//! return the segment velocity at a given location. A combination of all velocity modifier
class SegmentGrowthVelocityManager : public FEMaterial
{
public:
	//!constructor
	explicit SegmentGrowthVelocityManager(FEModel* pfem) : FEMaterial(pfem) { AddClassProperty(this, &seg_vel_modifiers, "velocity_modifier", FEProperty::Optional); }
	virtual ~SegmentGrowthVelocityManager() {}
	//! Apply all of the modifier to calculate the velocity at a location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_elem, FEMesh* mesh);
	std::vector<SegmentGrowthVelocity*>	seg_vel_modifiers;	//!< pointers to elastic materials
	//SegmentGrowthVelocity* seg_vel_modifiers;	
};

//! a fixed or load curve value for velocity
class SegmentVelocityModifier : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SegmentVelocityModifier(FEModel* pfem) : SegmentGrowthVelocity(pfem) {}
	//! modifies the velocity by the velocity over time parameter
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	double segment_velocity_over_time = 1;
};

//! scales the segment velocity based on the ecm density
class SegmentVelocityDensityScaleModifier : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SegmentVelocityDensityScaleModifier(FEModel* pfem) : SegmentGrowthVelocity(pfem) { AddClassProperty(this, &interpolation_prop, "interpolation_prop"); }
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	FEVariableInterpolation* interpolation_prop;
	vec3d m_density_scale_factor = vec3d(-0.016, 5.1605, 0.5112);
};

class SegmentVelocityRefDensityScaleModifier : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SegmentVelocityRefDensityScaleModifier(FEModel* pfem) : SegmentGrowthVelocity(pfem) { AddClassProperty(this, &interpolation_prop, "interpolation_prop"); }
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	FEVariableInterpolation* interpolation_prop;
	vec3d m_density_scale_factor = vec3d(-0.016, 5.1605, 0.5112);
};

//! scales the segment velocity based on the ecm density
class SegmentVelocity3PModifier : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SegmentVelocity3PModifier(FEModel* pfem) : SegmentGrowthVelocity(pfem) { AddClassProperty(this, &interpolation_prop, "interpolation_prop"); }
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	FEVariableInterpolation* interpolation_prop;
	double scale = 1;
	double threshold = 1;
};

//! scales the segment velocity based on the fractional anisotropy
class SegmentVelocityFAModifier : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SegmentVelocityFAModifier(FEModel* pfem) : SegmentGrowthVelocity(pfem) { AddClassProperty(this, &interpolation_prop, "interpolation_prop"); }
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	FEVariableInterpolation* interpolation_prop;
	double scale = 1;
	double threshold = 1;
};

class SigmoidSegmentVelocity : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SigmoidSegmentVelocity(FEModel* pfem) : SegmentGrowthVelocity(pfem) {}
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	double scale = 1;
	double a = 100;
	double b = 1; 
	double c = 5;
};

class GompertzSegmentVelocity : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit GompertzSegmentVelocity(FEModel* pfem) : SegmentGrowthVelocity(pfem) {}
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	double scale = 1;
	double a = 284;
	double b = 0.5;
	double c = 1;
	double d = 5;
};