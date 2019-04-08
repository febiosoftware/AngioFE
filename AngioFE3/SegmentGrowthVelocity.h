#pragma once
#include <FECore/FEMaterial.h>
#include <FEBioMech/FESPRProjection.h>
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
};

//! return the segment velocity at a given location. A combination of all velocity modifier
class SegmentGrowthVelocityManager : public FEMaterial
{
public:
	//!constructor
	explicit SegmentGrowthVelocityManager(FEModel* pfem) : FEMaterial(pfem) { AddProperty(&seg_vel_modifiers, "velocity_modifier"); seg_vel_modifiers.m_brequired = false; }
	virtual ~SegmentGrowthVelocityManager() {}
	//! Apply all of the modifier to calculate the velocity at a location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_elem, FEMesh* mesh);
private:
	FEVecPropertyT<SegmentGrowthVelocity> seg_vel_modifiers;
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
protected:
	//! parameter list
	DECLARE_PARAMETER_LIST();
private:
	double segment_velocity_over_time = 1;
};

//! scales the segment velocity based on the ecm density
class SegmentVelocityDensityScaleModifier : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SegmentVelocityDensityScaleModifier(FEModel* pfem) : SegmentGrowthVelocity(pfem) { AddProperty(&interpolation_prop, "interpolation_prop"); }
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
protected:
	//! parameter list
	DECLARE_PARAMETER_LIST();
private:
	FEPropertyT<FEVariableInterpolation> interpolation_prop;
	vec3d m_density_scale_factor = vec3d(-0.016, 5.1605, 0.5112);
};

class SegmentVelocityRefDensityScaleModifier : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SegmentVelocityRefDensityScaleModifier(FEModel* pfem) : SegmentGrowthVelocity(pfem) { AddProperty(&interpolation_prop, "interpolation_prop"); }
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
protected:
	//! parameter list
	DECLARE_PARAMETER_LIST();
private:
	FEPropertyT<FEVariableInterpolation> interpolation_prop;
	vec3d m_density_scale_factor = vec3d(-0.016, 5.1605, 0.5112);
};