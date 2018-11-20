#pragma once
#include <FECore/FEMaterial.h>
#include <FEBioMech/FESPRProjection.h>
#include "AngioElement.h"
#include "VariableInterpolation.h"
class FEAngio;
class Tip;


// Get the length over one unit of time at the given position
class SegmentGrowthVelocity : public FEMaterial
{
public:
	explicit SegmentGrowthVelocity(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~SegmentGrowthVelocity() {}
	virtual double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_elem, FEMesh* mesh) = 0;
};

class SegmentGrowthVelocityManager : public FEMaterial
{
public:
	explicit SegmentGrowthVelocityManager(FEModel* pfem) : FEMaterial(pfem) { AddProperty(&seg_vel_modifiers, "velocity_modifier"); seg_vel_modifiers.m_brequired = false; }
	virtual ~SegmentGrowthVelocityManager() {}
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_elem, FEMesh* mesh);
private:
	FEVecPropertyT<SegmentGrowthVelocity> seg_vel_modifiers;
};

class SegmentVelocityModifier : public SegmentGrowthVelocity
{
public:
	explicit SegmentVelocityModifier(FEModel* pfem) : SegmentGrowthVelocity(pfem) {}
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
protected:
	DECLARE_PARAMETER_LIST();
private:
	double segment_velocity_over_time = 1;
};

class SegmentVelocityDensityScaleModifier : public SegmentGrowthVelocity
{
public:
	explicit SegmentVelocityDensityScaleModifier(FEModel* pfem) : SegmentGrowthVelocity(pfem) { AddProperty(&interpolation_prop, "interpolation_prop"); }
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh) override;
	bool Init() override;
protected:
	DECLARE_PARAMETER_LIST();
private:
	FEPropertyT<FEVariableInterpolation> interpolation_prop;
	vec3d m_density_scale_factor = vec3d(-0.016, 5.1605, 0.5112);
};