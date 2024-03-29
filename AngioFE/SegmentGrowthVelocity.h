#pragma once
#include <FECore/FEMaterial.h>
#include <FECore/FESPRProjection.h>
#include "AngioElement.h"
#include "VariableInterpolation.h"
class FEAngio;
class Tip;


//! Get the segment velocity at the given position
class SegmentGrowthVelocity : public FEMaterialProperty
{
	FECORE_BASE_CLASS(SegmentGrowthVelocity)
public:
	//! constructor
	explicit SegmentGrowthVelocity(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~SegmentGrowthVelocity() {}
	//! returns the velocity at the position after it has been modified by this modifier
	virtual double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_elem, double time_shift, FEMesh* mesh) = 0;
	virtual void UpdateScale() = 0;
};

//! return the segment velocity at a given location. A combination of all velocity modifier
class SegmentGrowthVelocityManager : public FEMaterialProperty
{
	FECORE_BASE_CLASS(SegmentGrowthVelocityManager)
public:
	//!constructor
	explicit SegmentGrowthVelocityManager(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~SegmentGrowthVelocityManager() {}
	//! Apply all of the modifier to calculate the velocity at a location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_elem, double time_shift, FEMesh* mesh);
	std::vector<SegmentGrowthVelocity*>	seg_vel_modifiers;	//!< pointers to elastic materials
	//SegmentGrowthVelocity* seg_vel_modifiers;	
protected:
	DECLARE_FECORE_CLASS()
};

//! a fixed or load curve value for velocity
class SegmentVelocityModifier : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SegmentVelocityModifier(FEModel* pfem) : SegmentGrowthVelocity(pfem) {}
	//! modifies the velocity by the velocity over time parameter
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, double time_shift, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
protected:
	DECLARE_FECORE_CLASS()
private:
	double segment_velocity_over_time = 1.0;
};

//! scales the segment velocity based on the ecm density
class SegmentVelocityDensityScaleModifier : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SegmentVelocityDensityScaleModifier(FEModel* pfem) : SegmentGrowthVelocity(pfem) {}
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, double time_shift, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
protected:
	DECLARE_FECORE_CLASS()
private:
	FEVariableInterpolation* interpolation_prop = nullptr;
	vec3d m_density_scale_factor = vec3d(-0.016, 5.1605, 0.5112);
};

class SegmentVelocityRefDensityScaleModifier : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SegmentVelocityRefDensityScaleModifier(FEModel* pfem) : SegmentGrowthVelocity(pfem) {}
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, double time_shift, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
protected:
	DECLARE_FECORE_CLASS()
private:
	FEVariableInterpolation* interpolation_prop = nullptr;
	vec3d m_density_scale_factor = vec3d(-0.016, 5.1605, 0.5112);
};

//! scales the segment velocity based on the ecm density
class SegmentVelocityDensityFAScaleModifier : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SegmentVelocityDensityFAScaleModifier(FEModel* pfem) : SegmentGrowthVelocity(pfem) {}
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, double time_shift, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
protected:
	DECLARE_FECORE_CLASS()
private:
	FEVariableInterpolation* interpolation_prop = nullptr;
	double m_rFA_a = 3.413;
	double m_rFA_b = 1.759;
	double m_rFA_c = 0.6155;
	double m_rFA_d = 0.1;
	double m_rFA_r0 = 2.0;
	double m_rFA_f0 = 0.85;
};

//! scales the segment velocity based on the ecm density
class SegmentVelocity3PModifier : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SegmentVelocity3PModifier(FEModel* pfem) : SegmentGrowthVelocity(pfem) {}
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, double time_shift, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
protected:
	DECLARE_FECORE_CLASS()
private:
	FEVariableInterpolation* interpolation_prop = nullptr;
	double scale = 1.0;
	double threshold = 1.0;
};

//! scales the segment velocity based on the fractional anisotropy
class SegmentVelocityFAModifier : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SegmentVelocityFAModifier(FEModel* pfem) : SegmentGrowthVelocity(pfem) {}
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, double time_shift, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
	vec3d m_density_scale_factor = vec3d(-0.016, 5.1605, 0.5112);
protected:
	DECLARE_FECORE_CLASS()
private:
	FEVariableInterpolation* interpolation_prop = nullptr;
	double scale = 1.0;
	double threshold = 1.0;
};

class SigmoidSegmentVelocity : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SigmoidSegmentVelocity(FEModel* pfem) : SegmentGrowthVelocity(pfem) {}
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, double time_shift, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS()
private:
	double scale = 1.0;
	double a = 100.0;
	double b = 1.3; 
	double c = 5.0;
};

class SigmoidAdjustedSegmentVelocity : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit SigmoidAdjustedSegmentVelocity(FEModel* pfem) : SegmentGrowthVelocity(pfem) {}
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, double time_shift, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS()
private:
	double a = 50.0;
	double b = 2.572;
	double c = 9.368;
	double scale = 1.0;
};

class GompertzSegmentVelocity : public SegmentGrowthVelocity
{
public:
	//!constructor
	explicit GompertzSegmentVelocity(FEModel* pfem) : SegmentGrowthVelocity(pfem) {}
	//! Scales the velocity based on the ecm density at the location
	double ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, double time_shift, FEMesh* mesh) override;
	//! performs initialization
	bool Init() override;
	void UpdateScale() override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS()
private:
	double scale = 1.0;
	double a = 284.0;
	double b = 0.5;
	double c = 1.0;
	double d = 5.0;
};