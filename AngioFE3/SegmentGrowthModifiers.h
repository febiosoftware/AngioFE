#pragma once
#include <FECore/FEMaterial.h>
#include <FEBioMech/FESPRProjection.h>
#include "AngioElement.h"
#include "VariableInterpolation.h"
class FEAngio;

class Tip;

//the component of vessel growth that is dependent on position within the mesh
class PositionDependentDirection : public FEMaterial
{
public:
	explicit PositionDependentDirection(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~PositionDependentDirection() {}
	virtual vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, double& alpha, FEMesh* mesh, FEAngio* pangio) = 0;
	virtual void Update(FEMesh * mesh, FEAngio* angio){} //may be used to get values from loadcurves that modify the behavior as a whole
};

class PositionDependentDirectionManager : public FEMaterial
{
public:
	explicit PositionDependentDirectionManager(FEModel* pfem) : FEMaterial(pfem) { AddProperty(&pdd_modifiers, "pdd_modifier"); pdd_modifiers.m_brequired = false; }
	virtual ~PositionDependentDirectionManager() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, double& alpha, FEMesh* mesh, FEAngio* pangio);
	void Update(FEMesh * mesh, FEAngio* angio) {} //may be used to get values from loadcurves that modify the behavior as a whole
private:
	FEVecPropertyT<PositionDependentDirection> pdd_modifiers;
};

//The component of the mixture that represents the direction of the conttribution from the previous segment
class PreviousSegmentContribution : public FEMaterial
{
public:
	explicit PreviousSegmentContribution(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~PreviousSegmentContribution() {}
	virtual vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh)=0;
	virtual void Update(FEMesh * mesh) {}
};

//The component of the mixture that represents the direction of the conttribution from the previous segment
class PreviousSegmentContributionManager : public FEMaterial
{
public:
	explicit PreviousSegmentContributionManager(FEModel* pfem) : FEMaterial(pfem) { AddProperty(&psc_modifiers, "psc_modifier"); psc_modifiers.m_brequired = false; }
	virtual ~PreviousSegmentContributionManager() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh);
	void Update(FEMesh * mesh) {}
private:
	FEVecPropertyT<PreviousSegmentContribution> psc_modifiers;
};

//determines the mixture between the position dependent direction and previous segment contribution 
class ContributionMix : public FEMaterial
{
public:
	explicit ContributionMix(FEModel* pfem) : FEMaterial(pfem){}
	virtual ~ContributionMix(){}
	virtual double ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh)=0;
	virtual void Update(FEMesh * mesh) {}
};

class ContributionMixManager : public FEMaterial
{
public:
	explicit ContributionMixManager(FEModel* pfem) : FEMaterial(pfem) { AddProperty(&cm_modifiers, "psc_modifier"); cm_modifiers.m_brequired = false; }
	virtual ~ContributionMixManager() {}
	double ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh);
	void Update(FEMesh * mesh) {}
private:
	FEVecPropertyT<ContributionMix> cm_modifiers;
};

class PSCPDDContributionMix : public ContributionMix
{
public:
	explicit PSCPDDContributionMix(FEModel* pfem) : ContributionMix(pfem) {}
	virtual ~PSCPDDContributionMix() {}
	double ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh) override;
	void Update(FEMesh * mesh) override;
protected:
	DECLARE_PARAMETER_LIST();
private:
	double psc_weight = 0.5;
};

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

class FiberPDD : public PositionDependentDirection
{
public:
	explicit FiberPDD(FEModel* pfem) : PositionDependentDirection(pfem) { AddProperty(&interpolation_prop, "interpolation_prop"); }
	virtual ~FiberPDD() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, double& alpha, FEMesh* mesh, FEAngio* pangio) override;
	void Update(FEMesh * mesh, FEAngio* angio) override;
private:
	double contribution = 1.0;
	FEPropertyT<FEVariableInterpolation> interpolation_prop;
};

class ECMDensityGradientPDD : public PositionDependentDirection
{
public:
	explicit ECMDensityGradientPDD(FEModel* pfem) : PositionDependentDirection(pfem) { AddProperty(&interpolation_prop, "interpolation_prop"); }
	virtual ~ECMDensityGradientPDD() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, double& alpha, FEMesh* mesh, FEAngio* pangio) override;
	void Update(FEMesh * mesh, FEAngio* angio) override{}\
	DECLARE_PARAMETER_LIST();
private:
	double contribution = 1.0;
	double threshold = 0.00001;//vessels will deflect if above threshold
	FEPropertyT<FEVariableInterpolation> interpolation_prop;
};

//
class ConcentrationPDD : public PositionDependentDirection
{
public:
	explicit ConcentrationPDD(FEModel* pfem) : PositionDependentDirection(pfem) {}
	virtual ~ConcentrationPDD() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, double& alpha, FEMesh* mesh, FEAngio* pangio) override;
};

class DensityGradientPDD : public PositionDependentDirection
{
public:
	explicit DensityGradientPDD(FEModel* pfem) : PositionDependentDirection(pfem) {}
	virtual ~DensityGradientPDD() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, double& alpha, FEMesh* mesh, FEAngio* pangio) override;
};

class AnastomosisPDD : public PositionDependentDirection
{
public:
	explicit AnastomosisPDD(FEModel* pfem) : PositionDependentDirection(pfem) {}
	virtual ~AnastomosisPDD() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, double& alpha, FEMesh* mesh, FEAngio* pangio) override;
protected:
	DECLARE_PARAMETER_LIST();
private:
	double anastomosis_radius = 300;
	double contribution = 1.0;
};

class PreviousSegmentPSC : public PreviousSegmentContribution
{
	//overrides the previous segment contribution with the previous segment direction
public:
	explicit PreviousSegmentPSC(FEModel* pfem) : PreviousSegmentContribution(pfem) {}
	virtual ~PreviousSegmentPSC() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh) override;
};