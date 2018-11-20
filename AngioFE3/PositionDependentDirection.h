#pragma once
#include <FECore/FEMaterial.h>
#include <FEBioMech/FESPRProjection.h>
#include "AngioElement.h"
#include "VariableInterpolation.h"
class FEAngio;

class Tip;

//! The component of vessel growth that is dependent on position within the mesh
class PositionDependentDirection : public FEMaterial
{
public:
	explicit PositionDependentDirection(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~PositionDependentDirection() {}
	virtual vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continute_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) = 0;
	virtual void Update(FEMesh * mesh, FEAngio* angio) {} //may be used to get values from loadcurves that modify the behavior as a whole
	//! parameter list
	DECLARE_PARAMETER_LIST();
protected:
	double contribution = 1.0;
};

//! Contain the collection of Position Dependent Direction modifiers
class PositionDependentDirectionManager : public FEMaterial
{
public:
	explicit PositionDependentDirectionManager(FEModel* pfem) : FEMaterial(pfem) { AddProperty(&pdd_modifiers, "pdd_modifier"); pdd_modifiers.m_brequired = false; }
	virtual ~PositionDependentDirectionManager() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int buffer, bool& continue_growth, vec3d& tip_dir, double& alpha, FEMesh* mesh, FEAngio* pangio);
	void Update(FEMesh * mesh, FEAngio* angio) {} //may be used to get values from loadcurves that modify the behavior as a whole
private:
	FEVecPropertyT<PositionDependentDirection> pdd_modifiers;
};

//! Implements a position dependent modifier that modifies growth direction based on fiber direction
class FiberPDD : public PositionDependentDirection
{
public:
	explicit FiberPDD(FEModel* pfem) : PositionDependentDirection(pfem) { AddProperty(&interpolation_prop, "interpolation_prop"); }
	virtual ~FiberPDD() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continute_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	void Update(FEMesh * mesh, FEAngio* angio) override;
private:
	FEPropertyT<FEVariableInterpolation> interpolation_prop;
};

//! Implements a position dependent modifier that modifies growth direction if the density gradient is above a given threshold
class ECMDensityGradientPDD : public PositionDependentDirection
{
public:
	explicit ECMDensityGradientPDD(FEModel* pfem) : PositionDependentDirection(pfem) { AddProperty(&interpolation_prop, "interpolation_prop"); }
	virtual ~ECMDensityGradientPDD() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continute_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	void Update(FEMesh * mesh, FEAngio* angio) override {}
	//! parameter list
	DECLARE_PARAMETER_LIST();
private:
	double threshold = 0.00001;//vessels will deflect if above threshold
	bool alpha_override = true;//replace the alpha to have this take over
	FEPropertyT<FEVariableInterpolation> interpolation_prop;
};

//! The replacement for the bouncy boundary condition
class RepulsePDD : public PositionDependentDirection
{
public:
	explicit RepulsePDD(FEModel* pfem) : PositionDependentDirection(pfem) { AddProperty(&interpolation_prop, "interpolation_prop"); }
	virtual ~RepulsePDD() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continute_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	void Update(FEMesh * mesh, FEAngio* angio) override {}
	//! parameter list
	DECLARE_PARAMETER_LIST();
private:
	double threshold = 0.00001;//vessels will deflect if above threshold
	bool alpha_override = true;//replace the alpha to have this take over
	bool grad_threshold = true;//use a gradient to detect areas where repulsion should occour.
	FEPropertyT<FEVariableInterpolation> interpolation_prop;
};

//! Implements a chemical concentration based position direction modifier
class ConcentrationGradientPDD : public PositionDependentDirection
{
public:
	explicit ConcentrationGradientPDD(FEModel* pfem) : PositionDependentDirection(pfem) {}
	virtual ~ConcentrationGradientPDD() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continute_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	void Update(FEMesh * mesh, FEAngio* angio) override {}
	//! parameter list
	DECLARE_PARAMETER_LIST();
private:
	double threshold = 0.00001;//vessels will deflect if above threshold
	bool alpha_override = true;//replace the alpha to have this take over
	int sol_id = 0;
	FEPropertyT<FEVariableInterpolation> interpolation_prop;
};

//! Implements anastamosis as a position dependent direction growth modifier
class AnastamosisPDD : public PositionDependentDirection
{
public:
	explicit AnastamosisPDD(FEModel* pfem) : PositionDependentDirection(pfem) {}
	virtual ~AnastamosisPDD() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continute_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	void Update(FEMesh * mesh, FEAngio* angio) override {}
	Tip * FuseWith(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh, vec3d tip_pos, vec3d tip_dir, int exclude, double radius);
	bool ValidTip(Tip* tip, vec3d tip_dir, FEMesh * mesh);
	Tip * BestInElement(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh, vec3d tip_origin, vec3d tip_dir, int exclude, double& best_distance);
	double distance2(FESolidElement * se, vec3d local_pos, Tip * tip, FEMesh* mesh);
protected:
	//! parameter list
	DECLARE_PARAMETER_LIST();
private:
	double anastamosis_radius = 100;//! Radius at which the tip starts to grow towards another tip
	double fuse_radius = 30;
	double fuse_angle = 0.25;
	bool alpha_override = true;//! Replace the alpha to have this take over
};