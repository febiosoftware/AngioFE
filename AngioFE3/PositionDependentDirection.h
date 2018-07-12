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
	virtual void Update(FEMesh * mesh, FEAngio* angio) {} //may be used to get values from loadcurves that modify the behavior as a whole
	DECLARE_PARAMETER_LIST();
protected:
	double contribution = 1.0;
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

class FiberPDD : public PositionDependentDirection
{
public:
	explicit FiberPDD(FEModel* pfem) : PositionDependentDirection(pfem) { AddProperty(&interpolation_prop, "interpolation_prop"); }
	virtual ~FiberPDD() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, double& alpha, FEMesh* mesh, FEAngio* pangio) override;
	void Update(FEMesh * mesh, FEAngio* angio) override;
private:
	FEPropertyT<FEVariableInterpolation> interpolation_prop;
};

class ECMDensityGradientPDD : public PositionDependentDirection
{
public:
	explicit ECMDensityGradientPDD(FEModel* pfem) : PositionDependentDirection(pfem) { AddProperty(&interpolation_prop, "interpolation_prop"); }
	virtual ~ECMDensityGradientPDD() {}
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, double& alpha, FEMesh* mesh, FEAngio* pangio) override;
	void Update(FEMesh * mesh, FEAngio* angio) override {}
	DECLARE_PARAMETER_LIST();
private:
	double threshold = 0.00001;//vessels will deflect if above threshold
	bool alpha_override = true;//replace the alpha to have this take over
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
};