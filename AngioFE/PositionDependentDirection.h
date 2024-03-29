#pragma once
#include <FECore/FEMaterial.h>
#include <FECore/FESPRProjection.h>
#include "AngioElement.h"
#include "VariableInterpolation.h"
class FEAngio;

class Tip;

//! The component of vessel growth that is dependent on position within the mesh
class PositionDependentDirection : public FEMaterialProperty
{
	FECORE_BASE_CLASS(PositionDependentDirection)
public:
	//! constructor
	explicit PositionDependentDirection(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~PositionDependentDirection() {}
	//! return the direction given by this component
	virtual vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) = 0;
	//! may be used to get values from loadcurves that modify the behavior as a whole
	virtual void Update(FEMesh * mesh, FEAngio* angio) {} 
protected:
	//! how this modifier is mixed with the previous direction
	double contribution = 1.0;
	DECLARE_FECORE_CLASS()
};

//! Contain the collection of Position Dependent Direction modifiers
class PositionDependentDirectionManager : public FEMaterialProperty
{
	FECORE_BASE_CLASS(PositionDependentDirectionManager)
public:
	//! constructor
	explicit PositionDependentDirectionManager(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~PositionDependentDirectionManager() {}
	//! return the direction given by all direction modifiers
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int buffer, bool& continue_growth, vec3d& tip_dir, double& alpha, FEMesh* mesh, FEAngio* pangio);
	//! may be used to get values from loadcurves that modify the behavior as a whole
	void Update(FEMesh * mesh, FEAngio* angio); 
protected:
	DECLARE_FECORE_CLASS()
private:
	std::vector<PositionDependentDirection*>	pdd_modifiers;	//!< pointers to elastic materials
};

//! Implements a position dependent modifier that modifies growth direction based on fiber direction
class FiberPDD : public PositionDependentDirection
{
public:
	//! constructor
	explicit FiberPDD(FEModel* pfem) : PositionDependentDirection(pfem) {}
	virtual ~FiberPDD() {}
	//! return the direction given by the fibers at this location
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	//! may be used to get values from loadcurves that modify the behavior as a whole
	void Update(FEMesh * mesh, FEAngio* angio) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	FEVariableInterpolation* interpolation_prop = nullptr;
};

//! Implements a position dependent modifier that modifies growth direction based on fiber direction
class FractionalAnisotropyPDD : public PositionDependentDirection
{
public:
	//! constructor
	explicit FractionalAnisotropyPDD(FEModel* pfem) : PositionDependentDirection(pfem) {}
	virtual ~FractionalAnisotropyPDD() {}
	//! return the direction given by the fibers at this location
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	//! may be used to get values from loadcurves that modify the behavior as a whole
	void Update(FEMesh * mesh, FEAngio* angio) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	bool alpha_override = false;// replace alpha with the override
	double efd_exp = 1.0;
};

class LagrangePStrainPDD : public PositionDependentDirection
{
public:
	//! constructor
	explicit LagrangePStrainPDD(FEModel* pfem) : PositionDependentDirection(pfem) {}
	virtual ~LagrangePStrainPDD() {}
	//! return the direction given by the fibers at this location
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	//! may be used to get values from loadcurves that modify the behavior as a whole
	//void Update(FEMesh * mesh, FEAngio* angio) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	FEVariableInterpolation* interpolation_prop = nullptr;
	double beta = 0.5;
};

//! Implements a position dependent modifier that modifies growth direction if the density gradient is above a given threshold
class ECMDensityGradientPDD : public PositionDependentDirection
{
public:
	//! constructor
	explicit ECMDensityGradientPDD(FEModel* pfem) : PositionDependentDirection(pfem) {}
	virtual ~ECMDensityGradientPDD() {}
	//! return the direction given by the ecm density gradient
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	//! may be used to get values from loadcurves that modify the behavior as a whole
	void Update(FEMesh * mesh, FEAngio* angio) override {}
	//! parameter list
protected:
	DECLARE_FECORE_CLASS()
private:
	double threshold = 1.0e-5; //vessels will deflect if above threshold
	bool alpha_override = false; //replace the alpha to have this take over
	FEVariableInterpolation* interpolation_prop = nullptr;
};

//! The replacement for the bouncy boundary condition
class RepulsePDD : public PositionDependentDirection
{
public:
	//! constructor
	explicit RepulsePDD(FEModel* pfem) : PositionDependentDirection(pfem){}
	virtual ~RepulsePDD() {}
	//! return the direction given by the repulse component
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	//! may be used to get values from loadcurves that modify the behavior as a whole
	void Update(FEMesh * mesh, FEAngio* angio) override {}
protected:
	DECLARE_FECORE_CLASS()
private:
	double threshold = 1.0; //vessels will deflect if above threshold
	bool alpha_override = true; //replace the alpha to have this take over
	bool grad_threshold = false; //use a gradient to detect areas where repulsion should occur.
	FEVariableInterpolation* interpolation_prop = nullptr;
};

//! Implements a chemical concentration based position direction modifier
//! SL: Currently experimental; not thoroughly tested.
class ConcentrationGradientPDD : public PositionDependentDirection
{
public:
	//! constructor
	explicit ConcentrationGradientPDD(FEModel* pfem) : PositionDependentDirection(pfem) { }
	virtual ~ConcentrationGradientPDD() {}
	//! return the direction given by the concentration gradient
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	//! may be used to get values from loadcurves that modify the behavior as a whole
	void Update(FEMesh * mesh, FEAngio* angio) override {}
protected:
	DECLARE_FECORE_CLASS()
private:
	//! SL: will need to rethink this default
	double threshold = 1.0e-5; //vessels will deflect if above threshold
	bool alpha_override = false; //replace the alpha to have this take over
	int sol_id = 0;
	FEVariableInterpolation* interpolation_prop = nullptr;
};

//! Implements a chemical concentration based position direction modifier
//! SL: Currently experimental; not thoroughly tested.
class FisherConcentrationGradientPDD : public PositionDependentDirection
{
public:
	//! constructor
	explicit FisherConcentrationGradientPDD(FEModel* pfem) : PositionDependentDirection(pfem) { }
	virtual ~FisherConcentrationGradientPDD() {}
	//! return the direction given by the concentration gradient
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	//! may be used to get values from loadcurves that modify the behavior as a whole
	void Update(FEMesh * mesh, FEAngio* angio) override {}
protected:
	DECLARE_FECORE_CLASS()
private:
	//! SL: will need to rethink this default
	double threshold = 1.0e-5; //vessels will deflect if above threshold
	bool alpha_override = false; //replace the alpha to have this take over
	int sol_id = 0;
	FEVariableInterpolation* interpolation_prop = nullptr;
};

//! Implements anastamosis as a position dependent direction growth modifier
class AnastamosisPDD : public PositionDependentDirection
{
public:
	//! constructor
	explicit AnastamosisPDD(FEModel* pfem) : PositionDependentDirection(pfem) { }
	virtual ~AnastamosisPDD() {}
	//! return the direction given by the anastamosis modifier
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	//! may be used to get values from loadcurves that modify the behavior as a whole
	void Update(FEMesh * mesh, FEAngio* angio) override {}
	//! returns the tip that the position should fuse with
	Tip * FuseWith(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh, vec3d tip_pos, vec3d tip_dir, int exclude, double radius);
	//! return whether or not a tip can be fused with
	bool ValidTip(Tip* tip, vec3d tip_dir, FEMesh * mesh);
	//! returns the best tip within an element
	Tip * BestInElement(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh, vec3d tip_origin, vec3d tip_dir, int exclude, double& best_distance);
	//! return the distance squared between a tip and a local position
	double distance2(FESolidElement * se, vec3d local_pos, Tip * tip, FEMesh* mesh);
protected:
	DECLARE_FECORE_CLASS()
private:
	double anastamosis_radius = 100.0; //! Radius at which the tip starts to grow towards another tip
	double fuse_radius = 30.0;
	double fuse_angle = 0.25;
	bool alpha_override = true; //! Replace the alpha to have this take over
};

//! The component of vessel growth that is dependent on position within the mesh
class ProtoPositionDependentDirection : public FEMaterialProperty
{
	FECORE_BASE_CLASS(ProtoPositionDependentDirection)

public:
	//! constructor
	explicit ProtoPositionDependentDirection(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~ProtoPositionDependentDirection() {}
	//! return the direction given by this component
	virtual vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) = 0;
	//! may be used to get values from loadcurves that modify the behavior as a whole
	virtual void Update(FEMesh * mesh, FEAngio* angio) {}
	//! parameter list
protected:
	DECLARE_FECORE_CLASS()
	//! how this modifier is mixed with the previous direction
	double proto_contribution = 1.0;
};

//! Contain the collection of Position Dependent Direction modifiers
class ProtoPositionDependentDirectionManager : public FEMaterialProperty
{
	FECORE_BASE_CLASS(ProtoPositionDependentDirectionManager)

public:
	//! constructor
	explicit ProtoPositionDependentDirectionManager(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~ProtoPositionDependentDirectionManager() {}
	//! return the direction given by all direction modifiers
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int buffer, bool& continue_growth, vec3d& tip_dir, double& alpha, FEMesh* mesh, FEAngio* pangio);
	//! may be used to get values from loadcurves that modify the behavior as a whole
	void Update(FEMesh * mesh, FEAngio* angio);
protected:
	DECLARE_FECORE_CLASS()
private:
	std::vector<ProtoPositionDependentDirection*>	proto_pdd_modifiers;	//!< pointers to elastic materials
};

//! Implements a position dependent modifier that modifies growth direction based on fiber direction
class ProtoFiberPDD : public ProtoPositionDependentDirection
{
public:
	//! constructor
	explicit ProtoFiberPDD(FEModel* pfem) : ProtoPositionDependentDirection(pfem) {}
	virtual ~ProtoFiberPDD() {}
	//! return the direction given by the fibers at this location
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	//! may be used to get values from loadcurves that modify the behavior as a whole
	void Update(FEMesh * mesh, FEAngio* angio) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	FEVariableInterpolation* interpolation_prop = nullptr;
	double proto_alpha = 1.0;
};

//! Implements a position dependent modifier that modifies growth direction based on fiber direction
class ProtoFractionalAnisotropyPDD : public ProtoPositionDependentDirection
{
public:
	//! constructor
	explicit ProtoFractionalAnisotropyPDD(FEModel* pfem) : ProtoPositionDependentDirection(pfem) {}
	virtual ~ProtoFractionalAnisotropyPDD() {}
	//! return the direction given by the fibers at this location
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	//! may be used to get values from loadcurves that modify the behavior as a whole
	void Update(FEMesh * mesh, FEAngio* angio) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	FEVariableInterpolation* interpolation_prop = nullptr;
	bool alpha_override = false;// replace alpha with the override
	double proto_alpha = 1.0;
	mat3ds proto_efd = mat3ds(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
	double proto_efd_exp = 1.0;
};

//! The replacement for the bouncy boundary condition
class ProtoRepulsePDD : public ProtoPositionDependentDirection
{
public:
	//! constructor
	explicit ProtoRepulsePDD(FEModel* pfem) : ProtoPositionDependentDirection(pfem) {}
	virtual ~ProtoRepulsePDD() {}
	//! return the direction given by the repulse component
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	//! may be used to get values from loadcurves that modify the behavior as a whole
	void Update(FEMesh * mesh, FEAngio* angio) override {}
protected:
	DECLARE_FECORE_CLASS()
private:
	double threshold = 1.0; //vessels will deflect if above threshold
	bool alpha_override = true; //replace the alpha to have this take over
	bool grad_threshold = false; //use a gradient to detect areas where repulsion should occur.
	FEVariableInterpolation* interpolation_prop = nullptr;
};

//! Implements anastamosis as a position dependent direction growth modifier
class ProtoAnastamosisPDD : public ProtoPositionDependentDirection
{
public:
	//! constructor
	explicit ProtoAnastamosisPDD(FEModel* pfem) : ProtoPositionDependentDirection(pfem) { }
	virtual ~ProtoAnastamosisPDD() {}
	//! return the direction given by the anastamosis modifier
	vec3d ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) override;
	//! may be used to get values from loadcurves that modify the behavior as a whole
	void Update(FEMesh * mesh, FEAngio* angio) override {}
	//! returns the tip that the position should fuse with
	Tip * FuseWith(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh, vec3d tip_pos, vec3d tip_dir, int exclude, double radius);
	//! return whether or not a tip can be fused with
	bool ValidTip(Tip* tip, vec3d tip_dir, FEMesh * mesh);
	//! returns the best tip within an element
	Tip * BestInElement(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh, vec3d tip_origin, vec3d tip_dir, int exclude, double& best_distance);
	//! return the distance squared between a tip and a local position
	double distance2(FESolidElement * se, vec3d local_pos, Tip * tip, FEMesh* mesh);
protected:
	DECLARE_FECORE_CLASS()
private:
	double anastamosis_radius = 100.0; //! Radius at which the tip starts to grow towards another tip
	double fuse_radius = 30.0;
	double fuse_angle = 0.25;
	bool alpha_override = true; //! Replace the alpha to have this take over
};