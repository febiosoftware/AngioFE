#pragma once
#include <FECore/FEMaterial.h>
#include "AngioElement.h"
#include "FEProbabilityDistribution.h"
#include "VariableInterpolation.h"

class FEAngio;

//! Base class, derived classes store the infomation needed by the brancher on a per element basis
class BranchInfo
{
public:
	virtual ~BranchInfo(){}
};

//! Implements a class that will return the zentih angle for a given branch
class ZenithAngle : public FEMaterial
{
public:
	//! Constructor for zenith angle
	ZenithAngle(FEModel* pfem) : FEMaterial(pfem) {}
	//! Gives the zenith angle based on the parameters
	virtual double GetZenithAngle(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element) = 0;
	//! performs any modifications the zenith angle as time changes
	virtual void TimeStepUpdate(double current_time) = 0;
};

//! Implements a class that will return the azimuth angle for a given branch
class AzimuthAngle :public FEMaterial
{
public:
	//! Constructor for azimuth angle
	AzimuthAngle(FEModel* pfem) : FEMaterial(pfem) {}
	//! Gives the zenith angle based on the parameters
	virtual double GetAzimuthAngle(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element) = 0;
	//! performs any modifications the zenith angle as time changes
	virtual void TimeStepUpdate(double current_time) = 0;
};

//! Implements a probability distribution to determine the zenith angle
class ZenithAngleProbabilityDistribution : public ZenithAngle
{
public:
	//! constructor for class
	ZenithAngleProbabilityDistribution(FEModel* pfem) : ZenithAngle(pfem) { AddClassProperty(this, &angle, "angle"); }
	//! returns the zenith angle for a branch based on a probability distribution
	double GetZenithAngle(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element) override;
	//! performs initialization
	bool Init() override { return angle->Init(); }
	//! Updates the zenith angle to the given time step(may ajust probabilities based on time)
	void TimeStepUpdate(double current_time) override { angle->TimeStepUpdate(current_time); }

private:
	FEProbabilityDistribution* angle = nullptr;
};

//! Implements a probability distribution to determine the azimuth angle
class AzimuthAngleProbabilityDistribution :public AzimuthAngle
{
public:
	//! constructor for class
	AzimuthAngleProbabilityDistribution(FEModel* pfem) : AzimuthAngle(pfem) { AddClassProperty(this, &angle, "angle", FEProperty::Required); }
	//! Returns the azimuth angle based on a probability distribution
	double GetAzimuthAngle(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element) override;
	//! performs initialization
	bool Init() override {return angle->Init();}
	//! Updates the zenith angle to the given time step(may ajust probabilities based on time)
	void TimeStepUpdate(double current_time) override { angle->TimeStepUpdate(current_time); }
private:
	FEProbabilityDistribution* angle = nullptr;
};

//! A Branch Policy determines where and when branches occur
class BranchPolicy :public FEMaterial
{
public:
	//! constructor for class
	BranchPolicy(FEModel* pfem) : FEMaterial(pfem) { 
		AddClassProperty(this, &azimuth_angle, "azimuth_angle"); 
		AddClassProperty(this, &zenith_angle, "zenith_angle");
		AddClassProperty(this, &interpolation_prop, "interpolation_prop");
	}
	virtual ~BranchPolicy(){};
	//! adds the branches for a given element
	//! setup any data structures that are needed on a per element basis
	virtual void SetupBranchInfo(AngioElement * angio_elem) = 0;
	//! performs initialization
	bool Init() override { return azimuth_angle->Init() && zenith_angle->Init(); }
	//! Updates the branch policy to the given time step
	virtual void TimeStepUpdate(double current_time){ azimuth_angle->TimeStepUpdate(current_time); zenith_angle->TimeStepUpdate(current_time);}
	//! Return the direction of a branch
	void AddBranchTipEFD(AngioElement * angio_element, vec3d local_pos, Segment* parent_seg, double start_time, int vessel_id, int buffer_index, FEMesh* mesh);
	//! Return the direction of a branch
	vec3d GetBranchDirectionEFD(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element, FEMesh* mesh);
private:
	//angles are in radians
	AzimuthAngle* azimuth_angle = nullptr;
	ZenithAngle* zenith_angle = nullptr;
	FEVariableInterpolation* interpolation_prop = nullptr;
};

//! The data needed for a future branch
class FutureBranch
{
public:
	//! constructor for class
	explicit FutureBranch(vec3d local_pos, Segment * parent, double start_time) : _local_pos(local_pos), _parent(parent), _start_time(start_time) { assert(parent); }
	//! return the local position
	vec3d _local_pos;
	//! the segment this grew from
	Segment * _parent = nullptr;
	//! the time at which the branch starts to grow
	double _start_time;
};

//! Info needed by delayed branchin on a per element 
class DelayBranchInfo : public BranchInfo
{
public:
	//! the length before a branch occurs
	double length_to_branch = 0.0;
	//! all branches which may occur in the future
	std::list<FutureBranch> future_branches;
};

//! As growth occurs length to branch is calculated, once length to branch is hit a time to emerge is sampled from a probability distribution.Length to branch is calculated only when the last value of length to branch has been hit
class DelayedBranchingPolicyEFD :public BranchPolicy
{
public:
	//! constructor for class
	DelayedBranchingPolicyEFD(FEModel* pfem) : BranchPolicy(pfem) {
		AddClassProperty(this, &l2b, "length_to_branch"); AddClassProperty(this, &t2e, "time_to_emerge");
	}
	virtual ~DelayedBranchingPolicyEFD() {}
	//! performs initialization
	bool Init() override {
		return l2b->Init() && t2e->Init() && BranchPolicy::Init();
	}
	//! update to a given time
	void TimeStepUpdate(double current_time) override { BranchPolicy::TimeStepUpdate(current_time); l2b->TimeStepUpdate(current_time); t2e->TimeStepUpdate(current_time); }
	//! do the per element setup
	void SetupBranchInfo(AngioElement * angio_elem) override;
private:
	FEProbabilityDistribution* l2b = nullptr;//length to branch
	FEProbabilityDistribution* t2e = nullptr;//time to emerge
	double discretization_length = 1.0;

	//! Helper class for delayed branching policy
	class BranchPoint
	{
	public:
		explicit BranchPoint(double start_time, double end_time) : _start_time(start_time), _end_time(end_time) { assert(end_time > start_time); }
		//only needed while creating the collection of branch points
		double _start_time = 0.0;
		double _end_time = 0.0;

		double length = 0.0;
		double processed = 0.0;
		int discrete_sections = 0;
		std::vector<Segment *> current_segments;
	};
};