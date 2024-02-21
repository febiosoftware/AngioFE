#pragma once
#include <FECore/FEMaterial.h>
#include "AngioElement.h"
#include "FEProbabilityDistribution.h"
#include "VariableInterpolation.h"

class FEAngio;

//! Base class, derived classes store the infomation needed by the brancher on a
//!  per element basis
class BranchInfo
{
public:
	virtual ~BranchInfo(){}
};

//! Implements a class that will return the zentih angle (angle between parent 
//! and new segments)
class ZenithAngle : public FEMaterialProperty
{
	FECORE_BASE_CLASS(ZenithAngle)
public:
	//! Constructor for zenith angle
	ZenithAngle(FEModel* pfem) : FEMaterialProperty(pfem) {}
	//! Gives the zenith angle based on the parameters
	virtual double GetZenithAngle(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element) = 0;
	//! performs any modifications the zenith angle as time changes
	virtual void TimeStepUpdate(double current_time) = 0;
};

//! Implements a class that will return the azimuth angle (rotation about the 
//! parent segment)
class AzimuthAngle :public FEMaterialProperty
{
	FECORE_BASE_CLASS(AzimuthAngle)
public:
	//! Constructor for azimuth angle
	AzimuthAngle(FEModel* pfem) : FEMaterialProperty(pfem) {}
	//! Gives the zenith angle based on the parameters
	virtual double GetAzimuthAngle(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element) = 0;
	//! performs any modifications the zenith angle as time changes
	virtual void TimeStepUpdate(double current_time) = 0;
};

//! Implements a probability distribution to determine the zenith angle (angle 
//! between parent and new segments)
class ZenithAngleProbabilityDistribution : public ZenithAngle
{
public:
	//! constructor for class
	ZenithAngleProbabilityDistribution(FEModel* pfem) : ZenithAngle(pfem) { }
	//! returns the zenith angle for a branch based on a probability distribution
	double GetZenithAngle(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element) override;
	//! performs initialization
	bool Init() override { return angle->Init(); }
	//! Updates the zenith angle to the given time step (may ajust probabilities 
	//! based on time)
	void TimeStepUpdate(double current_time) override 
	{ 
		angle->TimeStepUpdate(current_time); 
	}
protected:
	DECLARE_FECORE_CLASS()
private:
	FEProbabilityDistribution* angle = nullptr;
};

//! Implements a probability distribution to determine the azimuth angle (rotation about the parent segment)
class AzimuthAngleProbabilityDistribution :public AzimuthAngle
{
public:
	//! constructor for class
	AzimuthAngleProbabilityDistribution(FEModel* pfem) : AzimuthAngle(pfem) { }
	//! Returns the azimuth angle based on a probability distribution
	double GetAzimuthAngle(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element) override;
	//! performs initialization
	bool Init() override {return angle->Init();}
	//! Updates the zenith angle to the given time step(may ajust probabilities based on time)
	void TimeStepUpdate(double current_time) override 
	{ 
		angle->TimeStepUpdate(current_time); 
	}
protected:
	DECLARE_FECORE_CLASS()
private:
	FEProbabilityDistribution* angle = nullptr;
};

//! Branch Policy class that determines where and when branches occur
class BranchPolicy :public FEMaterialProperty
{
	FECORE_BASE_CLASS(BranchPolicy)

public:
	//! constructor for class
	BranchPolicy(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~BranchPolicy(){};
	//! adds the branches for a given element
	//! setup any data structures that are needed on a per element basis
	virtual void SetupBranchInfo(AngioElement * angio_elem) = 0;
	//! performs initialization
	bool Init() override 
	{ 
		return azimuth_angle->Init() && zenith_angle->Init(); 
	}
	//! Updates the branch policy to the given time step
	virtual void TimeStepUpdate(double current_time){ azimuth_angle->TimeStepUpdate(current_time); zenith_angle->TimeStepUpdate(current_time);}
	//! Return the direction of a branch
	void AddBranchTipEFD(AngioElement * angio_element, vec3d local_pos, Segment* parent_seg, double start_time, int vessel_id, int buffer_index, FEMesh* mesh);
	//! Return the direction of a branch
	vec3d GetBranchDirectionEFD(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element, FEMesh* mesh);
protected:
	DECLARE_FECORE_CLASS()
private:
	//angles are in radians
	AzimuthAngle* azimuth_angle = nullptr;
	ZenithAngle* zenith_angle = nullptr;
	FEVariableInterpolation* interpolation_prop = nullptr;
};

//! The data needed for a future branch. These are branches whose position and orientation are determined but have not started growing.
class FutureBranch
{
public:
	//! constructor for class
	explicit FutureBranch(vec3d local_pos, Segment * parent, double start_time) : _local_pos(local_pos), _parent(parent), _start_time(start_time) 
	{ 
		assert(parent); 
	}
	//! return the local position
	vec3d _local_pos;
	//! the segment this grew from
	Segment * _parent = nullptr;
	//! the time at which the branch starts to grow
	double _start_time;
};

//! Info needed by delayed branching on a per element basis
class DelayBranchInfo : public BranchInfo
{
public:
	//! the length before a branch occurs
	double length_to_branch = std::numeric_limits<double>::max();
	//! all branches which may occur in the future
	std::list<FutureBranch> future_branches;
};

//! As growth occurs length to branch is calculated, once length to branch is hit a time to emerge is sampled from a probability distribution.Length to branch is calculated only when the last value of length to branch has been hit
class DelayedBranchingPolicyEFD :public BranchPolicy
{
	FECORE_BASE_CLASS(DelayedBranchingPolicyEFD)
public:
	//! constructor for class
	DelayedBranchingPolicyEFD(FEModel* pfem) : BranchPolicy(pfem) {}
	virtual ~DelayedBranchingPolicyEFD() {}
	//! performs initialization
	bool Init() override 
	{
		return l2b->Init() && t2e->Init() && BranchPolicy::Init();
	}
	//! update to a given time
	void TimeStepUpdate(double current_time) override 
	{ 
		BranchPolicy::TimeStepUpdate(current_time); 
		l2b->TimeStepUpdate(current_time); 
		t2e->TimeStepUpdate(current_time); 
	}
	//! do the per element setup
	void SetupBranchInfo(AngioElement * angio_elem) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	FEProbabilityDistribution* l2b = nullptr;//length to branch (distance between branches along vessel)
	FEProbabilityDistribution* t2e = nullptr;//time to emerge (time since the last branch occurred on this vessel)
	double discretization_length = 1.0;

	//! Helper class for delayed branching policy
	class BranchPoint
	{
	public:
		explicit BranchPoint(double start_time, double end_time) : _start_time(start_time), _end_time(end_time) 
		{ 
			assert(end_time > start_time); 
		}
		//only needed while creating the collection of branch points
		double _start_time = 0.0;
		double _end_time = 0.0;
		double length = 0.0;
		double processed = 0.0;
		int discrete_sections = 0;
		std::vector<Segment *> current_segments;
	};
};