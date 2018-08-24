#pragma once
#include <FECore/FEMaterial.h>
#include "AngioElement.h"
#include "FEProbabilityDistribution.h"
#include "VariableInterpolation.h"

class FEAngio;

//base class, derived classes store the infomation needed by the brancher on a per element basis
class BranchInfo
{
public:
	virtual ~BranchInfo(){}
};

class ZenithAngle : public FEMaterial
{
public:
	ZenithAngle(FEModel* pfem) : FEMaterial(pfem) {}
	virtual double GetZenithAngle(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element) = 0;
	virtual void TimeStepUpdate(double current_time) = 0;
};

class AzimuthAngle :public FEMaterial
{
public:
	AzimuthAngle(FEModel* pfem) : FEMaterial(pfem) {}
	virtual double GetAzimuthAngle(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element) = 0;
	virtual void TimeStepUpdate(double current_time) = 0;
};

class ZenithAngleProbabilityDistribution : public ZenithAngle
{
public:
	ZenithAngleProbabilityDistribution(FEModel* pfem) : ZenithAngle(pfem) { AddProperty(&angle, "angle"); }
	double GetZenithAngle(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element) override;
	bool Init() override { return angle->Init(); }
	void TimeStepUpdate(double current_time) override { angle->TimeStepUpdate(current_time); }

private:
	FEPropertyT<FEProbabilityDistribution> angle;
};

class AzimuthAngleProbabilityDistribution :public AzimuthAngle
{
public:
	AzimuthAngleProbabilityDistribution(FEModel* pfem) : AzimuthAngle(pfem) { AddProperty(&angle, "angle"); }
	double GetAzimuthAngle(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element) override;
	bool Init() override {return angle->Init();}
	void TimeStepUpdate(double current_time) override { angle->TimeStepUpdate(current_time); }
private:
	FEPropertyT<FEProbabilityDistribution> angle;
};
class BranchPolicy :public FEMaterial
{
public:
	BranchPolicy(FEModel* pfem) : FEMaterial(pfem) { 
		AddProperty(&azimuth_angle, "azimuth_angle"); 
		AddProperty(&zenith_angle, "zenith_angle");
		AddProperty(&interpolation_prop, "interpolation_prop");
	}
	virtual ~BranchPolicy(){};
	virtual void AddBranches(AngioElement * angio_elem, int buffer_index, double end_time, double final_time, double min_scale_factor, FEMesh* mesh, FEAngio* feangio)=0;
	virtual void SetupBranchInfo(AngioElement * angio_elem) = 0;
	bool Init() override { return azimuth_angle->Init() && zenith_angle->Init(); }
	virtual void TimeStepUpdate(double current_time){ azimuth_angle->TimeStepUpdate(current_time); zenith_angle->TimeStepUpdate(current_time);}
	//adds a branch tip at the given natural coordinates
	void AddBranchTip(AngioElement * angio_element, vec3d local_pos, vec3d parent_direction, double start_time, int vessel_id, int buffer_index, FEMesh* mesh);
	vec3d GetBranchDirection(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element, FEMesh* mesh);
private:
	//angles are in radians
	FEPropertyT<AzimuthAngle> azimuth_angle;
	FEPropertyT<ZenithAngle> zenith_angle;
	FEPropertyT<FEVariableInterpolation> interpolation_prop;
};
class FutureBranch
{
public:
	explicit FutureBranch(vec3d local_pos, Segment * parent, double start_time) : _local_pos(local_pos), _parent(parent), _start_time(start_time) { assert(parent); }
	vec3d _local_pos;
	Segment * _parent;
	double _start_time;
};

//info needed by delayed branchin on a per element 
class DelayBranchInfo : public BranchInfo
{
public:
	double length_to_branch = 0.0;
	std::list<FutureBranch> future_branches;
};


//length to branch is calculated only when the last value of length to branch has been hit
class DelayedBranchingPolicy :public BranchPolicy
{
public:
	DelayedBranchingPolicy(FEModel* pfem) : BranchPolicy(pfem) { AddProperty(&l2b, "length_to_branch"); AddProperty(&t2e, "time_to_emerge");
	}
	virtual ~DelayedBranchingPolicy(){}
	bool Init() override {
		return l2b.Init() && t2e.Init() && BranchPolicy::Init(); }
	void TimeStepUpdate(double current_time) override { BranchPolicy::TimeStepUpdate(current_time); l2b->TimeStepUpdate(current_time); t2e->TimeStepUpdate(current_time); }
	void AddBranches(AngioElement * angio_elem, int buffer_index, double end_time, double final_time, double min_scale_factor, FEMesh* mesh, FEAngio* feangio) override;
	void SetupBranchInfo(AngioElement * angio_elem) override;
private:
	FEPropertyT<FEProbabilityDistribution> l2b;//length to branch
	FEPropertyT<FEProbabilityDistribution> t2e;//time to emerge
	double discretization_length = 1.0;

	//helper class for delayed branching policy
	class BranchPoint
	{
	public:
		explicit BranchPoint(double start_time, double end_time) : _start_time(start_time), _end_time(end_time) { assert(end_time > start_time); }
		//only needed while creating the collection of branch points
		double _start_time = 0.0;
		double _end_time = 0.0;

		double length =0.0;
		double processed =0.0;
		int discrete_sections = 0;
		std::vector<Segment *> current_segments;
	};
};