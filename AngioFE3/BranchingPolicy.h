#pragma once
#include <FECore/FEMaterial.h>
#include "AngioElement.h"
#include "FEProbabilityDistribution.h"

class FEAngio;

//base class, derived classes store the infomation needed by the brancher on a per element basis
class BranchInfo
{
public:
	virtual ~BranchInfo(){}
};

class BranchPolicy :public FEMaterial
{
public:
	BranchPolicy(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~BranchPolicy(){};
	virtual void AddBranches(AngioElement * angio_elem, int buffer_index, FEMesh* mesh, FEAngio* feangio)=0;
	virtual void SetupBranchInfo(AngioElement * angio_elem) = 0;
	//adds a branch tip at the given natural coordinates
	void AddBranchTip(AngioElement * angio_element, vec3d local_pos, vec3d parent_direction, double start_time, int vessel_id, int buffer_index);
	vec3d GetBranchDirection(vec3d local_pos, vec3d parent_direction);
};

//info needed by delayed branchin on a per element 
class DelayBranchInfo : public BranchInfo
{
public:
	double length_to_branch;
};


//length to branch is calculated only when the last value of length to branch has been hit
class DelayedBranchingPolicy :public BranchPolicy
{
public:
	DelayedBranchingPolicy(FEModel* pfem) : BranchPolicy(pfem) { AddProperty(&l2b, "length_to_branch"); }
	virtual ~DelayedBranchingPolicy(){}
	bool Init() override {
		return l2b.Init(); }
	void AddBranches(AngioElement * angio_elem, int buffer_index, FEMesh* mesh, FEAngio* feangio) override;
	void SetupBranchInfo(AngioElement * angio_elem) override;
private:
	FEPropertyT<FEProbabilityDistribution> l2b;
	double discretization_length = 1.0;

	//helper class for delayed branching policy
	class BranchPoints
	{
	public:
		//only needed while creating the collection of branch points
		double start_time = 0.0;
		double end_time = 0.0;

		double length =0.0;
		double processed =0.0;
		int discrete_sections = 0;
		std::vector<Segment *> current_segments;
	};
};