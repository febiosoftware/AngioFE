#pragma once
#include <FECore/FEMaterial.h>
#include "AngioElement.h"
class BranchPolicy :public FEMaterial
{
public:
	virtual ~BranchPolicy() = 0;
	virtual void AddBranches(AngioElement * elem)=0;
};

class DelayedBranchinPolicy :public BranchPolicy
{
public:
	virtual ~DelayedBranchinPolicy(){}
	void AddBranches(AngioElement * elem) override;
};