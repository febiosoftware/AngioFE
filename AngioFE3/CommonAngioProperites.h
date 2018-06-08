#pragma once
#include <FECore/FEMaterial.h>
#include "FEProbabilityDistribution.h"
#include <FEBioMech/FESolidMaterial.h>
#include "BranchingPolicy.h"
#include "FragmentSeeder.h"


class CommonAngioProperties :public FEMaterial
{
public:
	CommonAngioProperties(FEModel * pfem);
	~CommonAngioProperties(){}

	//FEPropertyT<GrowDirectionModifiers> gdms;
	FEPropertyT<FragmentSeeder> fseeder;
	FEPropertyT<FESolidMaterial> vessel_material;
	

	void UpdateGDMs();
};
