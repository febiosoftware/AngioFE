#pragma once
#include <FECore/FEMaterial.h>
#include "FEProbabilityDistribution.h"
#include <FEBioMech/FESolidMaterial.h>
#include "BranchingPolicy.h"
#include "FragmentSeeder.h"

//! separate out components that may be shared if an angio material needs to inherit from another material
class CommonAngioProperties :public FEMaterial
{
public:
	CommonAngioProperties(FEModel * pfem);
	~CommonAngioProperties(){}

	//FEPropertyT<GrowDirectionModifiers> gdms;
	FEPropertyT<FragmentSeeder> fseeder;
	FEPropertyT<FESolidMaterial> vessel_material;
};
