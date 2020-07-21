#pragma once
#include <FECore/FEMaterial.h>
#include "FEProbabilityDistribution.h"
#include <FEBioMech/FESolidMaterial.h>
#include "BranchingPolicy.h"
#include "FragmentSeeder.h"
#include <FEBioMech/FEElasticMaterial.h>

//! separate out components that may be shared if an angio material needs to inherit from another material
class CommonAngioProperties :public FEMaterial
{
public:
	//! constructor
	CommonAngioProperties(FEModel * pfem);
	~CommonAngioProperties(){}

	//FEPropertyT<GrowDirectionModifiers> gdms;
	//! fragment seeder
	FragmentSeeder* fseeder;
	//! vessel material
	FEElasticMaterial* vessel_material;
};
