#pragma once
#include <FECore/FEMaterial.h>
#include "FEProbabilityDistribution.h"
#include <FEBioMech/FESolidMaterial.h>
#include "FiberManager.h"
#include "FragmentSeeder.h"
#include "GrowDirectionModifier.h"


class CommonAngioProperties :public FEMaterial
{
public:
	CommonAngioProperties(FEModel * pfem);
	~CommonAngioProperties(){}

	//FEPropertyT<GrowDirectionModifiers> gdms;
	//FEPropertyT<FragmentSeeder> fseeder;
	FEPropertyT<FESolidMaterial> vessel_material;
	FEPropertyT<FiberInitializer> fiber_initializer;

	void InitializeFibers(FiberManager * man);
	void UpdateGDMs();
};
