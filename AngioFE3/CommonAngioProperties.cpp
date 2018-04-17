#include "CommonAngioProperites.h"

CommonAngioProperties::CommonAngioProperties(FEModel * pfem) : FEMaterial(pfem)
{
	AddProperty(&vessel_material, "vessel");
	
	AddProperty(&fiber_initializer, "fiber_initializer");
	//AddProperty(&fseeder, "fragment_seeder");
	//AddProperty(&gdms, "grow_direction_modifiers");
}

void CommonAngioProperties::InitializeFibers(FiberManager * man)
{
	fiber_initializer->InitializeFibers(man);
}

void CommonAngioProperties::UpdateGDMs()
{
	//gdms->Update();
}
