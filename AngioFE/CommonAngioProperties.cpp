#include "CommonAngioProperites.h"

CommonAngioProperties::CommonAngioProperties(FEModel * pfem) : FEMaterial(pfem)
{
	AddClassProperty(this, &vessel_material, "vessel");
	//AddProperty(&fiber_initializer, "fiber_initializer");
	AddClassProperty(this, &fseeder, "fragment_seeder");
	//AddProperty(&gdms, "grow_direction_modifiers");
}