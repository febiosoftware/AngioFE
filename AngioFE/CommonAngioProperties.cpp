#include "CommonAngioProperites.h"

#pragma region FECoreClassDefs
BEGIN_FECORE_CLASS(CommonAngioProperties, FEMaterialProperty)
	ADD_PROPERTY(vessel_material, "vessel");
	ADD_PROPERTY(fseeder, "fragment_seeder");
END_FECORE_CLASS()
#pragma endregion FECoreClassDefs

		CommonAngioProperties::CommonAngioProperties(FEModel* pfem) : FEMaterialProperty(pfem)
{

}