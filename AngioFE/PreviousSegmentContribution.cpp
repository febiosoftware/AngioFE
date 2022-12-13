#include "PreviousSegmentContribution.h"
#include "Tip.h"
#include "angio3d.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FESolidDomain.h>
#include "FEAngioMaterialPoint.h"
#include "FEAngio.h"

#pragma region FECoreClassDefs
BEGIN_FECORE_CLASS(PreviousSegmentContributionManager, FEMaterialProperty)
	ADD_PROPERTY(psc_modifiers, "psc_modifier", FEProperty::Optional);
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(ProtoPreviousSegmentContributionManager, FEMaterialProperty)
	ADD_PROPERTY(proto_psc_modifiers, "proto_psc_modifier", FEProperty::Optional);
END_FECORE_CLASS()
#pragma endregion FECoreClassDefs

vec3d PreviousSegmentPSC::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh)
{
	return prev_direction;
}

// Managers
vec3d PreviousSegmentContributionManager::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh)
{
	for (int i = 0; i < psc_modifiers.size(); i++)
	{
		prev = psc_modifiers[i]->ApplyModifiers(prev, angio_element, local_pos, prev_direction, mesh);
	}
	return prev;
}

vec3d ProtoPreviousSegmentPSC::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh)
{
	return prev_direction;
}

// Managers
vec3d ProtoPreviousSegmentContributionManager::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh)
{
	for (int i = 0; i < proto_psc_modifiers.size(); i++)
	{
		prev = proto_psc_modifiers[i]->ApplyModifiers(prev, angio_element, local_pos, prev_direction, mesh);
	}
	return prev;
}
