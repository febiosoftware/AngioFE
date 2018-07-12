#include "ContributionMix.h"

#include "angio3d.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FESolidDomain.h>
#include "FEAngioMaterialPoint.h"
#include "FEAngio.h"

double PSCPDDContributionMix::ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh)
{
	return psc_weight*prev;
}

void PSCPDDContributionMix::Update(FEMesh * mesh)
{
	// TODO: Implement.
}

BEGIN_PARAMETER_LIST(PSCPDDContributionMix, ContributionMix)
ADD_PARAMETER(psc_weight, FE_PARAM_DOUBLE, "psc_weight");
END_PARAMETER_LIST();

double ContributionMixManager::ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh)
{
	for (int i = 0; i < cm_modifiers.size(); i++)
	{
		prev = cm_modifiers[i]->ApplyModifiers(prev, angio_element, local_pos, mesh);
	}
	return prev;
}