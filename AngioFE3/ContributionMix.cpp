#include "ContributionMix.h"

#include "angio3d.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FESolidDomain.h>
#include "FEAngioMaterialPoint.h"
#include "FEAngio.h"

double PSCPDDContributionMix::ApplyModifiers(double dt, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh)
{
		//! SL 8_21
		// doesn't seem to do anything anymore
		//return psc_weight*dt;
		return psc_weight;
}

// this is supposed to update the psc to a load curve? Probably also will allow it to be updated by fractional anisotropy, concentrations, gradients, etc.
void PSCPDDContributionMix::Update(FEMesh * mesh)
{
	// TODO: Implement.
}

BEGIN_FECORE_CLASS(PSCPDDContributionMix, ContributionMix)
ADD_PARAMETER(psc_weight, "psc_weight");
END_FECORE_CLASS();

double ContributionMixManager::ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh)
{
	for (int i = 0; i < cm_modifiers.size(); i++)
	{
		prev = cm_modifiers[i]->ApplyModifiers(prev, angio_element, local_pos, mesh);
	}
	return prev;
}