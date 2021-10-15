#include "ContributionMix.h"

#include "angio3d.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FESolidDomain.h>
#include "FEAngioMaterialPoint.h"
#include "FEAngio.h"
#include <iostream>

double PSCPDDContributionMix::ApplyModifiers(double dt, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh)
{
		return psc_weight*dt;  
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

double ProtoPSCPDDContributionMix::ApplyModifiers(double dt, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh)
{
	// dt is assumed to be 1 for the proto step as it only runs from -1 to 0. Initial segments will only grow straight if the dt is used.
	return proto_psc_weight;
}

void ProtoPSCPDDContributionMix::Update(FEMesh * mesh)
{
	// TODO: Implement.
}

BEGIN_FECORE_CLASS(ProtoPSCPDDContributionMix, ProtoContributionMix)
ADD_PARAMETER(proto_psc_weight, "proto_psc_weight");
END_FECORE_CLASS();

double ProtoContributionMixManager::ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh)
{
	for (int i = 0; i < proto_cm_modifiers.size(); i++)
	{
		prev = proto_cm_modifiers[i]->ApplyModifiers(prev, angio_element, local_pos, mesh);
	}
	return prev;
}