#include "ContributionMix.h"

#include "angio3d.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FESolidDomain.h>
#include "FEAngioMaterialPoint.h"
#include "FEAngio.h"
#include <iostream>
#include <FECore/mathalg.h>

struct EigComp {
	bool operator()(const std::pair<double, vec3d>& x, std::pair<double, vec3d>& y) const {
		if ((x.first) > (y.first)) return true;
		else return false;
	}
};

double PSCPDDContributionMix::ApplyModifiers(double dt, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh)
{
		return psc_weight;
}

// this is supposed to update the psc to a load curve? Probably also will allow it to be updated by fractional anisotropy, concentrations, gradients, etc.
void PSCPDDContributionMix::Update(FEMesh * mesh)
{
	// TODO: Implement.
}

double DensFAContributionMix::ApplyModifiers(double dt, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh)
{
	// vector containing the SPD for each gauss point in the element
	std::vector<mat3ds> SPDs_gausspts;

	// get each gauss point's SPD
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		// get the angio point
		FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		SPDs_gausspts.push_back(angio_mp->angioSPD);
	}

	// array containing the SPD for each node in the element
	mat3ds SPDs_nodes[FESolidElement::MAX_NODES];
	// array for the shape function values
	double H[FESolidElement::MAX_NODES];
	// project the spds from integration points to the nodes
	angio_element->_elem->project_to_nodes(&SPDs_gausspts[0], SPDs_nodes);
	// determine shape function value for the local position
	angio_element->_elem->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
	// Get the interpolated SPD from the shape function-weighted Average Structure Tensor
	mat3ds SPD_int = weightedAverageStructureTensor(SPDs_nodes, H, angio_element->_elem->Nodes());
	// get the vectors of the principal directions and sort in descending order
	// the FA can be calculated as std/rms of the ODF. We can assume each direction is one sample and perform the calculation on the eigenvalues.
	double d[3]; vec3d r[3];
	SPD_int.eigen2(d, r);
	std::vector<std::pair<double, vec3d>> v = { {d[0],r[0]}, {d[1],r[1]}, {d[2],r[2]} };
	// Sort based on absolute value of eigenvalues
	std::sort(v.begin(), v.end(), EigComp());
	double angioFA_int = 1 - (v[1].first / v[0].first);
	double alpha = a0 + a / (1 + exp(b * (angioFA_int-c)));
	return alpha;
}

BEGIN_FECORE_CLASS(DensFAContributionMix, ContributionMix)
ADD_PARAMETER(a0, "a0");
ADD_PARAMETER(a, "a");
ADD_PARAMETER(b, "b");
ADD_PARAMETER(c, "c");
END_FECORE_CLASS();

// this is supposed to update the psc to a load curve? Probably also will allow it to be updated by fractional anisotropy, concentrations, gradients, etc.
void DensFAContributionMix::Update(FEMesh* mesh)
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