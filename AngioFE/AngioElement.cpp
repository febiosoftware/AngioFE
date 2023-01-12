#include "FEAngioMaterial.h"

#ifndef PI
#define PI 3.14159265358979
#endif

double AngioElement::GetLengthAtTime(FEMesh* mesh, double time) const
{
	double len = 0.0;
	for (int i = 0; i < grown_segments.size(); i++)
	{
		len += grown_segments[i]->LengthAtTime(mesh, time);
	}
	return len;
}

void AngioElement::GetNatLengths(	double gr, double gs, double gt, 
									vec3d& er, vec3d& es, vec3d& et) const
{
	double Gr[FESolidElement::MAX_NODES];
	double Gs[FESolidElement::MAX_NODES];
	double Gt[FESolidElement::MAX_NODES];
	auto mesh = _angio_mat->m_pangio->GetMesh();
	_elem->shape_deriv(Gr, Gs, Gt, gr, gs, gt);

	// Calculate the length of each side of the element.
	for (int j = 0; j < _elem->Nodes(); j++)
	{
		er += mesh->Node(_elem->m_node[j]).m_rt * Gr[j];
	}
	for (int j = 0; j < _elem->Nodes(); j++)
	{
		es += mesh->Node(_elem->m_node[j]).m_rt * Gs[j];
	}
	for (int j = 0; j < _elem->Nodes(); j++)
	{
		et += mesh->Node(_elem->m_node[j]).m_rt * Gt[j];
	}
}