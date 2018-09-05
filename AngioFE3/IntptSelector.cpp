#include "IntptSelector.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FEMesh.h>
#include "AngioElement.h"
#include "FEAngio.h"

void NarrowIntptSelector::SelectIntPtr(std::vector<FEMaterialPoint*> & int_pts, FESolidElement * se, vec3d natc, FEAngio* feangio, FEMesh* mesh)
{
	//find the search radius for intpts
	double Gr[FEElement::MAX_NODES];
	double Gs[FEElement::MAX_NODES];
	double Gt[FEElement::MAX_NODES];
	se->shape_deriv(Gr, Gs, Gt, natc.x, natc.y, natc.z);

	// Calculate the length of each side of the element.
	vec3d er, es, et; // Basis vectors of the natural coordinates
	for (int j = 0; j < se->Nodes(); j++)
	{
		er += mesh->Node(se->m_node[j]).m_rt * Gr[j];
	}
	for (int j = 0; j < se->Nodes(); j++)
	{
		es += mesh->Node(se->m_node[j]).m_rt * Gs[j];
	}
	for (int j = 0; j < se->Nodes(); j++)
	{
		et += mesh->Node(se->m_node[j]).m_rt * Gt[j];
	}
	double search_radius = std::max(er.norm(), std::max( es.norm(), et.norm()));
	search_radius *= safety_factor;

	vec3d pos(0, 0, 0);
	double arr[FEElement::MAX_NODES];
	se->shape_fnc(arr, natc.x, natc.y, natc.z);
	for (int j = 0; j < se->Nodes(); j++)
	{
		pos += mesh->Node(se->m_node[j]).m_r0* arr[j];
	}


	AngioElement * angio_element = feangio->se_to_angio_element[se];
	std::vector<AngioElement*> & adj_angio_elements = feangio->angio_elements_to_all_adjacent_elements[angio_element];
	for(int i=0; i < se->GaussPoints();i++)
	{
		FEMaterialPoint * mp = se->GetMaterialPoint(i);
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		assert(emp);
		if((emp->m_rt - pos).norm() < search_radius)
		{
			int_pts.push_back(mp);
		}
	}

	for(int i=0; i < adj_angio_elements.size();i++)
	{
		AngioElement * cur_angio_element = adj_angio_elements[i];
		for(int j=0; j < cur_angio_element->_elem->GaussPoints();j++)
		{
			FEMaterialPoint * mp = cur_angio_element->_elem->GetMaterialPoint(j);
			FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
			assert(emp);
			if ((emp->m_rt - pos).norm() < search_radius)
			{
				int_pts.push_back(mp);
			}
		}
	}

}

void AdjacentElementIntptSelector::SelectIntPtr(std::vector<FEMaterialPoint*> & int_pts, FESolidElement * se, vec3d natc, FEAngio* feangio, FEMesh* mesh)
{
	AngioElement * angio_element = feangio->se_to_angio_element[se];
	std::vector<AngioElement*> & adj_angio_elements = feangio->angio_elements_to_all_adjacent_elements[angio_element];
	for (int i = 0; i < se->GaussPoints(); i++)
	{
		FEMaterialPoint * mp = se->GetMaterialPoint(i);
		int_pts.push_back(mp);
	}

	for (int i = 0; i < adj_angio_elements.size(); i++)
	{
		AngioElement * cur_angio_element = adj_angio_elements[i];
		for (int j = 0; j < cur_angio_element->_elem->GaussPoints(); j++)
		{
			FEMaterialPoint * mp = cur_angio_element->_elem->GetMaterialPoint(i);
			int_pts.push_back(mp);
		}
	}
}
