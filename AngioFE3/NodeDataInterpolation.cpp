#include "NodeDataInterpolation.h"
#include "FEAngio.h"
#include <FECore/FENodeDataMap.h>

BEGIN_PARAMETER_LIST(NodeDataInterpolation, FEMaterial)
ADD_PARAMETER(node_set_id, FE_PARAM_INT, "node_set_id");
ADD_PARAMETER(interpolation_mode, FE_PARAM_INT, "interpolation_mode");
END_PARAMETER_LIST();
double & DensityValuesNodeDataInterpolation::ValueReference(FEMaterialPoint * mp)
{
	return FEAngioMaterialPoint::FindAngioMaterialPoint(mp)->ref_ecm_density;
}

const char * DensityValuesNodeDataInterpolation::GetDataName() const
{
	return "ref_ecm_density";
}

void NodeDataInterpolation::PrepValues(FEAngio* angio, FEMesh* mesh, FEAngioMaterial* angio_mat)
{
	for(int i=0;i < angio->elements_by_material[angio_mat].size();i++)
	{
		FESolidElement * se = angio->elements_by_material[angio_mat][i]->_elem;
		double nodal_values[FEElement::MAX_NODES];
		std::vector<double> values_at_gauss_points;
		for(int j=0; j < se->GaussPoints();j++)
		{
			values_at_gauss_points.push_back(ValueReference(se->GetMaterialPoint(j)));
		}
		se->project_to_nodes(&values_at_gauss_points[0], nodal_values);
		for(int j=0; j < se->Nodes();j++)
		{
			node_id_to_values[se->m_node[j]] = nodal_values[j];
		}
		
	}
	//replace the values using the nodeset and properties

	FENodeSet * node_set = mesh->NodeSet(node_set_id);
	assert(node_set);
	FEDataArray* da = angio->m_fem->FindDataArray(GetDataName());
	assert(da);
	FENodeDataMap * ndm = dynamic_cast<FENodeDataMap*>(da);
	assert(ndm);
	//the number of items is valid
	assert(ndm->DataCount() <= node_set->size());
	auto node_list = node_set->GetNodeList();
	//not working
	//for(int i=0; i<ndm->DataCount();i++)
	for (int i = 0; i<node_list.size(); i++)
	{
		int node_id = node_list[i];
		node_id_to_values[node_id] = ndm->getValue(i);
	}

	//now map the values back to the gauss points

	for(int i=0; i < angio->elements_by_material[angio_mat].size();i++)
	{
		AngioElement * angio_element = angio->elements_by_material[angio_mat][i];
		switch (interpolation_mode)
		{
		case 0:
			//get nodal values
			double nodal_values[FEElement::MAX_NODES];
			for(int j=0; j < angio_element->_elem->Nodes();j++)
			{
				nodal_values[j] = node_id_to_values.at(angio_element->_elem->m_node[j]);
			}
			//set gauss values
			for (int j = 0; j < angio_element->_elem->GaussPoints(); j++)
			{
				FEMaterialPoint *mp = angio_element->_elem->GetMaterialPoint(j);
				double H[FEElement::MAX_NODES];
				angio_element->_elem->shape_fnc(H, angio_element->_elem->gr(j), angio_element->_elem->gs(j), angio_element->_elem->gt(j));
				double val = 0;
				for(int k=0; k < angio_element->_elem->Nodes();k++)
				{
					val += H[k]*nodal_values[k];
				}
				ValueReference(mp) = val;
			}
			break;
		case 1:
		default:
			assert(false);
		}
	}
}


void NodeDataInterpolationManager::DoInterpolations(FEAngio * angio, FEMesh* mesh, FEAngioMaterial* angio_mat)
{
	for(int i=0; i < node_data_interpolation_vals.size();i++)
	{
		node_data_interpolation_vals[i]->PrepValues(angio, mesh, angio_mat);
	}
}