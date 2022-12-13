#include "NodeDataInterpolation.h"
#include "FEAngio.h"
#include <FECore/FENodeDataMap.h>
#include <iostream>

#pragma region FECoreClassDefs
BEGIN_FECORE_CLASS(NodeDataInterpolationManager, FEMaterialProperty)
	ADD_PROPERTY(node_data_interpolation_vals, "node_interpolation_value", FEProperty::Optional);
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(NodeDataInterpolation, FEMaterialProperty)
	ADD_PARAMETER(node_set_id, "node_set_id");
	ADD_PARAMETER(interpolation_mode, "interpolation_mode");
END_FECORE_CLASS()
#pragma endregion FECoreClassDefs

double & DensityValuesNodeDataInterpolation::ValueReference(FEMaterialPoint * mp)
{
	return FEAngioMaterialPoint::FindAngioMaterialPoint(mp)->ref_ecm_density;
}

const char * DensityValuesNodeDataInterpolation::GetDataName() const
{
	return "ref_ecm_density";
}

double & RepulseValuesNodeDataInterpolation::ValueReference(FEMaterialPoint * mp)
{
	return FEAngioMaterialPoint::FindAngioMaterialPoint(mp)->repulse_value;
}

const char * RepulseValuesNodeDataInterpolation::GetDataName() const
{
	return "repulse_value";
}

void NodeDataInterpolation::PrepValues(FEAngio* angio, FEMesh* mesh, FEAngioMaterial* angio_mat)
{
	// convert the node values (lids) to ids

	// for each element in the angio material
	for(int i=0;i < angio->elements_by_material[angio_mat].size();i++)
	{
		//get the solid element
		FESolidElement * se = angio->elements_by_material[angio_mat][i]->_elem;
		// make a vector to store nodal values.
		double nodal_values[FESolidElement::MAX_NODES];
		// store the values at the gauss points
		std::vector<double> values_at_gauss_points;
		// get the values from the solid element gauss points
		for(int j=0; j < se->GaussPoints();j++)
		{
			values_at_gauss_points.push_back(ValueReference(se->GetMaterialPoint(j)));
		}
		// project the solid element values to the nodes from the gauss points
		se->project_to_nodes(&values_at_gauss_points[0], &nodal_values[0]);
		// for each node in the solid elements 
		for(int j=0; j < se->Nodes();j++)
		{
			node_id_to_values[se->m_node[j]] = nodal_values[j];
		}
		
	}
	//replace the values using the nodeset and properties

	// copy the nodeset from the id
	FENodeSet * node_set = mesh->NodeSet(node_set_id);
	assert(node_set);
	// create the data array for the desired type of data
	FEDataArray* da = mesh->FindDataMap(GetDataName());
	assert(da);
	// create a node data map from the data array
	FENodeDataMap * ndm = dynamic_cast<FENodeDataMap*>(da);
	assert(ndm);
	//the number of items is valid
	assert(ndm->DataCount() <= node_set->Size());
	// get the node list
	auto node_list = node_set->GetNodeList();
	// for each node
	for (int i = 0; i<node_list.Size(); i++)
	{
		// get the node id
		int node_id = node_list[i];
		// convert the node id to the values
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
			double nodal_values[FESolidElement::MAX_NODES];
			for(int j=0; j < angio_element->_elem->Nodes();j++)
			{
				nodal_values[j] = node_id_to_values.at(angio_element->_elem->m_node[j]);
			}
			//set gauss values
			for (int j = 0; j < angio_element->_elem->GaussPoints(); j++)
			{
				FEMaterialPoint *mp = angio_element->_elem->GetMaterialPoint(j);
				double H[FESolidElement::MAX_NODES];
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