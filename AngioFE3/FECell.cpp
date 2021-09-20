#include "Tip.h"
#include "FECell.h"
#include "Segment.h"
#include <iostream>
#include <algorithm>
#include "FEAngio.h"
#include "TipSpecies.h"
#include "FEProbabilityDistribution.h"

vec3d FECell::GetPosition(FEMesh * mesh) const
{
	double arr[FESolidElement::MAX_NODES];
	assert(angio_element);
	assert(angio_element->_elem);
	angio_element->_elem->shape_fnc(arr, local_pos.x, local_pos.y, local_pos.z);
	vec3d rc(0, 0, 0);

	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		rc += mesh->Node(angio_element->_elem->m_node[j]).m_rt* arr[j];
	}
	return rc;
}

vec3d FECell::GetRefPosition(FEMesh * mesh) const
{
	double arr[FESolidElement::MAX_NODES];
	assert(angio_element);
	assert(angio_element->_elem);
	angio_element->_elem->shape_fnc(arr, local_pos.x, local_pos.y, local_pos.z);
	vec3d rc(0, 0, 0);

	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		rc += mesh->Node(angio_element->_elem->m_node[j]).m_r0* arr[j];
	}
	return rc;
}

void FECell::PrintCellInfo(FEMesh *mesh, std::string title) const
{
#ifndef NDEBUG
	std::cout << title << std::endl;
	PrintCellInfo(mesh);
#endif
}

void FECell::PrintCellInfo(FEMesh *mesh) const
{
#ifndef NDEBUG
	/*std::cout << "local position: " << local_pos.x << " , " << local_pos.y << " , " << local_pos.z << std::endl;
	vec3d pos = GetPosition(mesh);
	std::cout << "global position: " << pos.x << " , " << pos.y << " , " << pos.z << std::endl;
	std::cout << "element id: " << angio_element->_elem->GetID() << std::endl;*/
	//if (face)
	//{
	//	std::cout << "face id: " << face->_elem->GetID() << std::endl;
	//}
	//else
	//{
	//	std::cout << "no next face" << std::endl;
	//}

	//vec3d dir = GetDirection(mesh);
	/*std::cout << "direction: " << dir.x << " , " << dir.y << " , " << dir.z << std::endl;

	std::cout << "is branch: " << is_branch << std::endl;
	std::cout << "use direction: " << use_direction << std::endl;
	std::cout << std::endl;*/
#endif
}

// create a new cell based on the pre-existing cell.
FECell::FECell(FECell * other, FEMesh * mesh)
{
	// get the other tip's parameters.
	angio_element = other->angio_element;
	local_pos = other->local_pos;
	// deparent the new tip
	initial_cell_id = other->initial_cell_id;
	time = other->time;
	FEModel* fem = angio_element->_mat->GetFEModel();
	Species = other->Species;
	// inherit species in the new tip and remove them in the parent.
	//Species = other->Species;
	//other->Species.clear();
}

void FECell::SetLocalPosition(vec3d pos, FEMesh* mesh)
{
	assert(angio_element);
	assert(angio_element->_elem->Type() == FE_Element_Type::FE_TET4G4 ? (pos.x + pos.y + pos.z) < 1.01 : true);
	assert(angio_element->_elem->Type() == FE_Element_Type::FE_TET4G4 ? (pos.x < 1.01) && (pos.x >= -0.01) : true);
	assert(angio_element->_elem->Type() == FE_Element_Type::FE_TET4G4 ? (pos.y < 1.01) && (pos.y >= -0.01) : true);
	assert(angio_element->_elem->Type() == FE_Element_Type::FE_TET4G4 ? (pos.z < 1.01) && (pos.z >= -0.01) : true);

	assert(angio_element->_elem->Type() == FE_Element_Type::FE_TET4G1 ? (pos.x + pos.y + pos.z) < 1.01 : true);
	assert(angio_element->_elem->Type() == FE_Element_Type::FE_TET4G1 ? (pos.x < 1.01) && (pos.x >= -0.01) : true);
	assert(angio_element->_elem->Type() == FE_Element_Type::FE_TET4G1 ? (pos.y < 1.01) && (pos.y >= -0.01) : true);
	assert(angio_element->_elem->Type() == FE_Element_Type::FE_TET4G1 ? (pos.z < 1.01) && (pos.z >= -0.01) : true);

	local_pos = pos;
	local_pos.x = std::max(std::min(FEAngio::NaturalCoordinatesUpperBound_r(angio_element->_elem->Type()), local_pos.x), FEAngio::NaturalCoordinatesLowerBound_r(angio_element->_elem->Type()));
	local_pos.y = std::max(std::min(FEAngio::NaturalCoordinatesUpperBound_s(angio_element->_elem->Type()), local_pos.y), FEAngio::NaturalCoordinatesLowerBound_s(angio_element->_elem->Type()));
	local_pos.z = std::max(std::min(FEAngio::NaturalCoordinatesUpperBound_t(angio_element->_elem->Type()), local_pos.z), FEAngio::NaturalCoordinatesLowerBound_t(angio_element->_elem->Type()));
	vec3d GlobalPos = GetPosition(mesh);
	for (unsigned it = 0; it < Species.size(); it++)
	{
		if (Species[it] != nullptr) {
			Species[it]->SetPosition(GlobalPos);
			Species[it]->Update();
		}
	}
}

vec3d FECell::GetLocalPosition() const
{
	return local_pos;
}

void FECell::InitSBM(FEMesh* mesh)
{
	FEModel* fem = angio_element->_mat->GetFEModel();
	FEDomain* dom = &mesh->Domain(0);
	FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());
	// assign the SBM properties for each property
	TipSpeciesManager* m_species = angio_element->_angio_mat->tip_species_manager;
	if (m_species) {
		for (int i = 0; i < m_species->tip_species_prop.size(); i++)
		{
			// get the sbm id
			int SBMID = m_species->tip_species_prop[i]->GetID();
			// get the production rate/concentration
			double prod_rate = m_species->tip_species_prop[i]->GetPR();
			// Create new source
			Species[SBMID] = new FESBMPointSource(fem);
			// update the id, rate, and position of the new species 
			Species[SBMID]->SetSBM(SBMID, prod_rate);
			Species[SBMID]->SetPosition(GetPosition(mesh));
			// initialize and activate the bc
			Species[SBMID]->Init();
			Species[SBMID]->Activate();
		}
	}
}

void FECell::UpdateSBM(FEMesh* mesh)
{
	for (unsigned it = 0; it < Species.size(); it++)
	{
		Species[it]->SetPosition(GetPosition(mesh));
		Species[it]->Update();
	}
};