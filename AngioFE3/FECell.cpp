#include "Tip.h"
#include "FECell.h"
#include "Segment.h"
#include <iostream>
#include <algorithm>
#include "FEAngio.h"
#include "CellSpecies.h"
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
	Solutes = other->Solutes;
	SBMs = other->SBMs;
	Species = other->Species;
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
	for (unsigned it = 0; it < SBMs.size(); it++)
	{
		if (SBMs[it] != nullptr) {
			SBMs[it]->SetPosition(GlobalPos);
			SBMs[it]->Update();
		}
	}
	for (unsigned it = 0; it < Solutes.size(); it++)
	{
		if (Solutes[it] != nullptr) {
			Solutes[it]->SetPosition(GlobalPos);
			Solutes[it]->Update();
		}
	}
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

void FECell::InitSBMs(FEMesh* mesh)
{
	FEModel* fem = angio_element->_mat->GetFEModel();
	FEDomain* dom = &mesh->Domain(0);
	FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());
	// assign the SBM properties for each property
	CellSBMManager* m_SBMs = angio_element->_angio_mat->cell_SBM_manager;
	if (m_SBMs) {
		for (int i = 0; i < m_SBMs->cell_SBM_prop.size(); i++)
		{
			// get the sbm id
			int SBMID = m_SBMs->cell_SBM_prop[i]->GetID();
			// get the production rate/concentration
			double prod_rate = m_SBMs->cell_SBM_prop[i]->GetPR();
			// Create new source
			SBMs[SBMID] = new FESBMPointSource(fem);
			// update the id, rate, and position of the new species 
			SBMs[SBMID]->SetSBMID(SBMID);
			SBMs[SBMID]->SetValue(prod_rate);
			SBMs[SBMID]->SetPosition(GetPosition(mesh));
			// initialize and activate the bc
			SBMs[SBMID]->Init();
			SBMs[SBMID]->Activate();
		}
	}
}

void FECell::UpdateSBMs(FEMesh* mesh)
{
	for (unsigned it = 0; it < SBMs.size(); it++)
	{
		SBMs[it]->SetPosition(GetPosition(mesh));
		SBMs[it]->Update();
	}
};

void FECell::InitSolutes(FEMesh* mesh)
{
	FEModel* fem = angio_element->_mat->GetFEModel();
	FEDomain* dom = &mesh->Domain(0);
	FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());
	// assign the SBM properties for each property
	CellSoluteManager* m_sol = angio_element->_angio_mat->cell_Sol_manager;
	if (m_sol) {
		for (int i = 0; i < m_sol->cell_sol_prop.size(); i++)
		{
			// get the sbm id
			int SolID = m_sol->cell_sol_prop[i]->GetID();
			// get the production rate/concentration
			double prod_rate = m_sol->cell_sol_prop[i]->GetPR();
			// Create new source
			Solutes[SolID] = new FESolutePointSource(fem);
			// update the id, rate, and position of the new species 
			Solutes[SolID]->SetSoluteID(SolID);
			Solutes[SolID]->SetRate(prod_rate);
			Solutes[SolID]->SetPosition(GetPosition(mesh));
			// initialize and activate the bc
			Solutes[SolID]->Init();
			Solutes[SolID]->Activate();
		}
	}
}

void FECell::UpdateSolutes(FEMesh* mesh)
{
	for (unsigned it = 0; it < Solutes.size(); it++)
	{
		Solutes[it]->SetPosition(GetPosition(mesh));
		Solutes[it]->Update();
	}
};

void FECell::InitSpecies(FEMesh* mesh)
{
	FEModel* fem = angio_element->_mat->GetFEModel();
	FEDomain* dom = &mesh->Domain(0);
	FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());
	// assign the SBM properties for each property
	CellSpeciesManager* m_species = angio_element->_angio_mat->cell_species_manager;
	if (m_species) {
		for (int i = 0; i < m_species->cell_species_prop.size(); i++)
		{
			// get the sbm id
			int SpeciesID = m_species->cell_species_prop[i]->GetSpeciesID()-1;
			// get the production rate/concentration
			double prod_rate = m_species->cell_species_prop[i]->GetPR();
			// Create new source
			Species[SpeciesID] = new FESolutePointSource(fem);
			// update the id, rate, and position of the new species 
			Species[SpeciesID]->SetSoluteID(SpeciesID);
			Species[SpeciesID]->SetRate(prod_rate);
			Species[SpeciesID]->SetPosition(GetPosition(mesh));
			// initialize and activate the bc
			Species[SpeciesID]->Init();
			Species[SpeciesID]->Activate();
		}
	}
}

void FECell::UpdateSpecies(FEMesh* mesh)
{
	for (unsigned it = 0; it < Species.size(); it++)
	{
		Species[it]->SetPosition(GetPosition(mesh));
		Species[it]->Update();
	}
};
