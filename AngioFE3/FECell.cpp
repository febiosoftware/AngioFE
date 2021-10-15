#include "Tip.h"
#include "FECell.h"
#include "Segment.h"
#include <iostream>
#include <algorithm>
#include "FEAngio.h"
#include "CellSpecies.h"
#include "FEProbabilityDistribution.h"
#include <FEBioMix\FEMultiphasicStandard.h>

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
}

vec3d FECell::GetLocalPosition() const
{
	return local_pos;
}

void FECell::InitSpecies(FEMesh* mesh)
{
	FEModel* fem = angio_element->_mat->GetFEModel();
	FEDomain* dom = &mesh->Domain(0);
	FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());
	CellSpeciesManager* m_species = angio_element->_angio_mat->cell_species_manager;
	// assign the properties to each solute. Initialize body loads and add to the solutes container.
	if (m_species) {
		for (int i = 0; i < m_species->cell_solute_prop.size(); i++) {
			CellSolute* cell_solute = new CellSolute(fem);
			//CellSolute* cell_solute = m_species->cell_solute_prop[i];
			//cell_solute = m_species->cell_solute_prop[i];
			CellSolute* ref_solute = m_species->cell_solute_prop[i];
			// get the sbm id
			int SoluteID = ref_solute->GetSoluteID();
			// get the production rate/concentration
			double prod_rate = ref_solute->GetPR();
			// Create new source
			cell_solute->CellSolutePS = new FESolutePointSource(fem);
			cell_solute->SetSoluteID(SoluteID);
			// update the id, rate, and position of the new species 
			cell_solute->CellSolutePS->SetSoluteID(SoluteID);
			cell_solute->CellSolutePS->SetRate(prod_rate);
			cell_solute->CellSolutePS->SetPosition(GetPosition(mesh));
			// initialize and activate the bc
			if (cell_solute->CellSolutePS->Init()) {
				mesh->GetFEModel()->AddBodyLoad(cell_solute->CellSolutePS);
				Solutes.emplace_back(cell_solute);
			}
		}
		// assign the properties to each SBM. Initialize body loads and add to the SBMs container.
		for (int i = 0; i < m_species->cell_SBM_prop.size(); i++)
		{
			CellSBM* cell_sbm = new CellSBM(fem);
			CellSBM* ref_sbm = m_species->cell_SBM_prop[i];
			//CellSBM* cell_sbm = m_species->cell_SBM_prop[i];
			// get the sbm id
			int SBMID = ref_sbm->GetSBMID();
			// get the production rate/concentration
			double prod_rate = ref_sbm->GetPR();
			// Create new source
			cell_sbm->CellSBMPS = new FESBMPointSource(fem);
			// update the id, rate, and position of the new species 
			cell_sbm->CellSBMPS->SetSBMID(SBMID);
			cell_sbm->CellSBMPS->SetValue(prod_rate);
			cell_sbm->CellSBMPS->SetPosition(GetPosition(mesh));
			// initialize and activate the bc
			if (cell_sbm->CellSBMPS->Init()) {
				mesh->GetFEModel()->AddBodyLoad(cell_sbm->CellSBMPS);
				SBMs.emplace_back(cell_sbm);
			}
		}
		CellReactionManager* m_reaction = angio_element->_angio_mat->cell_reaction_manager;
		// assign the properties to each solute. Initialize body loads and add to the solutes container.
		if (m_reaction) {
			for (int i = 0; i < m_reaction->cell_reaction.size(); i++) {
				FECellReaction* m_R = m_reaction->cell_reaction[i];
				FECellChemicalReaction* m_CR = m_reaction->cell_reaction[i];

				//m_CR->m_pMP = mat;
				////mat->AddChemicalReaction(m_CR);
				//if (m_CR->Init()) {
				//	Reactions.emplace_back(m_CR);
				//}
			}
		}
		FEMaterialPoint& mp = *angio_element->_elem->GetMaterialPoint(0);
	}
}

void FECell::UpdateSpecies(FEMesh* mesh)
{
	// Update body loads for secreted species
	for (int isol = 0; isol < Solutes.size(); isol++)
	{
		//! Set the new position
		Solutes[isol]->CellSolutePS->SetPosition(GetPosition(mesh));
		//! Call the Body Load update
		Solutes[isol]->CellSolutePS->Update();
	}
	for (int isbm = 0; isbm < SBMs.size(); isbm++)
	{
		//! Set the new position
		SBMs[isbm]->CellSBMPS->SetPosition(GetPosition(mesh));
		//! Call the Body Load update
		SBMs[isbm]->CellSBMPS->Update();
	}

	// update SBMs
	double ctime = mesh->GetFEModel()->GetTime().currentTime;
	double dt = ctime - eval_time;
	eval_time = ctime;
	//double dt = mesh->GetFEModel()->GetTime().timeIncrement;
	int nsbm = this->SBMs.size();
	int nsol = this->Solutes.size();
	FEModel* fem = angio_element->_mat->GetFEModel();
	FEDomain* dom = &mesh->Domain(0);
	FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());

	FEMaterialPoint& mp = *angio_element->_elem->GetMaterialPoint(0);
	//! Solutes
		
		for (int isol = 0; isol < nsol; isol++) {
			CellSolute* Sol = Solutes[isol];
			int SolID = Sol->GetSoluteID();
			Sol->SetSolhatp(Sol->GetSolhat());
			Sol->SetSolhat(0.0);
			// combine the molar supplies from all the reactions
			for (int k = 0; k < Reactions.size(); ++k) {
				//! Get the reaction supply for each reaction that the solute is involved in
				double zetahat = GetReactionSupply(Reactions[k]);
				//std::cout << "zetahat is " << zetahat << endl;
				// Get the net stoichiometric ratio for each solute
				double v = Reactions[k]->m_v[isol];
				//Add the product of the stoichiometric ratio with the supply for the reaction
				Sol->AddSolhat(v*zetahat);
			}
			// perform the time integration (midpoint rule)
			double newc = Sol->GetInt() + (Sol->GetSolhat() + Sol->GetSolhatp()) / 2 * dt;
			Sol->SetInt(std::max(newc, 0.0));
		}

		for (int isbm = 0; isbm < nsbm; isbm++) {
			CellSBM* SBM = SBMs[isbm];
			int SolID = SBM->GetSBMID();
			SBM->SetSBMhatp(SBM->GetSBMhat());
			SBM->SetSBMhat(0.0);
			// combine the molar supplies from all the reactions
			for (int k = 0; k < Reactions.size(); ++k) {
				double zetahat = GetReactionSupply(Reactions[k]);
				//! Get the net stoichiometric ratio for each sbm
				double v = Reactions[k]->m_v[nsol+isbm];
				SBM->AddSBMhat(v*zetahat);
			}
			// perform the time integration (midpoint rule)
			double newc = SBM->GetInt() + (SBM->GetSBMhat() + SBM->GetSBMhatp()) / 2 * dt;
			SBM->SetInt(std::max(newc, 0.0));
		}
};

double FECell::GetReactionSupply(FEChemicalReaction* m_R) {
	//! Determine the reaction supply for the desired reaction
	FEMaterialPoint& mp = *angio_element->_elem->GetMaterialPoint(0);
	double kfwd = m_R->m_pFwd->ReactionRate(mp);
	double zhat = kfwd;
	// get the produduct of the reactants exponentiated by their reactant stoichiometric coefficient
	//! Final form is zhat = kF*Pi_a(pow(c_a,v_a,R))
	int nsol = Solutes.size();
	for (int i = 0; i < nsol; ++i) {
		int vR = m_R->m_vR[i];
		if (vR > 0) {
			double c = Solutes[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}
	// add contribution of solid-bound molecules
	int nsbm = SBMs.size();
	for (int i = 0; i < nsbm; ++i) {
		//! m_vR appends sbms to solutes
		int vR = m_R->m_vR[nsol + i];
		if (vR > 0) {
			double c = SBMs[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}
	return zhat;
}