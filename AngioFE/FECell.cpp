#include "Tip.h"
#include "FECell.h"
#include "Segment.h"
#include <iostream>
#include <algorithm>
#include "FEAngio.h"
#include "CellSpecies.h"
#include "FEProbabilityDistribution.h"
#include <FEBioMix/FEMultiphasicStandard.h>
#include <FEBioMix/FESolute.h>
#include <FECore/log.h>
#include <FECore/FEMaterial.h>

vec3d FECell::GetPosition(FEMesh* mesh) const
{
	double arr[FESolidElement::MAX_NODES];
	assert(angio_element);
	assert(angio_element->_elem);
	angio_element->_elem->shape_fnc(arr, local_pos.x, local_pos.y, local_pos.z);
	vec3d rc(0, 0, 0);

	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		rc += mesh->Node(angio_element->_elem->m_node[j]).m_rt * arr[j];
	}
	return rc;
}

vec3d FECell::GetRefPosition(FEMesh* mesh) const
{
	double arr[FESolidElement::MAX_NODES];
	assert(angio_element);
	assert(angio_element->_elem);
	angio_element->_elem->shape_fnc(arr, local_pos.x, local_pos.y, local_pos.z);
	vec3d rc(0, 0, 0);

	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		rc += mesh->Node(angio_element->_elem->m_node[j]).m_r0 * arr[j];
	}
	return rc;
}

void FECell::PrintCellInfo(FEMesh* mesh, std::string title) const
{
#ifndef NDEBUG
	std::cout << title << std::endl;
	PrintCellInfo(mesh);
#endif
}

void FECell::PrintCellInfo(FEMesh* mesh) const
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
FECell::FECell(FECell* other, FEMesh* mesh)
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

bool FECell::Init() {
	return true;
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
	CellSpeciesManager* m_species = angio_element->_angio_mat->cell_species_manager;
	// assign the properties to each solute. Initialize body loads and add to the solutes container.
	if (m_species) {
		for (int i = 0; i < m_species->cell_solute_prop.size(); i++) {
			CellSolute* cell_solute = new CellSolute(fem);
			CellSolute* ref_solute = m_species->cell_solute_prop[i];
			// get the sbm id
			int SoluteID = ref_solute->GetSoluteID();
			double int_c = ref_solute->GetInt();
			cell_solute->SetInt(int_c);
			// get the production rate/concentration
			// Create new source
			cell_solute->CellSolutePS = new FESolutePointSource(fem);
			cell_solute->SetSoluteID(SoluteID);
			// update the id, rate, and position of the new species 
			cell_solute->CellSolutePS->SetSoluteID(SoluteID);
			cell_solute->CellSolutePS->SetPosition(GetPosition(mesh));
			FESoluteData* psd = cell_solute->FindSoluteData(SoluteID);
			cell_solute->SetDensity(psd->m_rhoT);
			cell_solute->SetMolarMass(psd->m_M);
			cell_solute->SetCharge(psd->m_z);
			cell_solute->CellSolutePS->SetRadius(cell_radius);

			// initialize and activate the bc
			if (cell_solute->CellSolutePS->Init()) {
				cell_solute->CellSolutePS->SetAccumulateFlag(false);
				mesh->GetFEModel()->AddBodyLoad(cell_solute->CellSolutePS);
				Solutes.emplace_back(cell_solute);
			}
		}
		// assign the properties to each SBM. Initialize body loads and add to the SBMs container.
		for (int i = 0; i < m_species->cell_SBM_prop.size(); i++)
		{
			CellSBM* cell_sbm = new CellSBM(fem);
			CellSBM* ref_sbm = m_species->cell_SBM_prop[i];
			// get the sbm id
			int SBMID = ref_sbm->GetSBMID();
			// get the production rate/concentration
			// Create new source
			cell_sbm->CellSBMPS = new FESBMPointSource(fem);
			// update the id, rate, and position of the new species 
			cell_sbm->SetSBMID(SBMID);
			cell_sbm->CellSBMPS->SetSBMID(SBMID);
			cell_sbm->CellSBMPS->SetPosition(GetPosition(mesh));
			FESBMData* psd = cell_sbm->FindSBMData(SBMID);
			cell_sbm->SetDensity(psd->m_rhoT);
			cell_sbm->SetMolarMass(psd->m_M);
			cell_sbm->SetCharge(psd->m_z);
			cell_sbm->CellSBMPS->SetRadius(cell_radius);
			// initialize and activate the bc
			if (cell_sbm->CellSBMPS->Init()) {
				//cell_sbm->CellSBMPS->SetResetFlag(false);
				//cell_sbm->CellSBMPS->SetWeighVolume(true);
				cell_sbm->CellSBMPS->SetAccumulateFlag(false);
				mesh->GetFEModel()->AddBodyLoad(cell_sbm->CellSBMPS);
				SBMs.emplace_back(cell_sbm);
			}
		}
		CellReactionManager* m_reaction = angio_element->_angio_mat->cell_reaction_manager;
		// assign the properties to each solute. Initialize body loads and add to the solutes container.
		if (m_reaction) {
			for (int i = 0; i < m_reaction->cell_reaction.size(); i++) {
				FECellChemicalReaction* m_CR = m_reaction->cell_reaction[i];
				/*FECellChemicalReaction* m_CR = m_reaction->cell_reaction[i];*/
				m_CR->SetCell(this);
				if (m_CR->m_pFwd) {
					m_CR->m_pFwd->Init();
					m_CR->m_pFwd->SetCell(this);
				}
				if (m_CR->m_pRev) {
					m_CR->m_pRev->Init();
					m_CR->m_pRev->SetCell(this);
				}
				m_CR->Init();

				this->Reactions.emplace_back(m_CR);
			}
		}
		FEMaterialPoint& mp = *angio_element->_elem->GetMaterialPoint(0);
	}
}

void FECell::UpdateSpecies(FEMesh* mesh)
{
	for (int isol = 0; isol < Solutes.size(); isol++)
	{
		//! Set the new position
		if (this->eval_time >= 0.0) {
			Solutes[isol]->CellSolutePS->SetPosition(GetPosition(mesh));
		}
		//! Update the body load
		Solutes[isol]->SetSolPRhat(0.0);
	}
	for (int isbm = 0; isbm < SBMs.size(); isbm++)
	{
		//! Set the new position
		if (this->eval_time >= 0.0) {
			SBMs[isbm]->CellSBMPS->SetPosition(this->GetPosition(mesh));
		}
		//! Update the body load
		SBMs[isbm]->SetSBMPRhat(0.0);
	}
	//! total timestep to be taken
	double dt = time - eval_time;
	double t0 = eval_time;
	eval_time = time;
	if (dt == 0) { return; }
	//! Get the maximum timestep before an internal species is negative
	double max_dt = 0.0;
	double t_r = dt - max_dt;
	while (t_r > 0.0)
	{
		max_dt = GetMaxSpeciesDT(mesh, t0, dt);
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
			Sol->SetSolPRhat(0.0);
			Sol->CellSolutePS->SetRate(0.0);
			// combine the molar supplies from all the reactions
			for (int k = 0; k < Reactions.size(); ++k)
			{
				Reactions[k]->SetCell(this);
				if (Reactions[k]->m_pFwd) { Reactions[k]->m_pFwd->SetCell(this); }
				if (Reactions[k]->m_pRev) { Reactions[k]->m_pRev->SetCell(this); }
				//! Get the reaction supply for each reaction that the solute is involved in
				double zetahat;
				zetahat = Reactions[k]->ReactionSupply(this);
				// Get the net stoichiometric ratio for each solute
				double v = Reactions[k]->m_v[isol];
				//Add the product of the stoichiometric ratio with the supply for the reaction
				Sol->AddSolhat(v * zetahat);
				//Sol->AddSolhat(v * zetahat * Sol->GetTolScale());
				std::cout << "tol scale during eval is " << Sol->GetTolScale() << endl;
				// accumulate unless there's no dt. Solutes have the dt evaluated elsewhere.
				Sol->CellSolutePS->Accumulate(Sol->GetSolPRhat() * max_dt / dt);// no dt since the rate is handled elsewhere so this is different than SBMs.
			}
		}

		for (int isbm = 0; isbm < nsbm; isbm++) {
			CellSBM* SBM = SBMs[isbm];
			int SBMID = SBM->GetSBMID();
			SBM->SetSBMhatp(SBM->GetSBMhat());
			SBM->SetSBMhat(0.0);
			SBM->SetSBMPRhat(0.0);
			SBM->CellSBMPS->SetRate(0.0);
			// combine the molar supplies from all the reactions
			for (int k = 0; k < Reactions.size(); ++k)
			{
				Reactions[k]->SetCell(this);
				if (Reactions[k]->m_pFwd) { Reactions[k]->m_pFwd->SetCell(this); }
				if (Reactions[k]->m_pRev) { Reactions[k]->m_pRev->SetCell(this); }
				//! Get the reaction supply for each reaction that the SBM is involved in
				double zetahat;
				zetahat = Reactions[k]->ReactionSupply(this);
				//! Get the net stoichiometric ratio for each sbm
				double v = Reactions[k]->m_v[nsol + isbm];
				//Add the product of the stoichiometric ratio with the supply for the reaction
				SBM->AddSBMhat(v * zetahat);
				SBM->CellSBMPS->Accumulate(SBM->GetSBMPRhat() * max_dt); //need to include dt for SBMs since they are not handled the same as solutes.
			}
		}
		SetInternalSpecies(t0, max_dt);
		t_r = t_r - max_dt;
	}
};

void FECell::SetInternalSpecies(double t0, double dt)
{
	int nsbm = this->SBMs.size();
	int nsol = this->Solutes.size();
	for (int isol = 0; isol < nsol; isol++)
	{
		CellSolute* Sol = Solutes[isol];
		double newc;
		if (t0 != 0.0)
		{
			newc = std::max((Sol->GetInt() + ((Sol->GetSolhat() + Sol->GetSolhatp()) / 2.0) * dt), 0.0);
		}
		else
		{
			newc = std::max((Sol->GetInt() + Sol->GetSolhat() * dt), 0.0);
		}
		Sol->SetInt(newc);
		Sol->CellSolutePS->Update();
	}
	for (int isbm = 0; isbm < nsbm; isbm++)
	{
		CellSBM* SBM = SBMs[isbm];
		double newc;
		if (t0 != 0.0)
		{
			newc = std::max((SBM->GetInt() + ((SBM->GetSBMhat() + SBM->GetSBMhatp()) / 2.0) * dt), 0.0);
		}
		else
		{
			newc = std::max((SBM->GetInt() + SBM->GetSBMhat() * dt), 0.0);
		}
		SBM->SetInt(newc);
		SBM->CellSBMPS->Update();
	}
}

double FECell::GetMaxSpeciesDT(FEMesh* mesh, double t0, double dt)
{
	//! SL: Currently catches the case where cells are constantly producing species. Unsure if this will cause bugs down the line.
	int nsbm = this->SBMs.size();
	int nsol = this->Solutes.size();
	double dt_max = dt;
	for (int isol = 0; isol < nsol; isol++)
	{
		CellSolute* Sol = Solutes[isol];
		double newc;
		double dt_try;
		if (t0 != 0.0)
		{
			newc = Sol->GetInt() + ((Sol->GetSolhat() + Sol->GetSolhatp()) / 2.0) * dt;
			if (newc < 0.0)
			{
				dt_try = std::min((-2.0 * Sol->GetInt()) / (Sol->GetSolhat() + Sol->GetSolhatp()), dt_max);
				if (dt_try > 0.0) { dt_max = dt_try; }
			}
		}
		else
		{
			newc = Sol->GetInt() + Sol->GetSolhat() * dt;
			if (newc < 0.0)
			{
				dt_try = std::min((-2.0 * Sol->GetInt() / Sol->GetSolhat()), dt_max);
				if (dt_try > 0.0) { dt_max = dt_try; }
			}
		}
	}
	for (int isbm = 0; isbm < nsbm; isbm++)
	{
		CellSBM* SBM = SBMs[isbm];
		double newc;
		double dt_try;
		if (t0 != 0.0)
		{
			newc = SBM->GetInt() + ((SBM->GetSBMhat() + SBM->GetSBMhatp()) / 2.0) * dt;
			if (newc < 0.0)
			{
				dt_try = std::min((-2.0 * SBM->GetInt()) / (SBM->GetSBMhat() + SBM->GetSBMhatp()), dt_max);
				if (dt_try > 0.0) { dt_max = dt_try; }
			}
		}
		else
		{
			newc = SBM->GetInt() + SBM->GetSBMhat() * dt;
			if (newc < 0.0)
			{
				dt_try = std::min((-2.0 * SBM->GetInt() / SBM->GetSBMhat()), dt_max);
				if (dt_try > 0.0) { dt_max = dt_try; }
			}
		}
	}
	return dt_max;
}

void FECell::ProtoUpdateSpecies(FEMesh* mesh)
{
	// Update body loads for secreted species
	for (int isol = 0; isol < Solutes.size(); isol++)
	{
		//! Set the new position
		Solutes[isol]->CellSolutePS->SetPosition(GetPosition(mesh));
	}
	for (int isbm = 0; isbm < SBMs.size(); isbm++)
	{
		//! Set the new position
		SBMs[isbm]->CellSBMPS->SetPosition(GetPosition(mesh));
	}
};