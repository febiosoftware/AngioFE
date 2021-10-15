#include "CellSpecies.h"
#include "FEAngio.h"
#include "FEAngio.h"
#include "FEAngioMaterial.h"
#include <iostream>
#include <FECore/log.h>

//BEGIN_FECORE_CLASS(CellSpecies, FEMaterial)
//ADD_PARAMETER(species_ID, "species_ID");
//ADD_PARAMETER(prod_rate, "prod_rate");
//ADD_PARAMETER(n_species, "n_species");
//END_FECORE_CLASS();

BEGIN_FECORE_CLASS(CellSBM, FEMaterial)
ADD_PARAMETER(SBM_ID, "SBM_ID");
ADD_PARAMETER(SBM_prod_rate, "SBM_prod_rate");
ADD_PARAMETER(n_SBM, "n_SBM");
END_FECORE_CLASS();

BEGIN_FECORE_CLASS(CellSolute, FEMaterial)
ADD_PARAMETER(Solute_ID, "Solute_ID");
ADD_PARAMETER(Solute_prod_rate, "Solute_prod_rate");
ADD_PARAMETER(n_Solute, "n_Solute");
END_FECORE_CLASS();

BEGIN_FECORE_CLASS(FECellReactionRateConst, FECellReactionRate)
ADD_PARAMETER(m_k, FE_RANGE_GREATER_OR_EQUAL(0.0), "k");
END_FECORE_CLASS();

FECellReaction::FECellReaction(FEModel* pfem) : FEMaterial(pfem) 
{
	m_cell = 0;
}

//! SL: What should this do?
bool FECellReaction::Init() {
	return true;
}

BEGIN_FECORE_CLASS(FECellChemicalReaction, FECellReaction)
ADD_PARAMETER(m_Vbar, "Vbar");
ADD_PARAMETER(m_vRtmp, "vR");
ADD_PARAMETER(m_vPtmp, "vP");

// set material properties
ADD_PROPERTY(m_pFwd, "forward_rate", FEProperty::Optional);
ADD_PROPERTY(m_pRev, "reverse_rate", FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FECellChemicalReaction::FECellChemicalReaction(FEModel* pfem) : FECellReaction(pfem)
{
	// additional initializations
	m_Vovr = false;

	m_pFwd = m_pRev = 0;
}

//-----------------------------------------------------------------------------
bool FECellChemicalReaction::Init()
{
	// initialize base class
	FECellReaction::Init();

	// set the parents for the reaction rates
	if (m_pFwd) m_pFwd->m_pReact = this;
	if (m_pRev) m_pRev->m_pReact = this;

	// initialize the reaction coefficients
	int isol, isbm, itot;

	int nsol, nsbm, ntot;
	//nsol = m_pMP->Solutes();
	//nsbm = m_pMP->SBMs();
	nsol = m_cell->Solutes.size();
	nsbm = m_cell->SBMs.size();
	ntot = nsol + nsbm;

		// initialize the stoichiometric coefficients to zero
	m_nsol = nsol;
	m_vR.assign(ntot, 0);
	m_vP.assign(ntot, 0);
	m_v.assign(ntot, 0);

	// cycle through all the solutes in the mixture and determine
	// if they participate in this reaction
	itrmap it;
	intmap solR = m_solR;
	intmap solP = m_solP;
	for (isol = 0; isol<nsol; ++isol) {
		int sid = isol;
		/*sid = m_pMP->GetSolute(isol)->GetSoluteID() - 1;*/
		sid = m_cell->Solutes[isol]->GetSoluteID() - 1;
		it = solR.find(sid);
		if (it != solR.end()) m_vR[isol] = it->second;
		it = solP.find(sid);
		if (it != solP.end()) m_vP[isol] = it->second;
	}

	// cycle through all the solid-bound molecules in the mixture
	// and determine if they participate in this reaction
	intmap sbmR = m_sbmR;
	intmap sbmP = m_sbmP;
	for (isbm = 0; isbm<nsbm; ++isbm) {
		/*int sid = m_pMP->GetSBM(isbm)->GetSBMID() - 1;*/
		int sid = m_cell->SBMs[isbm]->GetSBMID() - 1;
		it = sbmR.find(sid);
		if (it != sbmR.end()) m_vR[nsol + isbm] = it->second;
		it = sbmP.find(sid);
		if (it != sbmP.end()) m_vP[nsol + isbm] = it->second;
	}

	// evaluate the net stoichiometric coefficient
	for (itot = 0; itot<ntot; ++itot) {
		m_v[itot] = m_vP[itot] - m_vR[itot];
	}

	// evaluate the weighted molar volume of reactants and products
	//if (!m_Vovr) {
	//	m_Vbar = 0;
	//	/*m_pMP = cell->angio_elem-*/
	//	for (isol = 0; isol<nsol; ++isol){
	//		m_Vbar += m_v[isol] * m_cell->Solutes[isol]->MolarMass() / cell->Solutes[isol]->Density();
	//		/*m_Vbar += m_v[isol] * m_pMP->GetSolute(isol)->MolarMass() / m_pMP->GetSolute(isol)->Density();*/
	//	}
	//	for (isbm = 0; isbm<nsbm; ++isbm)
	//		/*m_Vbar += m_v[nsol + isbm] * m_pMP->GetSBM(isbm)->MolarMass() / m_pMP->GetSBM(isbm)->Density();*/
	//		m_Vbar += m_v[nsol + isbm] * m_cell->SBMs[isbm]->MolarMass() / cell->SBMs[isbm]->Density();
	//}

	// check that the chemical reaction satisfies electroneutrality
	//int znet = 0;
	//for (isol = 0; isol<nsol; ++isol){
	//		/*znet += m_v[isol] * m_pMP->GetSolute(isol)->ChargeNumber();*/
	//		znet += m_v[isol] * m_cell->Solutes[isol]->ChargeNumber();
	//}
	//for (isbm = 0; isbm<nsbm; ++isbm)
	//	/*znet += m_v[nsol + isbm] * m_pMP->GetSBM(isbm)->ChargeNumber();*/
	//	znet += m_v[nsol + isbm] * m_cell->SBMs[isbm]->ChargeNumber();
	//if (znet != 0) {
	//	feLogError("chemical reaction must satisfy electroneutrality");
	//	return false;
	//}

	return true;
}

//-----------------------------------------------------------------------------
void FECellChemicalReaction::SetParameter(FEParam& p)
{
	if (strcmp(p.name(), "Vbar") == 0)
	{
		m_Vovr = true;
	}
}

//-----------------------------------------------------------------------------
bool FECellChemicalReaction::SetParameterAttribute(FEParam& p, const char* szatt, const char* szval)
{
	// get number of DOFS
	DOFS& fedofs = GetFEModel()->GetDOFS();
	int MAX_CDOFS = fedofs.GetVariableSize("concentration");

	if (strcmp(p.name(), "vR") == 0)
	{
		if (strcmp(szatt, "sbm") == 0)
		{
			int id = atoi(szval) - 1;
			if (id < 0) return false;
			SetStoichiometricCoefficient(m_sbmR, id, m_vRtmp);
			return true;
		}
		if (strcmp(szatt, "sol") == 0)
		{
			int id = atoi(szval) - 1;
			if ((id < 0) || (id >= MAX_CDOFS)) return false;
			SetStoichiometricCoefficient(m_solR, id, m_vRtmp);
			return true;
		}
	}
	else if (strcmp(p.name(), "vP") == 0)
	{
		if (strcmp(szatt, "sbm") == 0)
		{
			int id = atoi(szval) - 1;
			if (id < 0) return false;
			SetStoichiometricCoefficient(m_sbmP, id, m_vPtmp);
			return true;
		}
		if (strcmp(szatt, "sol") == 0)
		{
			int id = atoi(szval) - 1;
			if ((id < 0) || (id >= MAX_CDOFS)) return false;
			SetStoichiometricCoefficient(m_solP, id, m_vPtmp);
			return true;
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Data serialization
void FECellChemicalReaction::Serialize(DumpStream& ar)
{
	FEMaterial::Serialize(ar);

	if (ar.IsShallow() == false)
	{
		if (ar.IsSaving())
		{
			itrmap p;
			ar << m_nsol << m_vR << m_vP << m_v << m_Vovr;
			ar << (int)m_solR.size();
			for (p = m_solR.begin(); p != m_solR.end(); ++p) { ar << p->first; ar << p->second; }
			ar << (int)m_solP.size();
			for (p = m_solP.begin(); p != m_solP.end(); ++p) { ar << p->first; ar << p->second; }
			ar << (int)m_sbmR.size();
			for (p = m_sbmR.begin(); p != m_sbmR.end(); ++p) { ar << p->first; ar << p->second; }
			ar << (int)m_sbmP.size();
			for (p = m_sbmP.begin(); p != m_sbmP.end(); ++p) { ar << p->first; ar << p->second; }
		}
		else
		{
			// restore pointers
			if (m_pFwd) m_pFwd->m_pReact = this;
			if (m_pRev) m_pRev->m_pReact = this;

			ar >> m_nsol >> m_vR >> m_vP >> m_v >> m_Vovr;
			int size, id, vR;
			ar >> size;
			for (int i = 0; i<size; ++i)
			{
				ar >> id; ar >> vR;
				SetStoichiometricCoefficient(m_solR, id, vR);
			}
			ar >> size;
			for (int i = 0; i<size; ++i)
			{
				ar >> id; ar >> vR;
				SetStoichiometricCoefficient(m_solP, id, vR);
			}
			ar >> size;
			for (int i = 0; i<size; ++i)
			{
				ar >> id; ar >> vR;
				SetStoichiometricCoefficient(m_sbmR, id, vR);
			}
			ar >> size;
			for (int i = 0; i<size; ++i)
			{
				ar >> id; ar >> vR;
				SetStoichiometricCoefficient(m_sbmP, id, vR);
			}
		}
	}
}

double FECellMassActionForward::ReactionSupply(FECell* cell)
{
	// get reaction rate
	double kF = m_pFwd->ReactionRate(cell);

	// evaluate the reaction molar supply
	double zhat = kF;

	// start with contribution from solutes

	/*int nsol = (int)spt.m_ca.size();*/
	int nsol = cell->Solutes.size();
	for (int i = 0; i < nsol; ++i) {
		int vR = m_vR[i];
		if (vR > 0) {
			/*double c = spt.m_ca[i];*/
			double c = cell->Solutes[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}
	// add contribution of solid-bound molecules
	/*const int nsbm = (int)spt.m_sbmr.size();*/
	const int nsbm = cell->SBMs.size();
	for (int i = 0; i < nsbm; ++i) {
		int vR = m_vR[nsol + i];
		if (vR > 0) {
			/*double c = m_pMP->SBMConcentration(pt, i);*/
			double c = cell->SBMs[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}

	return zhat;
}

double FECellMassActionForwardEffective::ReactionSupply(FECell* cell)
{
	// get reaction rate
	double kF = m_pFwd->ReactionRate(cell);

	// evaluate the reaction molar supply
	double zhat = kF;

	// start with contribution from solutes
	int nsol = 0;
	/*nsol = (int)spt.m_c.size();*/
	nsol = cell->Solutes.size();
	for (int i = 0; i<nsol; ++i) {
		int vR = m_vR[i];
		if (vR > 0) {
			double c = 0;
				/*c = spt.m_c[i];*/
			c = cell->Solutes[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}

	// add contribution of solid-bound molecules
	/*const int nsbm = (int)spt.m_sbmr.size();*/
	const int nsbm = cell->SBMs.size();
	for (int i = 0; i<nsbm; ++i) {
		int vR = m_vR[nsol + i];
		if (vR > 0) {
			/*double c = m_pMP->SBMConcentration(pt, i);*/
			double c = cell->SBMs[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}

	return zhat;
}

double FECellMassActionReversible::FwdReactionSupply(FECell* cell)
{

	// get forward reaction rate
	double k = m_pFwd->ReactionRate(cell);

	// evaluate the reaction molar supply
	double zhat = k;

	int nsol = 0;

	// start with contribution from solutes
	/*nsol = (int)spt.m_ca.size();*/
	nsol = cell->Solutes.size();
	for (int i = 0; i<nsol; ++i) {
		int vR = m_vR[i];
		if (vR > 0) {
			/*double c = spt.m_ca[i];*/
			double c = cell->Solutes[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}

	// add contribution of solid-bound molecules
	/*const int nsbm = (int)spt.m_sbmr.size();*/
	const int nsbm = cell->SBMs.size();
	for (int i = 0; i<nsbm; ++i) {
		int vR = m_vR[nsol + i];
		if (vR > 0) {
			double c = cell->SBMs[i]->GetInt();
			/*double c = m_pMP->SBMConcentration(pt, i);*/
			zhat *= pow(c, vR);
		}
	}

	return zhat;
}

double FECellMassActionReversible::RevReactionSupply(FECell* cell)
{

	// get forward reaction rate
	double k = m_pRev->ReactionRate(cell);

	// evaluate the reaction molar supply
	double zhat = k;

	// start with contribution from solutes
	int nsol = 0;
	/*nsol = (int)spt.m_ca.size();*/
	nsol = cell->Solutes.size();
	for (int i = 0; i<nsol; ++i) {
		int vP = m_vP[i];
		if (vP > 0) {
			/*double c = spt.m_ca[i];*/
			double c = cell->Solutes[i]->GetInt();
			zhat *= pow(c, vP);
		}
	}
	// add contribution of solid-bound molecules
	/*const int nsbm = (int)spt.m_sbmr.size();*/
	const int nsbm = cell->SBMs.size();
	for (int i = 0; i<nsbm; ++i) {
		int vP = m_vP[nsol + i];
		if (vP > 0) {
			/*double c = m_pMP->SBMConcentration(pt, i);*/
			double c = cell->SBMs[i]->GetInt();
			zhat *= pow(c, vP);
		}
	}

	return zhat;
}

double FECellMassActionReversible::ReactionSupply(FECell* cell)
{
	double zhatF = FwdReactionSupply(cell);
	double zhatR = RevReactionSupply(cell);
	return zhatF - zhatR;
}

double FECellMassActionReversibleEffective::FwdReactionSupply(FECell* cell)
{
	// get forward reaction rate
	double k = m_pFwd->ReactionRate(cell);

	// evaluate the reaction molar supply
	double zhat = k;

	// start with contribution from solutes
	int nsol = 0;
		/*nsol = (int)spt.m_c.size();*/
	nsol = cell->Solutes.size();
	for (int i = 0; i<nsol; ++i) {
		int vR = m_vR[i];
		if (vR > 0) {
			double c = cell->Solutes[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}
	// add contribution of solid-bound molecules
		// add contribution of solid-bound molecules
	/*const int nsbm = (int)spt.m_sbmr.size();*/
	const int nsbm = cell->SBMs.size();
	for (int i = 0; i<nsbm; ++i) {
		int vR = m_vR[nsol + i];
		if (vR > 0) {
			/*double c = m_pMP->SBMConcentration(pt, i);*/
			double c = cell->SBMs[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}

	return zhat;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FECellMassActionReversibleEffective::RevReactionSupply(FECell* cell)
{

	// get forward reaction rate
	double k = m_pRev->ReactionRate(cell);

	// evaluate the reaction molar supply
	double zhat = k;

	// start with contribution from solutes
	int nsol = 0;
	nsol = cell->Solutes.size();
	/*nsol = (int)spt.m_c.size();*/
	for (int i = 0; i<nsol; ++i) {
		int vP = m_vP[i];
		if (vP > 0) {
			double c = cell->Solutes[i]->GetInt();
			/*double c = spt.m_c[i];*/
			zhat *= pow(c, vP);
		}
	}

	// add contribution of solid-bound molecules
		// add contribution of solid-bound molecules

	/*const int nsbm = (int)spt.m_sbmr.size();*/
	const int nsbm = cell->SBMs.size();
	for (int i = 0; i<nsbm; ++i) {
		int vP = m_vP[nsol + i];
		if (vP > 0) {
			double c = cell->SBMs[i]->GetInt();
			/*double c = m_pMP->SBMConcentration(pt, i);*/
			zhat *= pow(c, vP);
		}
	}

	return zhat;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FECellMassActionReversibleEffective::ReactionSupply(FECell* cell)
{
	double zhatF = FwdReactionSupply(cell);
	double zhatR = RevReactionSupply(cell);
	return zhatF - zhatR;
}

// define the material parameters
BEGIN_FECORE_CLASS(FECellMichaelisMenten, FECellChemicalReaction)
ADD_PARAMETER(m_Km, "Km");
ADD_PARAMETER(m_c0, "c0");
END_FECORE_CLASS();

//! data initialization and checking
bool FECellMichaelisMenten::Init()
{
	// Initialize base class
	if (FECellChemicalReaction::Init() == false) return false;

	// there is only one reactant and one product in a Michaelis-Menten reaction
	if (m_solR.size() + m_sbmR.size() > 1) {
		feLogError("Provide only one vR for this reaction");
		return false;
	}

	if (m_solP.size() + m_sbmP.size() > 1) {
		feLogError("Provide only one vP for this reaction");
		return false;
	}

	if (m_c0 < 0) {
		feLogError("c0 must be positive");
		return false;
	}

	const int ntot = (int)m_v.size();
	for (int itot = 0; itot<ntot; itot++) {
		if (m_vR[itot] > 0) m_Rid = itot;
		if (m_vP[itot] > 0) m_Pid = itot;
	}

	if (m_Rid == -1) {
		feLogError("Provide vR for the reactant");
		return false;
	}

	// check if reactant is a solute or a solid-bound molecule
	if (m_Rid >= m_nsol) m_Rtype = true;

	return true;
}

//! molar supply at material point
double FECellMichaelisMenten::ReactionSupply(FECell* cell)
{

	// get reaction rate
	double Vmax = m_pFwd->ReactionRate(cell);
	double c = 0.0;
	if (m_Rtype) {
		c = cell->SBMs[m_Rid]->GetInt();
		/*c = m_pMP->SBMConcentration(cell, m_Rid);*/
	}
	else {
		c = cell->SBMs[m_Rid]->GetInt();
		/*c = spt.m_ca[m_Rid];*/
	}

	double zhat = 0;
	if (c > m_c0) zhat = Vmax*c / (m_Km + c);

	return zhat;
}

