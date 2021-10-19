#include "CellSpecies.h"
#include "FEAngio.h"
#include "FEAngio.h"
#include "FEAngioMaterial.h"
#include <iostream>
#include <FECore/log.h>

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

}

//! SL: What should this do?
bool FECellReaction::Init() {
	return true;
}

bool CellReactionManager::Init() {

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

FESBMData* CellSBM::FindSBMData(int nid) {
	FEModel& fem = *GetFEModel();
	int N = GetFEModel()->GlobalDataItems();
	for (int i = 0; i < N; ++i)
	{
		FESBMData* psd = dynamic_cast<FESBMData*>(fem.GetGlobalData(i));
		if (psd && (psd->GetID() == nid)) return psd;
	}
	return 0;
}

FESoluteData* CellSolute::FindSoluteData(int nid) {
	FEModel& fem = *GetFEModel();
	int N = GetFEModel()->GlobalDataItems();
	for (int i = 0; i < N; ++i)
	{
		FESoluteData* psd = dynamic_cast<FESoluteData*>(fem.GetGlobalData(i));
		if (psd && (psd->GetID() == nid)) return psd;
	}
	return 0;
}

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
	FEMesh* mesh = &this->GetFEModel()->GetMesh();
	int nsol, nsbm, ntot;
	if (m_cell) {
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
		for (isol = 0; isol < nsol; ++isol) {
			int sid = isol;
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
		for (isbm = 0; isbm < nsbm; ++isbm) {
			int sid = isbm;
			sid = m_cell->SBMs[isbm]->GetSBMID() - 1;
			it = sbmR.find(sid);
			if (it != sbmR.end()) m_vR[nsol + isbm] = it->second;
			it = sbmP.find(sid);
			if (it != sbmP.end()) m_vP[nsol + isbm] = it->second;
		}

		// evaluate the net stoichiometric coefficient
		for (itot = 0; itot < ntot; ++itot) {
			m_v[itot] = m_vP[itot] - m_vR[itot];
		}

		 //evaluate the weighted molar volume of reactants and products
		if (!m_Vovr) {
			m_Vbar = 0;
			for (isol = 0; isol<nsol; ++isol){
				m_Vbar += m_v[isol] * m_cell->Solutes[isol]->MolarMass() / m_cell->Solutes[isol]->Density();
			}
			for (isbm = 0; isbm<nsbm; ++isbm)
				m_Vbar += m_v[nsol + isbm] * m_cell->SBMs[isbm]->MolarMass() / m_cell->SBMs[isbm]->Density();
		}

		 //check that the chemical reaction satisfies electroneutrality
		int znet = 0;
		for (isol = 0; isol<nsol; ++isol){
				znet += m_v[isol] * m_cell->Solutes[isol]->ChargeNumber();
		}
		for (isbm = 0; isbm<nsbm; ++isbm)
			znet += m_v[nsol + isbm] * m_cell->SBMs[isbm]->ChargeNumber();
		if (znet != 0) {
			feLogError("chemical reaction must satisfy electroneutrality");
			return false;
		}
	}

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

double FECellMassActionForward::ReactionSupply()
{
	// get reaction rate
	double kF = m_pFwd->ReactionRate();
	m_cell = m_pFwd->m_cell;
	// evaluate the reaction molar supply
	double zhat = kF;

	// start with contribution from solutes

	int nsol = m_cell->Solutes.size();
	for (int i = 0; i < nsol; ++i) {
		int vR = m_vR[i];
		if (vR > 0) {
			double c = m_cell->Solutes[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}
	// add contribution of solid-bound molecules
	const int nsbm = m_cell->SBMs.size();
	for (int i = 0; i < nsbm; ++i) {
		int vR = m_vR[nsol + i];
		if (vR > 0) {
			double c = m_cell->SBMs[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}

	return zhat;
}

double FECellMassActionForwardEffective::ReactionSupply()
{
	// get reaction rate
	double kF = m_pFwd->ReactionRate();
	m_cell = m_pFwd->m_cell;

	// evaluate the reaction molar supply
	double zhat = kF;

	// start with contribution from solutes
	int nsol = 0;
	nsol = m_cell->Solutes.size();
	for (int i = 0; i<nsol; ++i) {
		int vR = m_vR[i];
		if (vR > 0) {
			double c = 0;
			c = m_cell->Solutes[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}

	// add contribution of solid-bound molecules
	const int nsbm = m_cell->SBMs.size();
	for (int i = 0; i<nsbm; ++i) {
		int vR = m_vR[nsol + i];
		if (vR > 0) {
			double c = m_cell->SBMs[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}

	return zhat;
}

double FECellMassActionReversible::FwdReactionSupply()
{

	// get forward reaction rate
	double k = m_pFwd->ReactionRate();
	m_cell = m_pRev->m_cell;

	// evaluate the reaction molar supply
	double zhat = k;

	int nsol = 0;

	m_cell = m_pFwd->m_cell;

	// start with contribution from solutes
	nsol = this->m_cell->Solutes.size();
	for (int i = 0; i<nsol; ++i) {
		int vR = m_vR[i];
		if (vR > 0) {
			double c = this->m_cell->Solutes[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}

	// add contribution of solid-bound molecules
	const int nsbm = this->m_cell->SBMs.size();
	for (int i = 0; i<nsbm; ++i) {
		int vR = m_vR[nsol + i];
		if (vR > 0) {
			double c = this->m_cell->SBMs[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}

	return zhat;
}

double FECellMassActionReversible::RevReactionSupply()
{

	// get forward reaction rate
	double k = m_pRev->ReactionRate();
	m_cell = m_pRev->m_cell;

	// evaluate the reaction molar supply
	double zhat = k;

	// start with contribution from solutes
	int nsol = 0;
	nsol = m_cell->Solutes.size();
	for (int i = 0; i<nsol; ++i) {
		int vP = m_vP[i];
		if (vP > 0) {
			double c = m_cell->Solutes[i]->GetInt();
			zhat *= pow(c, vP);
		}
	}
	// add contribution of solid-bound molecules
	const int nsbm = m_cell->SBMs.size();
	for (int i = 0; i<nsbm; ++i) {
		int vP = m_vP[nsol + i];
		if (vP > 0) {
			double c = m_cell->SBMs[i]->GetInt();
			zhat *= pow(c, vP);
		}
	}

	return zhat;
}

double FECellMassActionReversible::ReactionSupply()
{
	double zhatF = FwdReactionSupply();
	double zhatR = RevReactionSupply();
	return zhatF - zhatR;
}

double FECellMassActionReversibleEffective::FwdReactionSupply()
{
	// get forward reaction rate
	double k = m_pFwd->ReactionRate();
	m_cell = m_pFwd->m_cell;

	// evaluate the reaction molar supply
	double zhat = k;

	// start with contribution from solutes
	int nsol = 0;
		/*nsol = (int)spt.m_c.size();*/
	nsol = m_cell->Solutes.size();
	for (int i = 0; i<nsol; ++i) {
		int vR = m_vR[i];
		if (vR > 0) {
			double c = m_cell->Solutes[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}
	// add contribution of solid-bound molecules
	// add contribution of solid-bound molecules
	const int nsbm = m_cell->SBMs.size();
	for (int i = 0; i<nsbm; ++i) {
		int vR = m_vR[nsol + i];
		if (vR > 0) {
			double c = m_cell->SBMs[i]->GetInt();
			zhat *= pow(c, vR);
		}
	}

	return zhat;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FECellMassActionReversibleEffective::RevReactionSupply()
{

	// get forward reaction rate
	double k = m_pRev->ReactionRate();
	m_cell = m_pRev->m_cell;

	// evaluate the reaction molar supply
	double zhat = k;

	// start with contribution from solutes
	int nsol = 0;
	nsol = m_cell->Solutes.size();
	/*nsol = (int)spt.m_c.size();*/
	for (int i = 0; i<nsol; ++i) {
		int vP = m_vP[i];
		if (vP > 0) {
			double c = m_cell->Solutes[i]->GetInt();
			zhat *= pow(c, vP);
		}
	}

	// add contribution of solid-bound molecules
		// add contribution of solid-bound molecules

	/*const int nsbm = (int)spt.m_sbmr.size();*/
	const int nsbm = m_cell->SBMs.size();
	for (int i = 0; i<nsbm; ++i) {
		int vP = m_vP[nsol + i];
		if (vP > 0) {
			double c = m_cell->SBMs[i]->GetInt();
			zhat *= pow(c, vP);
		}
	}

	return zhat;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FECellMassActionReversibleEffective::ReactionSupply()
{
	double zhatF = FwdReactionSupply();
	double zhatR = RevReactionSupply();
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
double FECellMichaelisMenten::ReactionSupply()
{

	// get reaction rate
	double Vmax = m_pFwd->ReactionRate();
	m_cell = m_pFwd->m_cell;
	double c = 0.0;
	if (m_Rtype) {
		c = m_cell->SBMs[m_Rid]->GetInt();
	}
	else {
		c = m_cell->Solutes[m_Rid]->GetInt();
	}

	double zhat = 0;
	if (c > m_c0) zhat = Vmax*c / (m_Km + c);

	return zhat;
}

double FECellInternalization::ReactionSupply()
{
	// get reaction rate
	double m_Ki = m_pFwd->ReactionRate();
	m_cell = m_pFwd->m_cell;
	
	//double zhat = m_Ki;
	// start with contribution from solutes

	int nsol = m_cell->Solutes.size();
	for (int isol = 0; isol < nsol; ++isol) {
		int vP = m_vP[isol];
		if (vP > 0) {
			// Cells take away from the matrix
			m_cell->Solutes[isol]->AddPR(-1.0*m_Ki);
			return m_Ki;
		}
	}
	// add contribution of solid-bound molecules
	const int nsbm = m_cell->SBMs.size();
	for (int isbm = 0; isbm < nsbm; ++isbm) {
		int vP = m_vP[nsol + isbm];
		if (vP > 0) {
			m_cell->SBMs[isbm]->AddPR(-1.0*m_Ki);

			// Cells take away from the matrix
			return m_Ki;
		}
	}
	//! Add the species to the cell
	//return zhat;
}

BEGIN_FECORE_CLASS(FECellInternalization, FECellChemicalReaction)
ADD_PARAMETER(m_Ki, "Ki");
END_FECORE_CLASS();

double FECellSecretion::ReactionSupply()
{
	// get reaction rate
	double m_s = m_pFwd->ReactionRate();
	m_cell = m_pFwd->m_cell;

	double zhat = m_s;
	// start with contribution from solutes

	int nsol = m_cell->Solutes.size();
	for (int isol = 0; isol < nsol; ++isol) {
		int vR = m_vR[isol];
		if (vR > 0) {

			double m_c = m_cell->Solutes[isol]->GetInt();
			// Cells take away from the matrix
			zhat *= pow(m_c, vR);
			//double pr = m_cell->Solutes[isol]->GetPR();
			//m_cell->Solutes[isol]->SetPR(pr+zhat);
		}
	}
	// add contribution of solid-bound molecules
	const int nsbm = m_cell->SBMs.size();
	for (int isbm = 0; isbm < nsbm; ++isbm) {
		int vR = m_vR[nsol + isbm];
		if (vR > 0) {
			double m_c = m_cell->SBMs[isbm]->GetInt();
			zhat *= pow(m_c, vR);
		}
	}
	//! Remove the species from the cell
	return zhat;
}

BEGIN_FECORE_CLASS(FECellSecretion, FECellChemicalReaction)
ADD_PARAMETER(m_s, "s");
END_FECORE_CLASS();