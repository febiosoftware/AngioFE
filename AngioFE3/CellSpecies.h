#pragma once
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FEMaterial.h>
#include "AngioElement.h"
#include <FEBioMix\FESolutePointSource.h>
#include <FEBioMix\FESBMPointSource.h>
#include <FEBioMix\FEChemicalReaction.h>
#include <FEBioMix\FEReaction.h>
#include <map>

//! (Incomplete) Prescribed dof to allow tips to deposit chemicals within a multiphasic material
class CellSBM: public FEMaterial
{
public:
	//! constructor
	explicit CellSBM(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~CellSBM() {}
	int GetSBMID() { return SBM_ID; }
	void SetSBMID(int i) { SBM_ID = i; }
	double GetPR() { return SBM_prod_rate; }
	void SetPR(double r) { SBM_prod_rate = r; }
	double GetInt() { return n_SBM; }
	void SetInt(double r) { n_SBM = r; }
	double GetSBMhat() { return SBMhat; }
	void SetSBMhat(double r) { SBMhat = r; }
	void SetSBMhatp(double r) { SBMhatp = r; }
	double GetSBMhatp() { return SBMhatp; }
	void AddSBMhat(double r) { SBMhat = SBMhat + r; }
	void AddSBMhatp(double r) { SBMhatp = SBMhatp + r; }
	FESBMPointSource* CellSBMPS;
	//! idealy changes how the concentrations are updated to be based on the location of tips
	//void Update();
protected:
	DECLARE_FECORE_CLASS();
private:
	int SBM_ID = -1;
	double SBM_prod_rate = -1.0;
	// initial number of moles
	double n_SBM = 0;
	double SBMhatp = 0;
	double SBMhat = 0;
};

//! (Incomplete) Prescribed dof to allow tips to deposit chemicals within a multiphasic material
class CellSolute : public FEMaterial
{
public:
	//! constructor
	explicit CellSolute(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~CellSolute() {}
	int GetSoluteID() { return Solute_ID; }
	void SetSoluteID(int i) { Solute_ID = i; }
	double GetPR() { return Solute_prod_rate; }
	void SetPR(double r) { Solute_prod_rate = r; }
	FESolutePointSource* CellSolutePS = nullptr;
	double GetInt() { return n_Solute; }
	double GetSolhat() { return Solhat; }
	double GetSolhatp() { return Solhatp; }
	void SetInt(double r) { n_Solute = r; }
	void SetSolhat(double r) { Solhat = r; }
	void SetSolhatp(double r) { Solhatp = r; }
	void AddSolhat(double r) { Solhat = Solhat + r; }
	void AddSolhatp(double r) { Solhatp = Solhatp + r; }

	//! idealy changes how the concentrations are updated to be based on the location of tips
	//void Update();
protected:
	DECLARE_FECORE_CLASS();
private:
	int Solute_ID = -1;
	double Solute_prod_rate = -1.0;
	double n_Solute = 0;
	double Solhat = 0;
	double Solhatp = 0;
};

class CellSpeciesManager : public FEMaterial
{
public:
	explicit CellSpeciesManager(FEModel* pfem) : FEMaterial(pfem) 
	{/*	AddClassProperty(this, &cell_species_prop, "cell_species_prop", FEProperty::Optional); */
		AddClassProperty(this, &cell_solute_prop, "cell_solute_prop", FEProperty::Optional);
		AddClassProperty(this, &cell_SBM_prop, "cell_SBM_prop", FEProperty::Optional);
	}
	virtual ~CellSpeciesManager() {}
	/*std::vector<CellSpecies*>	cell_species_prop;*/
	std::vector<CellSolute*> cell_solute_prop;
	std::vector<CellSBM*> cell_SBM_prop;
protected:
private:
};

class FECellReaction : public FEMaterial
{
public:
	//! constructor
	FECellReaction(FEModel* pfem);

	//! initialization
	bool Init() override;

public:
	//! set stoichiometric coefficients
	void SetStoichiometricCoefficient(intmap& RP, int id, int v) { RP.insert(std::pair<int, int>(id, v)); }
public:
	FECell* m_cell;	//!< pointer to the cell where reaction occurs
};

class FECellReactionRate : public FEMaterial
{
public:
	//! constructor
	FECellReactionRate(FEModel* pfem) : FEMaterial(pfem) {}

	//! reaction rate 
	virtual double ReactionRate(FECell* cell) = 0;

public:
	FECellReaction* m_pReact;	//!< pointer to parent reaction

};

class FECellReactionRateConst : public FECellReactionRate
{
public:
	//! constructor
	FECellReactionRateConst(FEModel* pfem) : FECellReactionRate(pfem) { m_k = 0; }

	//! reaction rate at material point
	double ReactionRate(FECell* cell) override { return m_k; }

public:
	double m_k;		//!< reaction rate

	DECLARE_FECORE_CLASS();

};

class FECellChemicalReaction : public FECellReaction
{
public:
	//! constructor
	FECellChemicalReaction(FEModel* pfem);

	//! initialization
	bool Init() override;

public:
	void SetParameter(FEParam& p) override;

	bool SetParameterAttribute(FEParam& p, const char* szatt, const char* szval) override;

	//! set the forward reaction rate
	void SetForwardReactionRate(FECellReactionRate* pfwd) { m_pFwd = pfwd; }

	//! set the reverse reaction rate
	void SetReverseReactionRate(FECellReactionRate* prev) { m_pRev = prev; }

public:
	//! molar supply at material point
	virtual double ReactionSupply(FECell* cell) = 0;

public:
	//! Serialization
	void Serialize(DumpStream& ar) override;

public:
	FECellReactionRate*    m_pFwd;        //!< pointer to forward reaction rate
	FECellReactionRate*    m_pRev;        //!< pointer to reverse reaction rate

public:
	intmap			m_solR;		//!< stoichiometric coefficients of solute reactants (input)
	intmap			m_solP;		//!< stoichiometric coefficients of solute products (input)
	intmap			m_sbmR;		//!< stoichiometric coefficients of solid-bound reactants (input)
	intmap			m_sbmP;		//!< stoichiometric coefficients of solid-bound products (input)

public:
	int				m_nsol;		//!< number of solutes in the mixture
	vector<int>		m_vR;		//!< stoichiometric coefficients of reactants
	vector<int>		m_vP;		//!< stoichiometric coefficients of products
	vector<int>		m_v;		//!< net stoichiometric coefficients of reactants and products
	double          m_Vbar;     //!< weighted molar volume of reactants and products
	bool            m_Vovr;     //!< override flag for m_Vbar
	int				m_vRtmp;	//!< helper variable for reading in stoichiometric coefficients for reactants
	int				m_vPtmp;	//!< helper variable for reading in stoichiometric coefficients for products

	DECLARE_FECORE_CLASS();
};

class FECellMassActionForward : public FECellChemicalReaction
{
public:
	//! constructor
	FECellMassActionForward(FEModel* pfem) : FECellChemicalReaction(pfem) {}

	//! molar supply
	double ReactionSupply(FECell* cell);
};

class FECellMassActionForwardEffective : public FECellChemicalReaction
{
public:
	//! constructor
	FECellMassActionForwardEffective(FEModel* pfem) : FECellChemicalReaction(pfem) {}

	//! molar supply
	double ReactionSupply(FECell* cell);
};

class FECellMassActionReversible : public FECellChemicalReaction
{
public:
	//! constructor
	FECellMassActionReversible(FEModel* pfem) : FECellChemicalReaction(pfem) {}

	//! molar supply
	double ReactionSupply(FECell* cell);

	//! molar supply of Fwd
	double FwdReactionSupply(FECell* cell);

	//! molar supply of Rev
	double RevReactionSupply(FECell* cell);
};

class FECellMassActionReversibleEffective : public FECellChemicalReaction
{
public:
	//! constructor
	FECellMassActionReversibleEffective(FEModel* pfem) : FECellChemicalReaction(pfem) {}

	//! molar supply
	double ReactionSupply(FECell* cell);

	//! molar supply of Fwd
	double FwdReactionSupply(FECell* cell);

	//! molar supply of Rev
	double RevReactionSupply(FECell* cell);
};

class FECellMichaelisMenten : public FECellChemicalReaction
{
public:
	//! constructor
	FECellMichaelisMenten(FEModel* pfem) : FECellChemicalReaction(pfem) { m_Rid = m_Pid = -1; m_Km = m_c0 = 0; m_Rtype = false; }

	//! data initializatino and checking
	bool Init() override;

	//! molar supply
	double ReactionSupply(FECell* cell) override;

public:
	double	m_Km;			//!< concentration at which half-maximum rate occurs
	int		m_Rid;			//!< local id of reactant
	int		m_Pid;			//!< local id of product
	bool	m_Rtype;		//!< flag for reactant type (solute = false, sbm = true)
	double	m_c0;			//!< minimum reactant concentration to trigger reaction

							// declare parameter list
	DECLARE_FECORE_CLASS();

};

class CellReactionManager : public FEMaterial
{
public:
	explicit CellReactionManager(FEModel* pfem) : FEMaterial(pfem)
	{
		AddClassProperty(this, &cell_reaction, "cell_reaction", FEProperty::Optional);
	}
	virtual ~CellReactionManager() {}
	std::vector<FECellChemicalReaction*>	cell_reaction;
protected:
private:
};

//class CellSBMManager : public FEMaterial
//{
//public:
//	explicit CellSBMManager(FEModel* pfem) : FEMaterial(pfem) { AddClassProperty(this, &cell_SBM_prop, "cell_SBM_prop", FEProperty::Optional); }
//	virtual ~CellSBMManager() {}
//	std::vector<CellSBM*>	cell_SBM_prop;	
//protected:
//private:
//};

//class CellSoluteManager : public FEMaterial
//{
//public:
//	explicit CellSoluteManager(FEModel* pfem) : FEMaterial(pfem) { AddClassProperty(this, &cell_solute_prop, "cell_solute_prop", FEProperty::Optional); }
//	virtual ~CellSoluteManager() {}
//	std::vector<CellSolute*>	cell_solute_prop;	//!< pointers to elastic materials
//													//TipSpecies* tip_species_prop;
//protected:
//private:
//};

typedef std::map<int, int> intmap;
typedef std::map<int, int>::iterator itrmap;

