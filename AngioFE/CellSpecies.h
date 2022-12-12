#pragma once
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FEMaterial.h>
#include "AngioElement.h"
#include <FEBioMix/FESolutePointSource.h>
#include <FEBioMix/FESBMPointSource.h>
#include <FEBioMix/FEChemicalReaction.h>
#include <FEBioMix/FEReaction.h>
#include <map>
#include <FEBioMix/FESolute.h>

//! Forward declaration
class FECellChemicalReaction;

//! (Incomplete) Prescribed dof to allow tips to deposit chemicals within a multiphasic material
class CellSBM : public FEMaterialProperty
{
	FECORE_BASE_CLASS(CellSBM)

public:
	//! constructor
	explicit CellSBM(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~CellSBM() {}
	int GetSBMID() { return SBM_ID; }
	void SetSBMID(int r) { SBM_ID = r; }
	double GetInt() { return c_SBM; }
	void SetInt(double r) { c_SBM = r; }
	double GetSBMhat() { return SBMhat; }
	void SetSBMhat(double r) { SBMhat = r; }
	void SetSBMhatp(double r) { SBMhatp = r; }
	double GetSBMhatp() { return SBMhatp; }
	void AddSBMhat(double r) { SBMhat = SBMhat + r; }
	void AddSBMhatp(double r) { SBMhatp = SBMhatp + r; }
	void SetSBMPRhat(double r) { SBMPRhat = r; }
	void AddSBMPRhat(double r) { SBMPRhat = SBMPRhat + r; }
	double GetSBMPRhat() { return SBMPRhat; }
	FESBMPointSource* CellSBMPS;
	double GetPR() { return CellSBMPS->GetRate(); }
	void SetPR(double r) { CellSBMPS->SetRate(r); }
	void AddPR(double r) { CellSBMPS->SetRate(CellSBMPS->GetRate() + r); }
	double c_flux = 0;
	double MolarMass() { return m_M; }
	double Density() { return m_rhoT; }
	int ChargeNumber() { return m_z; }
	void SetMolarMass(double r) { m_M = r; }
	void SetDensity(double r) { m_rhoT = r; }
	void SetCharge(int r) { m_z = r; }
	double GetTolScale() { return tol_scale; }
	void ScaleTolScale(double m) { tol_scale = tol_scale * m; }
	void ResetTolScale() { tol_scale = 1.0; }
	FESBMData* SBMData;
	FESBMData* FindSBMData(int nid);
	//! idealy changes how the concentrations are updated to be based on the location of tips
	//void Update();
protected:
	DECLARE_FECORE_CLASS();
private:
	int SBM_ID = -1;
	// initial number of moles
	double c_SBM = 0;			// SBM number/concentration
	double SBMhatp = 0;			// Previous SBM reaction supply
	double SBMhat = 0;			// New SBM reaction supply
	double SBMPRhat = 0;		// SBM point source reaction supply
	double m_M = 1;
	double m_rhoT = 1;
	int m_z = 0;
	double tol_scale = 1.0;
};

//! (Incomplete) Prescribed dof to allow tips to deposit chemicals within a multiphasic material
class CellSolute : public FEMaterialProperty
{
	FECORE_BASE_CLASS(CellSolute)

public:
	//! constructor
	explicit CellSolute(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~CellSolute() {}
	int GetSoluteID() { return Solute_ID; }
	void SetSoluteID(int i) { Solute_ID = i; }
	FESolutePointSource* CellSolutePS = nullptr;
	double GetPR() { return CellSolutePS->GetRate(); }
	void SetPR(double r) { CellSolutePS->SetRate(r); }
	void AddPR(double r) { CellSolutePS->SetRate(CellSolutePS->GetRate() + r); }
	double c_flux = 0;
	double GetInt() { return c_Sol; }
	double GetSolhat() { return Solhat; }
	double GetSolhatp() { return Solhatp; }
	void SetInt(double r) { c_Sol = r; }
	void SetSolhat(double r) { Solhat = r; }
	void SetSolhatp(double r) { Solhatp = r; }
	void AddSolhat(double r) { Solhat = Solhat + r; }
	void AddSolhatp(double r) { Solhatp = Solhatp + r; }
	void SetSolPRhat(double r) { SolPRhat = r; }
	void AddSolPRhat(double r) { SolPRhat = SolPRhat + r; }
	double GetSolPRhat() { return SolPRhat; }
	void SetMolarMass(double r) { m_M = r; }
	void SetDensity(double r) { m_rhoT = r; }
	void SetCharge(int r) { m_z = r; }
	double MolarMass() { return m_M; }
	double Density() { return m_rhoT; }
	int ChargeNumber() { return m_z; }
	FESoluteData* SolData;
	FESoluteData* FindSoluteData(int nid);
	double GetTolScale() { return tol_scale; }
	void ScaleTolScale(double m) { tol_scale = tol_scale * m; }
	void ResetTolScale() { tol_scale = 1.0; }
	//! idealy changes how the concentrations are updated to be based on the location of tips
	//void Update();
protected:
	DECLARE_FECORE_CLASS();
private:
	int Solute_ID = -1;
	double c_Sol = 0;				// solute concentration/number
	double Solhat = 0;				// Solute reaction supply
	double Solhatp = 0;				// Previous solute reaction supply
	double SolPRhat = 0;			// Solute point source reaction supply
	double m_M = 0;
	double m_rhoT = 1;
	int m_z = 1;
	double tol_scale = 1.0;
};

class CellSpeciesManager : public FEMaterialProperty
{
	FECORE_BASE_CLASS(CellSpeciesManager)

public:
	explicit CellSpeciesManager(FEModel* pfem) : FEMaterialProperty(pfem)
	{
		AddClassProperty(this, &cell_solute_prop, "cell_solute_prop", FEProperty::Optional);
		AddClassProperty(this, &cell_SBM_prop, "cell_SBM_prop", FEProperty::Optional);
	}
	virtual ~CellSpeciesManager() {}
	std::vector<CellSolute*> cell_solute_prop;
	std::vector<CellSBM*> cell_SBM_prop;
protected:
private:
};

class CellReactionManager : public FEMaterialProperty
{
	FECORE_BASE_CLASS(CellReactionManager)

public:
	explicit CellReactionManager(FEModel* pfem) : FEMaterialProperty(pfem)
	{
		AddClassProperty(this, &cell_reaction, "cell_reaction", FEProperty::Optional);
	}
	virtual ~CellReactionManager() {}
	std::vector<FECellChemicalReaction*>	cell_reaction;
public:
	bool Init();
protected:
private:
};

class FECellReaction : public FEMaterialProperty
{
	FECORE_BASE_CLASS(FECellReaction)

public:
	//! constructor
	FECellReaction(FEModel* pfem);

	//! initialization
	bool Init() override;

public:
	//! set stoichiometric coefficients
	void SetStoichiometricCoefficient(intmap& RP, int id, int v) { RP.insert(std::pair<int, int>(id, v)); }
	//! set the cell
	void SetCell(FECell* cell) { m_cell = cell; }
public:
	FECell* m_cell;	//!< pointer to the cell where reaction occurs
	FEModel* m_pfem;
};

class FECellReactionRate : public FEMaterialProperty
{
	FECORE_BASE_CLASS(FECellReactionRate)

public:
	//! constructor
	FECellReactionRate(FEModel* pfem) : FEMaterialProperty(pfem) {}

	//! data initializatino and checking
	bool Init() override;

	//! reaction rate 
	virtual double ReactionRate() = 0;

	void SetCell(FECell* cell) { m_cell = cell; }

public:
	FECellReaction* m_pReact;	//!< pointer to parent reaction
	FECell* m_cell;				//!< cell containing this

};

class FECellReactionRateConst : public FECellReactionRate
{
public:
	//! constructor
	FECellReactionRateConst(FEModel* pfem) : FECellReactionRate(pfem) { m_k = 0; }

	//! data initializatino and checking
	bool Init() override;

	//! reaction rate at material point
	double ReactionRate() override { return m_k; }

	//! set the cell
	void SetCell(FECell* cell) { m_cell = cell; }

public:
	double m_k;		//!< reaction rate
	FECell* m_cell;				//!< cell containing this

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
//	void SetParameter(FEParam& p) override;

//	bool SetParameterAttribute(FEParam& p, const char* szatt, const char* szval) override;

	//! set the forward reaction rate
	void SetForwardReactionRate(FECellReactionRate* pfwd) { m_pFwd = pfwd; }

	//! set the reverse reaction rate
	void SetReverseReactionRate(FECellReactionRate* prev) { m_pRev = prev; }

	//! set the cell
	void SetCell(FECell* cell) { m_cell = cell; }

public:
	//! molar supply at material point
	virtual double ReactionSupply(FECell* cell) = 0;

public:
	//! Serialization
	void Serialize(DumpStream& ar) override;

public:
	FECellReactionRate* m_pFwd;        //!< pointer to forward reaction rate
	FECellReactionRate* m_pRev;        //!< pointer to reverse reaction rate
	FECell* m_cell;							//!< cell containing this

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
//	bool            m_Vovr;     //!< override flag for m_Vbar
	int				m_vRtmp;	//!< helper variable for reading in stoichiometric coefficients for reactants
	int				m_vPtmp;	//!< helper variable for reading in stoichiometric coefficients for products

	DECLARE_FECORE_CLASS();
};

class FECellMassActionForward : public FECellChemicalReaction
{
public:
	//! constructor
	FECellMassActionForward(FEModel* pfem) : FECellChemicalReaction(pfem) {}

	//! data initializatino and checking
	bool Init() override;

	//! molar supply
	double ReactionSupply(FECell* cell);
	//! set the cell
	void SetCell(FECell* cell) { m_cell = cell; }

public:
	FECell* m_cell;				//!< cell containing this
};

class FECellMassActionForwardConstant : public FECellChemicalReaction
{
public:
	//! constructor
	FECellMassActionForwardConstant(FEModel* pfem) : FECellChemicalReaction(pfem) {}

	//! data initializatino and checking
	bool Init() override;

	//! molar supply
	double ReactionSupply(FECell* cell);
	//! set the cell
	void SetCell(FECell* cell) { m_cell = cell; }

public:
	FECell* m_cell;				//!< cell containing this
};

class FECellMassActionForwardEffective : public FECellChemicalReaction
{
public:
	//! constructor
	FECellMassActionForwardEffective(FEModel* pfem) : FECellChemicalReaction(pfem) {}

	//! data initializatino and checking
	bool Init() override;

	//! molar supply
	double ReactionSupply(FECell* cell);

	//! set the cell
	void SetCell(FECell* cell) { m_cell = cell; }

public:
	FECell* m_cell;				//!< cell containing this
};

class FECellMassActionReversible : public FECellChemicalReaction
{
public:
	//! constructor
	FECellMassActionReversible(FEModel* pfem) : FECellChemicalReaction(pfem) {}

	//! data initializatino and checking
	bool Init() override;

	//! molar supply
	double ReactionSupply(FECell* cell);

	//! molar supply of Fwd
	double FwdReactionSupply(FECell* cell);

	//! molar supply of Rev
	double RevReactionSupply(FECell* cell);

	//! set the cell
	void SetCell(FECell* cell) { m_cell = cell; }

public:
	FECell* m_cell;				//!< cell containing this
};

class FECellMassActionReversibleEffective : public FECellChemicalReaction
{
public:
	//! constructor
	FECellMassActionReversibleEffective(FEModel* pfem) : FECellChemicalReaction(pfem) {}

	//! data initializatino and checking
	bool Init() override;

	//! molar supply
	double ReactionSupply(FECell* cell);

	//! molar supply of Fwd
	double FwdReactionSupply(FECell* cell);

	//! molar supply of Rev
	double RevReactionSupply(FECell* cell);

	//! set the cell
	void SetCell(FECell* cell) { m_cell = cell; }

public:
	FECell* m_cell;				//!< cell containing this
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

	//! set the cell
	void SetCell(FECell* cell) { m_cell = cell; }

public:
	double	m_Km;			//!< concentration at which half-maximum rate occurs
	int		m_Rid;			//!< local id of reactant
	int		m_Pid;			//!< local id of product
	bool	m_Rtype;		//!< flag for reactant type (solute = false, sbm = true)
	double	m_c0;			//!< minimum reactant concentration to trigger reaction
	FECell* m_cell;				//!< cell containing this

							// declare parameter list
	DECLARE_FECORE_CLASS();

};

//class CellSBMManager : public FEMaterial

typedef std::map<int, int> intmap;
typedef std::map<int, int>::iterator itrmap;

class FECellInternalization : public FECellChemicalReaction
{
public:
	//! constructor
	FECellInternalization(FEModel* pfem) : FECellChemicalReaction(pfem) { m_Ki = 0; }

	bool Init() override;

	//! reaction rate at material point
	double ReactionSupply(FECell* cell);

	//! set the cell
	void SetCell(FECell* cell) { m_cell = cell; }

public:
	double m_Ki;					//!< reaction rate
	FECell* m_cell;					//!< cell containing this
	//FESolutesMaterialPoint* m_pMP;	//!< pointer to solutes material point for internalizing species from

	DECLARE_FECORE_CLASS();

};

class FECellInternalizationConstant : public FECellChemicalReaction
{
public:
	//! constructor
	FECellInternalizationConstant(FEModel* pfem) : FECellChemicalReaction(pfem) { m_Ki = 0; }

	bool Init() override;

	//! reaction rate at material point
	double ReactionSupply(FECell* cell);

	//! set the cell
	void SetCell(FECell* cell) { m_cell = cell; }

public:
	double m_Ki;					//!< reaction rate
	FECell* m_cell;					//!< cell containing this
	//FESolutesMaterialPoint* m_pMP;	//!< pointer to solutes material point for internalizing species from

	DECLARE_FECORE_CLASS();

};

class FECellSecretion : public FECellChemicalReaction
{
public:
	bool Init() override;

	//! constructor
	FECellSecretion(FEModel* pfem) : FECellChemicalReaction(pfem) { m_s = 0; }

	//! reaction rate at material point
	double ReactionSupply(FECell* cell);

	//! set the cell
	void SetCell(FECell* cell) { m_cell = cell; }

public:
	double m_s;					//!< reaction rate
	FECell* m_cell;					//!< cell containing this
	//FESolutesMaterialPoint* m_pMP;	//!< pointer to solutes material point for internalizing species from

	DECLARE_FECORE_CLASS();

};

class FECellSecretionConstant : public FECellChemicalReaction
{
public:
	bool Init() override;

	//! constructor
	FECellSecretionConstant(FEModel* pfem) : FECellChemicalReaction(pfem) { m_s = 0; }

	//! reaction rate at material point
	double ReactionSupply(FECell* cell);

	//! set the cell
	void SetCell(FECell* cell) { m_cell = cell; }

public:
	double m_s;					//!< reaction rate
	FECell* m_cell;					//!< cell containing this
	//FESolutesMaterialPoint* m_pMP;	//!< pointer to solutes material point for internalizing species from

	DECLARE_FECORE_CLASS();

};