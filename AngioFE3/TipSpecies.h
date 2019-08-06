#pragma once
#include <FECore/BC.h>
#include <FECore/FEMaterial.h>
#include "AngioElement.h"

/*class FEAngio;

//! (Incomplete) Prescribed dof to allow tips to deposit chemicals within a multiphasic material
class TipSpecies: public FEMaterial
{
public:
	//! constructor
	explicit TipSpecies(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~TipSpecies() {}
	int GetSBMID() { return SBM_ID; }
	double GetPR() { return prod_rate; }
	//! idealy changes how the concentrations are updated to be based on the location of tips
	//void Update();
protected:
	DECLARE_PARAMETER_LIST();
private:
	int SBM_ID = 0;
	double prod_rate = 0.0;
};

class TipSpeciesManager : public FEMaterial
{
public:
	explicit TipSpeciesManager(FEModel* pfem) : FEMaterial(pfem) { AddProperty(&tip_species_prop, "tip_species_prop"); tip_species_prop.m_brequired = false; }
	virtual ~TipSpeciesManager() {}
	FEVecPropertyT<TipSpecies> tip_species_prop;
protected:
private:
};
*/