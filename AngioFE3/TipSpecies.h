#pragma once
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FEMaterial.h>
#include "AngioElement.h"

class FEAngio;

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
	DECLARE_FECORE_CLASS();
private:
	int SBM_ID = 0;
	double prod_rate = 0.0;
};

class TipSpeciesManager : public FEMaterial
{
public:
	explicit TipSpeciesManager(FEModel* pfem) : FEMaterial(pfem) { AddClassProperty(this, &tip_species_prop, "tip_species_prop"), FEProperty::Optional; }
	virtual ~TipSpeciesManager() {}
	std::vector<TipSpecies*>	tip_species_prop;	//!< pointers to elastic materials
	//TipSpecies* tip_species_prop;
protected:
private:
};