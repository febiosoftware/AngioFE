#pragma once
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FEMaterial.h>
#include "AngioElement.h"

class FEAngio;

//! (Incomplete) Prescribed dof to allow tips to deposit chemicals within a multiphasic material
class CellSpecies: public FEMaterial
{
public:
	//! constructor
	explicit CellSpecies(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~CellSpecies() {}
	int GetSpeciesID() { return species_ID; }
	double GetPR() { return prod_rate; }
	//! idealy changes how the concentrations are updated to be based on the location of tips
	//void Update();
protected:
	DECLARE_FECORE_CLASS();
private:
	int species_ID = 0;
	double prod_rate = 0.0;
};

class CellSpeciesManager : public FEMaterial
{
public:
	explicit CellSpeciesManager(FEModel* pfem) : FEMaterial(pfem) { AddClassProperty(this, &cell_species_prop, "cell_species_prop"), FEProperty::Optional; }
	virtual ~CellSpeciesManager() {}
	std::vector<CellSpecies*>	cell_species_prop;
protected:
private:
};

class CellSBMManager : public FEMaterial
{
public:
	explicit CellSBMManager(FEModel* pfem) : FEMaterial(pfem) { AddClassProperty(this, &cell_SBM_prop, "cell_SBM_prop"), FEProperty::Optional; }
	virtual ~CellSBMManager() {}
	std::vector<CellSpecies*>	cell_SBM_prop;	
protected:
private:
};

class CellSoluteManager : public FEMaterial
{
public:
	explicit CellSoluteManager(FEModel* pfem) : FEMaterial(pfem) { AddClassProperty(this, &cell_sol_prop, "cell_sol_prop"), FEProperty::Optional; }
	virtual ~CellSoluteManager() {}
	std::vector<CellSpecies*>	cell_sol_prop;	//!< pointers to elastic materials
													//TipSpecies* tip_species_prop;
protected:
private:
};