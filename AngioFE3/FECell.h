#pragma once
#include <FECore/vec3d.h>
#include "AngioElement.h"
#include <FEBioMix/FESBMPointSource.h>
#include <FEBioMix/FESolutePointSource.h>
#include "FEProbabilityDistribution.h"
#include "Tip.h"
#include "CellSpecies.h"

class Segment;
class FEAngio;

//! a location where growth has occured, or can occur on the next growth step/substep.
class FECell
{
public:
	//! constructor
	FECell() {}
	//! used to copy another tip
	FECell(FECell * other, FEMesh * mesh);
	//! sets the local positions clamps the values to -1 to 1
	void SetLocalPosition(vec3d pos, FEMesh * mesh);
	//! returns the local position
	vec3d GetLocalPosition() const;
	//! The element that contains the local_pos coordinates
	AngioElement * angio_element = nullptr;
	//! Time at which the cell occurs
	double time = 0.0;
	//! Will be initialized to values greater than or equal to zero, unique values per vessel
	int initial_cell_id = -1;
	//! Return the global position of the tip
	vec3d GetPosition(FEMesh * mesh) const;
	//! returns the global position of the tip in the refernce frame
	vec3d GetRefPosition(FEMesh * mesh) const;
	//! Prints information about the tip to the console
	void PrintCellInfo(FEMesh *mesh) const;
	//! Prints information about the tip to the console
	void PrintCellInfo(FEMesh *mesh, std::string title) const;
	//! Cell Species
	/*std::unordered_map<int, CellSBM*> SBMs;
	std::unordered_map<int, CellSolute*> Solutes;*/
	std::vector<CellSBM*> SBMs;
	std::vector<CellSolute*> Solutes;
	std::vector<FEChemicalReaction*> Reactions;
	//std::unordered_map<int, FESolutePointSource*> Species;
	void InitSpecies(FEMesh* mesh);
	void UpdateSpecies(FEMesh * mesh);
	//void InitReactions(FEMesh* mesh);
	double GetReactionSupply(FEChemicalReaction* m_R);
	Tip* ParentTip;
	double eval_time = 0;
	double cell_vol = 1000; // cell volume is 1000 um3 = 1 uL
private:
	vec3d local_pos;
	friend class Tip;
};