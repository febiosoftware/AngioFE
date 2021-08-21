#pragma once
#include <FECore/vec3d.h>
#include "AngioElement.h"
#include <FEBioMix/FESBMPointSource.h>
#include "FEProbabilityDistribution.h"

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
	FESBMPointSource* FECellSBM = nullptr;
	//! The SBM point source
	std::unordered_map<int, FESBMPointSource*> Species;
	//! The species point source
	//std::unordered_map<int, double> Cell_Species;
	//std::unordered_map<int, double> Membrane_Species;
	void InitSBM(FEMesh* mesh);
	void UpdateSBM(FEMesh * mesh);
private:
	vec3d local_pos;
};