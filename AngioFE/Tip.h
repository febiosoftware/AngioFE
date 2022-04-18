#pragma once
#include <FECore/vec3d.h>
#include "AngioElement.h"
#include <FEBioMix/FESBMPointSource.h>
#include "FEProbabilityDistribution.h"
#include "FECell.h"

class Segment;
class FEAngio;

//! a location where growth has occured, or can occur on the next growth step/substep.
class Tip
{
public:
	//! constructor
	Tip(){}
	//! used to copy another tip
	Tip(Tip * other, FEMesh * mesh);
	//! sets the local positions clamps the values to -1 to 1
	void SetLocalPosition(vec3d pos, FEMesh * mesh);
	//! returns the local position
	vec3d GetLocalPosition() const;
	//! The element that contains the local_pos coordinates
	AngioElement * angio_element = nullptr;
	//! Time at which the tip occurs
	double time = 0.0;
	//! evaluation time
	double etime = 0.0;
	//! segment initial time
	double seg_time = 0.0;
	//! Velocity at which the tip grew
	double growth_velocity;
	//! The element where growth originated from
	AngioElement* face = nullptr;
	//! Segment that contains this this tip. May be nullptr for no parent segment 
	Segment * parent = nullptr;
	//! Will be initialized to values greater than or equal to zero, unique values per vessel
	int initial_fragment_id = -1;
	//! identifies whether or not a given tip is the base of a branch
	bool is_branch = false;
	//! Use the direction member or not
	bool use_direction = false;
	//! Only used if the tip is a branch
	vec3d direction;
	//! Use direction or the direction of the parent segment
	vec3d GetDirection(FEMesh* mesh) const;
	//! Return the global position of the tip
	vec3d GetPosition(FEMesh * mesh) const;
	//! returns the global position of the tip in the reference frame
	vec3d GetRefPosition(FEMesh * mesh) const;
	//! Prints information about the tip to the console
	void PrintTipInfo(FEMesh *mesh) const;
	//! Prints information about the tip to the console
	void PrintTipInfo(FEMesh *mesh, std::string title) const;
	FECell * TipCell = nullptr;
	void SetProtoGrowthLength(FEProbabilityDistribution* dist);
	void SetProtoGrowthLength(Tip* tip);
	void SetProtoGrowthLength(double initial_length);
	double GetProtoGrowthLength();
private:
	vec3d local_pos;
	//! Proto growth velocity
	double proto_growth_length;
};

