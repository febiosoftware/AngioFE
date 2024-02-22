#pragma once
#include <FEBioMech/FEElasticMaterial.h>
#include <FEBioMech/FEElasticFiberMaterial.h>
#include <FECore/FEDataArray.h>
#include <FECore/FESurface.h>
#include <FECore/FENormalProjection.h>
#include "FEAngio.h"
#include "FEProbabilityDistribution.h"
#include "FEAngioMaterialPoint.h"
#include <FEBioMix/FEBiphasic.h>
#include "CommonAngioProperites.h"
#include <FEBioMix/FEMultiphasic.h>
#include "ContributionMix.h"
#include "PositionDependentDirection.h"
#include "PreviousSegmentContribution.h"
#include "SegmentGrowthVelocity.h"
#include "AngioStressPolicy.h"
#include "InitialModifiers.h"
#include "NodeDataInterpolation.h"

class AngioElement;

//! the material that allows vascular growth within it
class FEAngioMaterial : public FEElasticMaterial
{
public:
	
	//! constructor
	explicit FEAngioMaterial(FEModel* pfem);
	//! destructor
	virtual ~FEAngioMaterial();

	friend class Fileout;
	friend class FEPlotMatrixStress;
	friend class FEPlotVesselStress;
	friend class FEPlotMatrixTangent;

	//! Calculate the stress from the stress policy(ussually tips)
	mat3ds AngioStress(FEAngioMaterialPoint& mp);

	//not threadsafe
	//! set the seed for a given element
	void SetSeeds(AngioElement* angio_elem);

	//should be const and threadsafe
	//! the smallest safe step an elment can take
	double GetMin_dt(AngioElement* angio_elem, FEMesh* mesh);

	//! vascular growth that occurs before time =0
	void ProtoGrowSegments(AngioElement * angio_elem, double end_time, int buffer_index);
	
	//! vascular growth that occurs after time =0
	void GrowSegments(AngioElement * angio_elem, double end_time, int buffer_index);

	//! performs vascular growth within an element
	void GrowthInElement(double end_time, Tip* active_tip, int source_index, int buffer_index);

	//! performs vascular growth within an element before time = 0
	void ProtoGrowthInElement(double end_time, Tip* active_tip, int source_index, int buffer_index);

	//! returns the index of the element in which the tip will grow next, this resolves when a tip can grow in multiple elements
	int SelectNextTip(std::vector<AngioElement*>& possible_locations, std::vector<vec3d>& possible_local_coordinates, Tip* tip, double dt, int buffer);

	//should be const and threadsafe
	//! sets buffers to the correct state and perfom branching as needed, run after each growth substep
	void PostGrowthUpdate(AngioElement* angio_elem, int buffer_index);

	//! sets buffers to the correct state and perfom branching as needed
	void ProtoPostGrowthUpdate(AngioElement* angio_elem, int buffer_index);

	//! clean up the buffers at the end of a growth step
	void Cleanup(AngioElement* angio_elem, double end_time, int buffer_index);

	//! prep the buffers before each growth substep
	static void PrepBuffers(AngioElement* angio_elem, double end_time, int buffer_index);

	//! seed the fragments within the given elements
	bool SeedFragments(std::vector<AngioElement *>& angio_elements, FEMesh* mesh);

	//! returns the density scale at the integration point
	double FindDensityScale(FEAngioMaterialPoint * mp);

	//! sets the internal pointer to the class that runs everything
	void SetFEAngio(FEAngio * ctl) {m_pangio = ctl; }
	
	//! returns the matrix material
	FEMaterial * GetMatrixMaterial()  { return matrix_material; }

	//! returns a pointer to the common angio properites
	CommonAngioProperties * GetCommonAngioProperties()  { return common_properties; };

	//! gets the febio id of the angio material
	int GetID_ang() 
	{ 
		FECoreBase* base = GetParent();
		if (base)
			return base->GetID();

		return FEElasticMaterial::GetID(); 
	};

	//! functions from FEMaterial
	//! material initialization
	bool Init() override;
	//! Calculate Cauchy-stress
	mat3ds Stress(FEMaterialPoint& mp) override;
	//! Calculate spatial elasticity tangent
	tens4ds Tangent(FEMaterialPoint& mp) override;
	//! create material point data for this material
	FEMaterialPointData* CreateMaterialPointData() override;
	//! calculates the strain energy density. Effects convergence rate
	double StrainEnergyDensity(FEMaterialPoint& mp) override;

	//! Angio material functions
	//! //! returns the segment velocity at a given location
	double GetSegmentVelocity(AngioElement* angio_element, vec3d local_pos, double time_shift);
	//! prevent vessels from getting stuck on faces
	vec3d CheckFaceProximity(vec3d pos, vec3d dir);
	//! find the possible places to move to when the tip is stuck on an exterior face/corner
	std::vector<double> CornerPossibleValues(vec3d local_pos, vec3d nat_dir);

	//! Public pointers to policies and managers
	FEAngio* m_pangio = nullptr;
	AngioStressPolicy* angio_stress_policy = nullptr;
	InitialModifierManager* im_manager = nullptr;
	BranchPolicy* branch_policy = nullptr;
	BranchPolicy* proto_branch_policy = nullptr;
	NodeDataInterpolationManager* nodedata_interpolation_manager = nullptr;
	SegmentGrowthVelocityManager* velocity_manager = nullptr;

	//! Properties
	double vessel_radius = 6.3;
	double thresh_vess_weight = 0.0215;
	FEMixMethod* mix_method = nullptr;
	
private:
	DECLARE_FECORE_CLASS()
	
	//! Private properties
	double initial_segment_velocity = 7.5;
	double dt_safety_multiplier = 1.0;

	//! Private pointers to policies and managers
	FEElasticMaterial* matrix_material = nullptr;
	ProtoPreviousSegmentContributionManager* proto_psc_manager = nullptr;
	PreviousSegmentContributionManager* psc_manager = nullptr;
	ProtoPositionDependentDirectionManager* proto_pdd_manager = nullptr;
	PositionDependentDirectionManager* pdd_manager = nullptr;
	ProtoContributionMixManager* proto_cm_manager = nullptr;
	ContributionMixManager* cm_manager = nullptr;
	CommonAngioProperties* common_properties = nullptr;
};