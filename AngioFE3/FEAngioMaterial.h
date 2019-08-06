#pragma once
#include <FEBioMech/FEElasticMaterial.h>
#include <FEBioMech/FEElasticFiberMaterial.h>
#include <FECore/FEDataArray.h>
#include <FECore/FESurface.h>
#include <FECore/FENormalProjection.h>
#include "FEAngio.h"
#include "FEProbabilityDistribution.h"
#include "FiberManager.h"
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
#include "Tip.h"

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

	//! vascular growth that occours before time =0
	void ProtoGrowSegments(AngioElement * angio_elem, double end_time, int buffer_index, double min_scale_factor, double bounds_tolerance, double min_angle);
	
	//! vascular growth that occours after time =0
	void GrowSegments(AngioElement * angio_elem, double end_time, int buffer_index, double min_scale_factor, double bounds_tolerance, double min_angle);

	//! performs vascular growth within an element
	void GrowthInElement(double end_time, Tip * active_tip, int source_index, int buffer_index, double min_scale_factor, double bounds_tolerance, double min_angle);

	//! performs vascular growth within an element before time = 0
	void ProtoGrowthInElement(double end_time, Tip * active_tip, int source_index, int buffer_index, double min_scale_factor, double bounds_tolerance, double min_angle);

	//! returns the index of the element in which the tip will grow next, this resolves when a tip can grow in multiple elements
	int SelectNextTip(std::vector<AngioElement*> & possible_locations, std::vector<vec3d> & possible_local_coordinates, Tip* tip, double dt, int buffer, FEMesh* mesh, double min_scale_factor, double min_angle);

	//should be const and threadsafe
	//! sets buffers to the correct state and perfom branching as needed, run after each growth substep
	void PostGrowthUpdate(AngioElement* angio_elem, double end_time, double final_time, double min_scale_factor, int buffer_index, FEMesh* mesh, FEAngio* feangio);

	//! sets buffers to the correct state and perfom branching as needed
	void ProtoPostGrowthUpdate(AngioElement* angio_elem, double end_time, double min_scale_factor, int buffer_index, FEMesh* mesh, FEAngio* feangio);

	//! clean up the buffers at the end of a growth step
	void Cleanup(AngioElement* angio_elem, double end_time, int buffer_index);

	//! prep the buffers before each growth substep
	static void PrepBuffers(AngioElement* angio_elem, double end_time, int buffer_index);

	//! seed the fragments within the given elements
	bool SeedFragments(std::vector<AngioElement *>& angio_elements, FEMesh* mesh);

	//! returns the density scale at the integration point
	double FindDensityScale(FEAngioMaterialPoint * mp);

	//! sets the internal pointer to the class that runs everything
	void SetFEAngio(FEAngio * ctl)
	{
		m_pangio = ctl;
	}
	
	//! returns the matrix material
	FEMaterial * GetMatrixMaterial()  { return matrix_material; }

	//! returns a pointer to the common angio properites
	CommonAngioProperties * GetCommonAngioProperties()  { return common_properties; };

	//! gets the febio id of the angio material
	int GetID_ang() { 
		FECoreBase  * base = GetParent();
		if(base)
		{
			return base->GetID();
		}
		return FEElasticMaterial::GetID(); };

	//begin functions from FEMaterial

	//! material initialization
	bool Init() override;

	//! Calculate Cauchy-stress
	mat3ds Stress(FEMaterialPoint& mp) override;

	//! returns the segment velocity at a given location
	double GetSegmentVelocity(AngioElement * angio_element, vec3d local_pos, FEMesh* mesh);

	//! returns the inital segment velocity
	double GetInitialVelocity() const;

	//! Calculate spatial elasticity tangent
	tens4ds Tangent(FEMaterialPoint& mp) override;

	//! create material point data for this material
	FEMaterialPoint* CreateMaterialPointData() override;

	//! Set the local coordinate system for a material point (overridden from FEMaterial)
	void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp) override;

	//! calculates the strain energy density. Effects convergence rate
	double StrainEnergyDensity(FEMaterialPoint& mp) override;
	//! pointer to controlling class
	FEAngio*	m_pangio = nullptr;
	//! policy to calculate stress
	FEPropertyT<AngioStressPolicy> angio_stress_policy;
	//! initialization of elements
	FEPropertyT<InitialModifierManager> im_manager;
	//! branching policy
	FEPropertyT<BranchPolicy> branch_policy;
	//! branching policy for time before t=0
	FEPropertyT<BranchPolicy> proto_branch_policy;
	//! interpolation and initialization of per node properties
	FEPropertyT<NodeDataInterpolationManager> nodedata_interpolation_manager;
	//! radius of vessels used to calculate mixture between matrix and vessel materials
	double vessel_radius = 7.0;
	FEPropertyT<FEMixMethod> mix_method;
	FEPropertyT<SegmentGrowthVelocityManager> velocity_manager;
	//FEPropertyT<TipSpeciesManager> tip_species_manager;
private:
	DECLARE_PARAMETER_LIST();

	double initial_segment_velocity = 40.0;
	double dt_safety_multiplier = 1.0;

	FEPropertyT<FESolidMaterial> matrix_material;
	FEPropertyT<PositionDependentDirectionManager> pdd_manager;
	FEPropertyT<PreviousSegmentContributionManager> psc_manager;
	//FEPropertyT<TipSpecies> tip_species;
	FEPropertyT<ContributionMixManager> cm_manager;
//	FEPropertyT<SegmentGrowthVelocityManager> velocity_manager;
	FEPropertyT<CommonAngioProperties> common_properties;
	
};