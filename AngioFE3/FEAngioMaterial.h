#pragma once
#include <FEBioMech/FEElasticMaterial.h>
#include <FEBioMech/FEElasticFiberMaterial.h>
#include <FECore/FEDataArray.h>
#include <FECore/FESurface.h>
#include <FECore/FENormalProjection.h>
#include "FEAngio.h"
#include "FEProbabilityDistribution.h"
#include "KDTree/kdtree.h"
#include "FiberManager.h"
#include "FEAngioMaterialPoint.h"
#include <FEBioMix/FEBiphasic.h>
#include "ECMInitializer.h"
#include "CommonAngioProperites.h"
#include <FEBioMix/FEMultiphasic.h>
#include "ContributionMix.h"
#include "PositionDependentDirection.h"
#include "PreviousSegmentContribution.h"
#include "SegmentGrowthVelocity.h"
#include "AngioStressPolicy.h"
#include "InitialModifiers.h"
#include "NodeDataInterpolation.h"
#include "TipDoping.h"

class AngioElement;

//-----------------------------------------------------------------------------
// Class implementing a stress induced by a non-local point force
class FEAngioMaterial : public FEElasticMaterial
{
public:
	

	explicit FEAngioMaterial(FEModel* pfem);
	virtual ~FEAngioMaterial();

	friend class Fileout;
	friend class FEPlotMatrixStress;
	friend class FEPlotVesselStress;
	friend class FEPlotMatrixTangent;

	//begin functions from FEAngioMaterialBase

	// Calculate the active Angio stress
	mat3ds AngioStress(FEAngioMaterialPoint& mp);

	//not threadsafe
	void SetSeeds(AngioElement* angio_elem);


	//should be const and threadsafe
	double GetMin_dt(AngioElement* angio_elem, FEMesh* mesh);

	void ProtoGrowSegments(AngioElement * angio_elem, double end_time, int buffer_index, double min_scale_factor, double bounds_tolerance, double min_angle);

	void GrowSegments(AngioElement * angio_elem, double end_time, int buffer_index, double min_scale_factor, double bounds_tolerance, double min_angle);

	void GrowthInElement(double end_time, Tip * active_tip, int source_index, int buffer_index, double min_scale_factor, double bounds_tolerance, double min_angle);

	void ProtoGrowthInElement(double end_time, Tip * active_tip, int source_index, int buffer_index, double min_scale_factor, double bounds_tolerance, double min_angle);

	//returns the index of the element in which the tip will grow next, this resolves when a tip can grow in multiple elements
	int SelectNextTip(std::vector<AngioElement*> & possible_locations, std::vector<vec3d> & possible_local_coordinates, Tip* tip, double dt, int buffer, FEMesh* mesh, double min_scale_factor, double min_angle);

	//should be const and threadsafe
	void PostGrowthUpdate(AngioElement* angio_elem, double end_time, double final_time, double min_scale_factor, int buffer_index, FEMesh* mesh, FEAngio* feangio);

	void ProtoPostGrowthUpdate(AngioElement* angio_elem, double end_time, double min_scale_factor, int buffer_index, FEMesh* mesh, FEAngio* feangio);

	void Cleanup(AngioElement* angio_elem, double end_time, int buffer_index);

	void PrepBuffers(AngioElement* angio_elem, double end_time, int buffer_index);

	bool SeedFragments(std::vector<AngioElement *>& angio_elements, FEMesh* mesh);

	vec3d ApplySegmentDirectionModifiers(Tip * tip, double dt);

	double FindDensityScale(FEAngioMaterialPoint * mp);

	void SetFEAngio(FEAngio * ctl)
	{
		m_pangio = ctl;
	}
	
	FEMaterial * GetMatrixMaterial()  { return matrix_material; }

	CommonAngioProperties * GetCommonAngioProperties()  { return common_properties; };

	int GetID_ang() { 
		FECoreBase  * base = GetParent();
		if(base)
		{
			return base->GetID();
		}
		return FEElasticMaterial::GetID(); };

	FEMaterial * GetMaterial() { return dynamic_cast<FEMaterial*>(this); }
	//begin functions from FEMaterial

	// material initialization
	bool Init() override;

	// Calculate Cauchy-stress
	mat3ds Stress(FEMaterialPoint& mp) override;

	double GetSegmentVelocity(AngioElement * angio_element, vec3d local_pos, FEMesh* mesh);

	double GetInitialVelocity();

	// Calculate spatial elasticity tangent
	tens4ds Tangent(FEMaterialPoint& mp) override;

	//! create material point data for this material
	FEMaterialPoint* CreateMaterialPointData() override;

	//! Set the local coordinate system for a material point (overridden from FEMaterial)
	void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp) override;

	double StrainEnergyDensity(FEMaterialPoint& mp) override;
	FEAngio*	m_pangio = nullptr;
	FEPropertyT<AngioStressPolicy> angio_stress_policy;
	FEPropertyT<InitialModifierManager> im_manager;
	FEPropertyT<BranchPolicy> branch_policy;
	FEPropertyT<BranchPolicy> proto_branch_policy;
	FEPropertyT<NodeDataInterpolationManager> nodedata_interpolation_manager;
	FEPropertyT<TipDopingManager> tip_doping_manager;
	double vessel_radius = 7.0;
private:
	DECLARE_PARAMETER_LIST();

	double initial_segment_velocity = 40.0;
	double dt_safety_multiplier = 1.0;

	FEPropertyT<FESolidMaterial> matrix_material;
	
	FEPropertyT<PositionDependentDirectionManager> pdd_manager;
	FEPropertyT<PreviousSegmentContributionManager> psc_manager;
	FEPropertyT<ContributionMixManager> cm_manager;
	FEPropertyT<SegmentGrowthVelocityManager> velocity_manager;
	FEPropertyT<CommonAngioProperties> common_properties;
	
};