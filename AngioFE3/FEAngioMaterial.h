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
#include "CultureParameters.h"
#include "SegmentGrowthModifiers.h"
#include "AngioStressPolicy.h"

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

	void GrowSegments(AngioElement * angio_elem, double end_time, int buffer_index);

	void GrowthInElement(double end_time, Tip * active_tip, int source_index, int buffer_index, double override_grow_length=10, bool grow_len_overrride=false);

	//should be const and threadsafe
	void PostGrowthUpdate(AngioElement* angio_elem, double end_time, int buffer_index);

	void Cleanup(AngioElement* angio_elem, double end_time, int buffer_index);

	void PrepBuffers(AngioElement* angio_elem, double end_time, int buffer_index);

	bool SeedFragments(std::vector<AngioElement *>& angio_elements);

	vec3d ApplySegmentDirectionModifiers(Tip * tip, double dt);

	void UpdateGDMs();

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

	double GetSegmentVelocity(AngioElement * angio_element, vec3d local_pos);

	// Calculate spatial elasticity tangent
	tens4ds Tangent(FEMaterialPoint& mp) override;

	//! create material point data for this material
	FEMaterialPoint* CreateMaterialPointData() override;

	//! Set the local coordinate system for a material point (overridden from FEMaterial)
	void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp) override;

	double StrainEnergyDensity(FEMaterialPoint& mp) override;
	FEAngio*	m_pangio = nullptr;
	FEPropertyT<AngioStressPolicy> angio_stress_policy;
private:
	DECLARE_PARAMETER_LIST();

	FEPropertyT<FESolidMaterial> matrix_material;
	
	FEPropertyT<PositionDependentDirectionManager> pdd_manager;
	FEPropertyT<PreviousSegmentContributionManager> psc_manager;
	FEPropertyT<ContributionMixManager> cm_manager;
	FEPropertyT<SegmentGrowthVelocityManager> velocity_manager;
	FEPropertyT<CommonAngioProperties> common_properties;
	CultureParameters m_cultureParams;
	
};