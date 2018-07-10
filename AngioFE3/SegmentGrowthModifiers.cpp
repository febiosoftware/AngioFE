#include "SegmentGrowthModifiers.h"
#include "Tip.h"
#include "angio3d.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FESolidDomain.h>
#include "FEAngioMaterialPoint.h"
#include "FEAngio.h"

vec3d PreviousSegmentPSC::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh)
{
	return prev_direction;
}

void FiberPDD::Update(FEMesh * mesh, FEAngio* angio)
{

}

double PSCPDDContributionMix::ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh)
{
	return psc_weight*prev;
}

void PSCPDDContributionMix::Update(FEMesh * mesh)
{
	// TODO: Implement.
}

BEGIN_PARAMETER_LIST(PSCPDDContributionMix, ContributionMix)
ADD_PARAMETER(psc_weight, FE_PARAM_DOUBLE, "psc_weight");
END_PARAMETER_LIST();

double SegmentVelocityModifier::ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh)
{
	return prev *= segment_velocity_over_time;
}

bool SegmentVelocityModifier::Init()
{
	return true;
}

BEGIN_PARAMETER_LIST(SegmentVelocityModifier, SegmentGrowthVelocity)
ADD_PARAMETER(segment_velocity_over_time, FE_PARAM_DOUBLE, "segment_velocity_over_time");
END_PARAMETER_LIST();

double SegmentVelocityDensityScaleModifier::ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh)
{
	std::vector<double> density_at_integration_points;

	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		FEElasticMaterialPoint* elastic_mp = gauss_point->ExtractData<FEElasticMaterialPoint>();

		density_at_integration_points.push_back(angio_mp->ref_ecm_density * (1.0 / elastic_mp->m_J));
	}

	double density_at_point = interpolation_prop->Interpolate(angio_element->_elem, density_at_integration_points, natural_coords, mesh);
	double density_scale = m_density_scale_factor.x + m_density_scale_factor.y * exp(-m_density_scale_factor.z * density_at_point);

	return density_scale * prev;
}

bool SegmentVelocityDensityScaleModifier::Init()
{
	return true;
}

BEGIN_PARAMETER_LIST(SegmentVelocityDensityScaleModifier, SegmentGrowthVelocity)
ADD_PARAMETER(m_density_scale_factor, FE_PARAM_VEC3D, "density_scale_factor");
END_PARAMETER_LIST();

vec3d FiberPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, double& alpha, FEMesh* mesh, FEAngio* pangio)
{
	std::vector<quatd> gauss_data;
	std::vector<quatd> nodal_data;
	for(int i=0; i< angio_element->_elem->GaussPoints();i++)
	{
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		vec3d axis(1, 0, 0);
		axis = emp->m_F * emp->m_Q * axis;
		gauss_data.push_back({ axis });
	}
	
	quatd rv = interpolation_prop->Interpolate(angio_element->_elem, gauss_data, local_pos, mesh);
	vec3d fiber_direction = rv.GetVector();
	return mix(prev, fiber_direction, contribution);
}

BEGIN_PARAMETER_LIST(ECMDensityGradientPDD, PositionDependentDirection)
ADD_PARAMETER(threshold, FE_PARAM_DOUBLE, "threshold");
END_PARAMETER_LIST();

vec3d ECMDensityGradientPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, double& alpha, FEMesh* mesh, FEAngio* pangio)
{
	std::vector<double> density_at_integration_points;
	
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		FEElasticMaterialPoint* elastic_mp = gauss_point->ExtractData<FEElasticMaterialPoint>();

		density_at_integration_points.push_back(angio_mp->ref_ecm_density * (1.0 / elastic_mp->m_J));
	}
	vec3d grad = pangio->gradient(angio_element->_elem, density_at_integration_points, local_pos);
	double gradnorm = grad.norm();
	if(gradnorm > threshold)
	{
		vec3d currentDirectionGradientPlane = grad ^ prev;
		currentDirectionGradientPlane.unit();
		vec3d perpendicularToGradient = currentDirectionGradientPlane ^ grad;
		perpendicularToGradient.unit();
		alpha = 1.0;
		return mix(prev, perpendicularToGradient, contribution);
	}
	return prev;
}

vec3d AnastomosisPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, double& alpha, FEMesh* mesh, FEAngio* pangio)
{
	return vec3d(0, 0, 0);
}

BEGIN_PARAMETER_LIST(AnastomosisPDD, PositionDependentDirection)
ADD_PARAMETER(anastomosis_radius, FE_PARAM_DOUBLE, "anastomosis_radius");
ADD_PARAMETER(contribution, FE_PARAM_DOUBLE, "contribution");
END_PARAMETER_LIST();

// Managers
vec3d PreviousSegmentContributionManager::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, vec3d prev_direction, FEMesh* mesh)
{
	for(int i=0; i < psc_modifiers.size();i++)
	{
		prev = psc_modifiers[i]->ApplyModifiers(prev, angio_element,local_pos, prev_direction, mesh);
	}
	return prev;
}


vec3d PositionDependentDirectionManager::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, double& alpha, FEMesh* mesh, FEAngio* pangio)
{
	for (int i = 0; i < pdd_modifiers.size(); i++)
	{
		prev = pdd_modifiers[i]->ApplyModifiers(prev, angio_element, local_pos, alpha, mesh, pangio);
	}
	return prev;
}

double ContributionMixManager::ApplyModifiers(double prev, AngioElement* angio_element, vec3d local_pos, FEMesh* mesh)
{
	for (int i = 0; i < cm_modifiers.size(); i++)
	{
		prev = cm_modifiers[i]->ApplyModifiers(prev, angio_element,local_pos, mesh);
	}
	return prev;
}

double SegmentGrowthVelocityManager::ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_elem, FEMesh* mesh)
{
	for (int i = 0; i < seg_vel_modifiers.size(); i++)
	{
		prev = seg_vel_modifiers[i]->ApplyModifiers(prev, natural_coords, angio_elem, mesh);
	}
	return prev;
}