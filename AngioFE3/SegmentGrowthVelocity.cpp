#include "SegmentGrowthVelocity.h"
#include "Tip.h"
#include "angio3d.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FESolidDomain.h>
#include "FEAngioMaterialPoint.h"
#include "FEAngio.h"
#include <iostream>

// scales prev (should be just value of 1 if this is placed first) by the constant growth length over time
double SegmentVelocityModifier::ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh)
{
	return prev *= segment_velocity_over_time;
}

bool SegmentVelocityModifier::Init()
{
	return true;
}

void SegmentVelocityModifier::UpdateScale() 
{

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

void SegmentVelocityDensityScaleModifier::UpdateScale() 
{

}

BEGIN_PARAMETER_LIST(SegmentVelocityDensityScaleModifier, SegmentGrowthVelocity)
ADD_PARAMETER(m_density_scale_factor, FE_PARAM_VEC3D, "density_scale_factor");
END_PARAMETER_LIST();

// SL: Added so that growth is scaled only by the referential density
double SegmentVelocityRefDensityScaleModifier::ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh)
{
	std::vector<double> density_at_integration_points;

	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		FEElasticMaterialPoint* elastic_mp = gauss_point->ExtractData<FEElasticMaterialPoint>();

		density_at_integration_points.push_back(angio_mp->ref_ecm_density);
	}

	double density_at_point = interpolation_prop->Interpolate(angio_element->_elem, density_at_integration_points, natural_coords, mesh);
	double density_scale = m_density_scale_factor.x + m_density_scale_factor.y * exp(-m_density_scale_factor.z * density_at_point);

	return density_scale * prev;
}

bool SegmentVelocityRefDensityScaleModifier::Init()
{
	return true;
}

void SegmentVelocityRefDensityScaleModifier::UpdateScale() 
{

}

BEGIN_PARAMETER_LIST(SegmentVelocityRefDensityScaleModifier, SegmentGrowthVelocity)
ADD_PARAMETER(m_density_scale_factor, FE_PARAM_VEC3D, "density_scale_factor");
END_PARAMETER_LIST();

double SigmoidSegmentVelocity::ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh)
{
	return scale*prev;
}

bool SigmoidSegmentVelocity::Init()
{
	return true;
}

void SigmoidSegmentVelocity::UpdateScale()
{
	double time = GetFEModel()->GetTime().currentTime;
	//double shift = GetFEModel()->GetGlobalConstant("max_angio_dt");
	//c = c - shift;
	// account for shift in curve due to the initial min angio dt
	double e_val = -((time - c) / b);
	// derivative of the sigmoid equation
	scale = (a*exp(e_val)) / (b*pow((1 + exp(e_val)), 2));
}

BEGIN_PARAMETER_LIST(SigmoidSegmentVelocity, SegmentGrowthVelocity)
ADD_PARAMETER(scale, FE_PARAM_DOUBLE, "scale");
ADD_PARAMETER(a, FE_PARAM_DOUBLE, "a");
ADD_PARAMETER(b, FE_PARAM_DOUBLE, "b");
ADD_PARAMETER(c, FE_PARAM_DOUBLE, "c");
END_PARAMETER_LIST();

double GompertzSegmentVelocity::ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh)
{
	return scale * prev;
}

bool GompertzSegmentVelocity::Init()
{
	return true;
}

void GompertzSegmentVelocity::UpdateScale()
{
	double time = GetFEModel()->GetTime().currentTime;
	// derivative of Gompertz eqn
	scale = a * b * c * exp(-b * exp(-c * (time - d)))*exp(-c * (time - d));
}

BEGIN_PARAMETER_LIST(GompertzSegmentVelocity, SegmentGrowthVelocity)
ADD_PARAMETER(scale, FE_PARAM_DOUBLE, "scale");
ADD_PARAMETER(a, FE_PARAM_DOUBLE, "a");
ADD_PARAMETER(b, FE_PARAM_DOUBLE, "b");
ADD_PARAMETER(c, FE_PARAM_DOUBLE, "c");
ADD_PARAMETER(d, FE_PARAM_DOUBLE, "d");
END_PARAMETER_LIST();

double SegmentGrowthVelocityManager::ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_elem, FEMesh* mesh)
{
	for (int i = 0; i < seg_vel_modifiers.size(); i++)
	{
		//seg_vel_modifiers[i]->UpdateScale();
		prev = seg_vel_modifiers[i]->ApplyModifiers(prev, natural_coords, angio_elem, mesh);
	}
	return prev;
}

