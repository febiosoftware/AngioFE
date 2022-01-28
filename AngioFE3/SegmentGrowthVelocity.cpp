#include "SegmentGrowthVelocity.h"
#include "Tip.h"
#include "angio3d.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FESolidDomain.h>
#include "FEAngioMaterialPoint.h"
#include "FEAngio.h"
#include <iostream>
#include <FECore/mathalg.h>

struct EigComp {
	bool operator()(const std::pair<double, vec3d>& x, std::pair<double, vec3d>& y) const {
		if ((x.first) > (y.first)) return true;
		else return false;
	}
};

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

BEGIN_FECORE_CLASS(SegmentVelocityModifier, SegmentGrowthVelocity)
ADD_PARAMETER(segment_velocity_over_time, "segment_velocity_over_time");
END_FECORE_CLASS();

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

BEGIN_FECORE_CLASS(SegmentVelocityDensityScaleModifier, SegmentGrowthVelocity)
ADD_PARAMETER(m_density_scale_factor, "density_scale_factor");
END_FECORE_CLASS();

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

BEGIN_FECORE_CLASS(SegmentVelocityRefDensityScaleModifier, SegmentGrowthVelocity)
ADD_PARAMETER(m_density_scale_factor, "density_scale_factor");
END_FECORE_CLASS();

double SegmentVelocityDensityFAScaleModifier::ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh)
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

	// vector containing the SPD for each gauss point in the element
	std::vector<mat3ds> SPDs_gausspts;

	// get each gauss point's SPD
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		// get the angio point
		FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		SPDs_gausspts.push_back(angio_mp->angioSPD);
	}

	// array containing the SPD for each node in the element
	mat3ds SPDs_nodes[FESolidElement::MAX_NODES];
	// array for the shape function values
	double H[FESolidElement::MAX_NODES];
	// project the spds from integration points to the nodes
	angio_element->_elem->project_to_nodes(&SPDs_gausspts[0], SPDs_nodes);
	// determine shape function value for the local position
	angio_element->_elem->shape_fnc(H, natural_coords.x, natural_coords.y, natural_coords.z);
	// Get the interpolated SPD from the shape function-weighted Average Structure Tensor
	mat3ds SPD_int = weightedAverageStructureTensor(SPDs_nodes, H, angio_element->_elem->Nodes());
	// get the vectors of the principal directions and sort in descending order
	// the FA can be calculated as std/rms of the ODF. We can assume each direction is one sample and perform the calculation on the eigenvalues.
	double d[3]; vec3d r[3];
	SPD_int.eigen2(d, r);
	std::vector<std::pair<double, vec3d>> v = { {d[0],r[0]}, {d[1],r[1]}, {d[2],r[2]} };
	// Sort based on absolute value of eigenvalues
	std::sort(v.begin(), v.end(), EigComp());
	double angioFA_int = 1 - (v[1].first / v[0].first);

	double L = m_rFA_a / ((1 + pow(((density_at_point - m_rFA_r0) / m_rFA_b), 2)) * (1 + pow(((angioFA_int - m_rFA_f0) / m_rFA_c), 2))) + m_rFA_d;

	return L * prev;
}

bool SegmentVelocityDensityFAScaleModifier::Init()
{
	return true;
}

void SegmentVelocityDensityFAScaleModifier::UpdateScale()
{

}

BEGIN_FECORE_CLASS(SegmentVelocityDensityFAScaleModifier, SegmentGrowthVelocity)
ADD_PARAMETER(m_rFA_a, "rFA_a");
ADD_PARAMETER(m_rFA_b, "rFA_b");
ADD_PARAMETER(m_rFA_c, "rFA_c");
ADD_PARAMETER(m_rFA_d, "rFA_d");
ADD_PARAMETER(m_rFA_r0, "rFA_r0");
ADD_PARAMETER(m_rFA_f0, "rFA_f0");
END_FECORE_CLASS();

double SegmentVelocity3PModifier::ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh)
{
	mat3ds E;
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		FEElasticMaterialPoint* elastic_mp = gauss_point->ExtractData<FEElasticMaterialPoint>();
		E += elastic_mp->Strain();
	}
	E = E / (angio_element->_elem->GaussPoints());
	// get eigenvalues and eigenvectors
	double d[3]; vec3d r[3];
	E.eigen(d, r);
	// Sort based on absolute value of eigenvalues
	double c_strain = std::min(d[0], std::min(d[1], d[2]));
	if (c_strain < threshold) { 
		return scale * prev; 
	}
	else {
		return prev;
	}
}

bool SegmentVelocity3PModifier::Init()
{
	return true;
}

void SegmentVelocity3PModifier::UpdateScale()
{

}

BEGIN_FECORE_CLASS(SegmentVelocity3PModifier, SegmentGrowthVelocity)
ADD_PARAMETER(scale, "scale");
ADD_PARAMETER(threshold, "threshold");
END_FECORE_CLASS();

double SegmentVelocityFAModifier::ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh)
{

	// vector containing the SPD for each gauss point in the element
	std::vector<mat3ds> SPDs_gausspts;

	// get each gauss point's SPD
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		// get the angio point
		FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		SPDs_gausspts.push_back(angio_mp->angioSPD);
	}
	
	
	// array containing the SPD for each node in the element
	mat3ds SPDs_nodes[FESolidElement::MAX_NODES];
	// array for the shape function values
	double H[FESolidElement::MAX_NODES];
	// project the spds from integration points to the nodes
	angio_element->_elem->project_to_nodes(&SPDs_gausspts[0], SPDs_nodes);
	// determine shape function value for the local position
	angio_element->_elem->shape_fnc(H, natural_coords.x, natural_coords.y, natural_coords.z);
	// Get the interpolated SPD from the shape function-weighted Average Structure Tensor
	mat3ds SPD_int = weightedAverageStructureTensor(SPDs_nodes, H, angio_element->_elem->Nodes());
	// get the vectors of the principal directions and sort in descending order
	// the FA can be calculated as std/rms of the ODF. We can assume each direction is one sample and perform the calculation on the eigenvalues.
	double d[3]; vec3d v[3];
	SPD_int.eigen2(d, v);
	double sum = d[0] + d[1] + d[2];
	double mean = sum / 3.0;
	double sq_sum = (d[0] - mean) * (d[0] - mean) + (d[1] - mean) * (d[1] - mean) + (d[2] - mean) * (d[2] - mean);
	double stdev = sqrt(sq_sum / 2.0);
	double rms = sqrt((d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) / 3.0);
	double angioFA_int = stdev / rms;
	
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
	double dens_point = std::min(std::max(0.5, density_at_point), 1.0);
	scale = 1.0 + (-0.25*density_at_point+1.75)*(density_scale * (1.0 / (1.0 + exp((-1.0 * angioFA_int + 0.47) / 0.005))));
	return prev*scale;
}

bool SegmentVelocityFAModifier::Init()
{
	return true;
}

void SegmentVelocityFAModifier::UpdateScale()
{

}

BEGIN_FECORE_CLASS(SegmentVelocityFAModifier, SegmentGrowthVelocity)
ADD_PARAMETER(scale, "scale");
ADD_PARAMETER(m_density_scale_factor, "density_scale_factor");
END_FECORE_CLASS();

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
	// account for shift in curve due to the initial min angio dt
	double e_val = -((time - c) / b);
	// derivative of the sigmoid equation
	scale = (a*exp(e_val)) / (b*pow((1 + exp(e_val)), 2));
}

BEGIN_FECORE_CLASS(SigmoidSegmentVelocity, SegmentGrowthVelocity)
ADD_PARAMETER(scale, "scale");
ADD_PARAMETER(a, "a");
ADD_PARAMETER(b, "b");
ADD_PARAMETER(c, "c");
END_FECORE_CLASS();

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

BEGIN_FECORE_CLASS(GompertzSegmentVelocity, SegmentGrowthVelocity)
ADD_PARAMETER(scale, "scale");
ADD_PARAMETER(a, "a");
ADD_PARAMETER(b, "b");
ADD_PARAMETER(c, "c");
ADD_PARAMETER(d, "d");
END_FECORE_CLASS();

double SegmentGrowthVelocityManager::ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_elem, FEMesh* mesh)
{
	for (int i = 0; i < seg_vel_modifiers.size(); i++)
	{
		prev = seg_vel_modifiers[i]->ApplyModifiers(prev, natural_coords, angio_elem, mesh);
	}
	return prev;
}

