#include "AngioStressPolicy.h"
#include "AngioElement.h"
#include "FEAngio.h"
#include "Tip.h"
#include <FECore/FELoadCurve.h>
#include <iostream>

#pragma region FECoreClassDefs
BEGIN_FECORE_CLASS(SigmoidAngioStressPolicy, AngioStressPolicy)
	ADD_PARAMETER(sprout_mag, "sprout_mag");
	ADD_PARAMETER(fan_exponential, "fan_exponential");
	ADD_PARAMETER(sprout_range, "sprout_range");
	ADD_PARAMETER(sprout_radius_multiplier, "sprout_radius_multiplier");
	ADD_PARAMETER(a, "a");
	ADD_PARAMETER(b, "b");
	ADD_PARAMETER(x0, "x0");
	ADD_PARAMETER(y0, "y0");
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(SigmoidDensAngioStressPolicy, AngioStressPolicy)
	ADD_PARAMETER(sprout_mag, "sprout_mag");
	ADD_PARAMETER(fan_exponential, "fan_exponential");
	ADD_PARAMETER(sprout_range, "sprout_range");
	ADD_PARAMETER(sprout_radius_multiplier, "sprout_radius_multiplier");
	ADD_PARAMETER(a, "a");
	ADD_PARAMETER(b, "b");
	ADD_PARAMETER(x0, "x0");
	ADD_PARAMETER(y0, "y0");
	ADD_PARAMETER(scale, "scale");
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(LoadCurveVelAngioStressPolicy, AngioStressPolicy)
	ADD_PARAMETER(sprout_mag, "sprout_mag");
	ADD_PARAMETER(fan_exponential, "fan_exponential");
	ADD_PARAMETER(sprout_range, "sprout_range");
	ADD_PARAMETER(sprout_radius_multiplier, "sprout_radius_multiplier");
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(LoadCurveAngioStressPolicy, AngioStressPolicy)
	ADD_PARAMETER(sprout_mag, "sprout_mag");
	ADD_PARAMETER(fan_exponential, "fan_exponential");
	ADD_PARAMETER(sprout_range, "sprout_range");
	ADD_PARAMETER(sprout_radius_multiplier, "sprout_radius_multiplier");
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(LoadCurveDenAngioStressPolicy, AngioStressPolicy)
	ADD_PARAMETER(sprout_mag, "sprout_mag");
	ADD_PARAMETER(fan_exponential, "fan_exponential");
	ADD_PARAMETER(sprout_range, "sprout_range");
	ADD_PARAMETER(sprout_radius_multiplier, "sprout_radius_multiplier");
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(LoadCurveRefDenAngioStressPolicy, AngioStressPolicy)
	ADD_PARAMETER(sprout_mag, "sprout_mag");
	ADD_PARAMETER(fan_exponential, "fan_exponential");
	ADD_PARAMETER(sprout_range, "sprout_range");
	ADD_PARAMETER(sprout_radius_multiplier, "sprout_radius_multiplier");
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(GrownSegmentsAngioStressPolicy, AngioStressPolicy)
	ADD_PARAMETER(sprout_mag, "sprout_mag");
	ADD_PARAMETER(fan_exponential, "fan_exponential");
	ADD_PARAMETER(sprout_range, "sprout_range");
	ADD_PARAMETER(sprout_radius_multiplier, "sprout_radius_multiplier");
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(GrownSegmentsVelAngioStressPolicy, AngioStressPolicy)
	ADD_PARAMETER(sprout_mag, "sprout_mag");
	ADD_PARAMETER(fan_exponential, "fan_exponential");
	ADD_PARAMETER(sprout_range, "sprout_range");
	ADD_PARAMETER(sprout_radius_multiplier, "sprout_radius_multiplier");
END_FECORE_CLASS()
#pragma endregion FECoreClassDefs

void AngioStressPolicy::UpdateToLoadCurve(const char* param_name, double& value)
{
	//if load curves are used they must use step interpolation
	FEParam* m = FindParameter(ParamString(param_name));
	assert(m);

	FEModel* model = GetFEModel();
	FELoadCurve* mlc = dynamic_cast<FELoadCurve*>(model->GetLoadController(m));
	assert(mlc);
	if (mlc) 
		value = mlc->GetValue(model->GetTime().currentTime);
}

double AngioStressPolicy::GetDensScale(AngioElement* angio_element, Tip* tip, FEMesh* mesh, int is_ref)
{
	std::vector<double> density_at_integration_points;
	auto se = angio_element->_elem;

	//this is really slow and it isn't apparent why.
	int nint = se->GaussPoints();
	for (int i = 0; i < nint; i++)
	{
		FEMaterialPoint* gauss_point = se->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		FEElasticMaterialPoint* elastic_mp = gauss_point->ExtractData<FEElasticMaterialPoint>();
		// determines whether to use the referential or current density
		double m_J = pow(elastic_mp->m_J, is_ref);
		density_at_integration_points.push_back(angio_mp->ref_ecm_density * (1.0 / m_J));
	}

	//Get interpolation method
	FEModel* m_pfem = this->GetFEModel();
	PerElementVI interp(m_pfem);
	double density_at_point = interp.Interpolate(se, density_at_integration_points, tip->GetLocalPosition(), mesh);
	double density_scale = m_density_scale_factor.x + (m_density_scale_factor.y * exp(-m_density_scale_factor.z * density_at_point));
	if (density_scale < 0) 
		return 0; 
	else 
		return density_scale; 
}

bool SigmoidAngioStressPolicy::Init()
{
	return true;
}

void SigmoidAngioStressPolicy::UpdateScale()
{
	double time = GetFEModel()->GetTime().currentTime;
	scale = y0 + a / (1.0 + exp(-(time - x0) / b));
}

void SigmoidAngioStressPolicy::AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh)
{
	std::vector<Tip*> final_active_tips;
	auto se = angio_element->_elem;
	// Get active tips within a radius
	double total_sprout_falloff = sprout_range * sprout_radius_multiplier;
	FEAngio::GetActiveFinalTipsInRadius(angio_element, total_sprout_falloff, pangio, final_active_tips);
	// for each integration point
	int nint = se->GaussPoints();
	for (int i = 0; i < nint; i++)
	{
		// get the material point of the gauss point
		FEMaterialPoint* mp = se->GetMaterialPoint(i);
		// get eh angio material point at the gauss point
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		// get the elastic material point at the gauss point
		FEElasticMaterialPoint* emp = mp->ExtractData<FEElasticMaterialPoint>();
		// assert that the material point and gauss point are in the same point
		assert(mp && angio_mp);
		// zero the angio point
		angio_mp->m_as.zero();
		// get global position of elastic material point?
		vec3d y = mp->m_rt;
		// for each active tip
		int ntips = final_active_tips.size();
		for (int j = 0; j < ntips; j++)
		{
			Tip* tip = final_active_tips[j];
			// get the tip position in the mesh 
			// consider moving this out and calling it less
			vec3d x = tip->GetPosition(mesh); 
			// determine vector from tip to integration point
			vec3d r = y - x;
			double l = r.unit();
			// get the angle between the tip and the sprout direction 
			// consider moving this out and calling it less
			double theta = acos(tip->GetDirection(mesh) * r);

			//sprout s mag replaced with giving correct coeficients for scale
			double p = scale * sprout_mag * pow(cos(theta / 2.0), fan_exponential) * exp(-l / sprout_range);
			angio_mp->m_as += dyad(r) * p;
		}
	}
}

bool SigmoidDensAngioStressPolicy::Init()
{
	return true;
}

void SigmoidDensAngioStressPolicy::UpdateScale()
{
	double time = GetFEModel()->GetTime().currentTime;
	scale = y0 + a / (1.0 + exp(-(time - x0) / b));
}

void SigmoidDensAngioStressPolicy::
	AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh)
{
	std::vector<Tip*> final_active_tips;
	auto se = angio_element->_elem;
	// Get active tips within a radius
	double total_sprout_falloff = sprout_range * sprout_radius_multiplier;
	FEAngio::GetActiveFinalTipsInRadius(angio_element, total_sprout_falloff, pangio, final_active_tips);
	// for each integration point
	int nint = se->GaussPoints();
	for (int i = 0; i < nint; i++)
	{
		// get the material point of the gauss point
		FEMaterialPoint* mp = se->GetMaterialPoint(i);
		// get eh angio material point at the gauss point
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		// get the elastic material point at the gauss point
		FEElasticMaterialPoint* emp = mp->ExtractData<FEElasticMaterialPoint>();
		// get the density scale for the integration point
		double den_scale = angio_element->_angio_mat->FindDensityScale(angio_mp);
		// assert that the material point and gauss point are in the same point
		assert(mp && angio_mp);
		// zero the angio point
		angio_mp->m_as.zero();
		// get global position of elastic material point?
		vec3d y = mp->m_rt;
		this->UpdateScale();
		// for each active tip
		int ntips = final_active_tips.size();
		for (int j = 0; j < ntips; j++)
		{
			Tip* tip = final_active_tips[j];
			// get the tip position in the mesh
			//consider moving this out and calling it less
			vec3d x = tip->GetPosition(mesh);
			// determine vector from tip to integration point
			vec3d r = y - x;
			double l = r.unit();
			// get the angle between the tip and the sprout direction
			//consider moving this out and calling it less
			double theta = acos(tip->GetDirection(mesh) * r);
			//sprout s mag replaced with giving correct coeficients for scale
			double p = den_scale * scale * sprout_mag * pow(cos(theta / 2.0), fan_exponential) * exp(-l / sprout_range);
			angio_mp->m_as += dyad(r) * p;
		}
	}
}

bool LoadCurveVelAngioStressPolicy::Init()
{
	return true;
}

void LoadCurveVelAngioStressPolicy::UpdateScale()
{
	//! Implementation blank
}

void LoadCurveVelAngioStressPolicy::AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh)
{
	std::vector<Tip*> final_active_tips;
	auto se = angio_element->_elem;
	// get all tips within a radius
	double total_sprout_falloff = sprout_range * sprout_radius_multiplier;
	FEAngio::GetActiveFinalTipsInRadius(angio_element, total_sprout_falloff, pangio, final_active_tips);
	// for each gauss point in the element
	int nint = se->GaussPoints();
	for (int i = 0; i < nint; i++)
	{
		// get the material point
		FEMaterialPoint* mp = se->GetMaterialPoint(i);
		// get the angio material point
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		// get the elastic material point
		FEElasticMaterialPoint* emp = mp->ExtractData<FEElasticMaterialPoint>();
		assert(mp && angio_mp);
		angio_mp->m_as.zero();
		// get the global position of the material point
		vec3d y = mp->m_rt;
		int ntips = final_active_tips.size();

		for (int j = 0; j < ntips; j++)
		{
			Tip* tip = final_active_tips[j];
			// get global position of tip
			//consider moving this out and calling it less
			vec3d x = tip->GetPosition(mesh);
			// get vector from tip to gauss point
			vec3d r = y - x;
			double l = r.unit();
			// get angle between tip direction and r
			//consider moving this out and calling it less
			double theta = acos(tip->GetDirection(mesh) * r);

			//sprout s mag replaced with giving correct coeficients for scale
			//! this growth velocity may need to be a scaling term that is a function of velocity
			double p = tip->growth_velocity * sprout_mag * pow(cos(theta / 2.0), fan_exponential) * exp(-l / sprout_range);
			angio_mp->m_as += dyad(r) * p;
		}
	}
}

bool LoadCurveAngioStressPolicy::Init()
{
	return true;
}

void LoadCurveAngioStressPolicy::UpdateScale()
{
	//! Implementation blank
}

void LoadCurveAngioStressPolicy::AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh)
{
	std::vector<Tip*> final_active_tips;
	auto se = angio_element->_elem;
	double total_sprout_falloff = sprout_range * sprout_radius_multiplier;
	FEAngio::GetActiveFinalTipsInRadius(angio_element, total_sprout_falloff, pangio, final_active_tips);

	int nint = se->GaussPoints();
	for (int i = 0; i < nint; i++)
	{
		FEMaterialPoint* mp = se->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		FEElasticMaterialPoint* emp = mp->ExtractData<FEElasticMaterialPoint>();
		assert(mp && angio_mp);
		angio_mp->m_as.zero();
		vec3d y = mp->m_rt;
		int ntips = final_active_tips.size();
		for (int j = 0; j < ntips; j++)
		{
			Tip* tip = final_active_tips[j];
			//consider moving this out and calling it less
			vec3d x = tip->GetPosition(mesh);
			vec3d r = y - x;
			double l = r.unit();
			//consider moving this out and calling it less
			double theta = acos(tip->GetDirection(mesh) * r);

			//sprout s mag replaced with giving correct coeficients for scale
			double p = sprout_mag * pow(cos(theta / 2.0), fan_exponential) * exp(-l / sprout_range);
			angio_mp->m_as += dyad(r) * p;
		}
	}
}

bool LoadCurveDenAngioStressPolicy::Init()
{
	return true;
}

void LoadCurveDenAngioStressPolicy::UpdateScale()
{
	//! Implementation blank
}

void LoadCurveDenAngioStressPolicy::AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh)
{
	std::vector<Tip*> final_active_tips;
	auto se = angio_element->_elem;
	double total_sprout_falloff = sprout_range * sprout_radius_multiplier;
	FEAngio::GetActiveFinalTipsInRadius(angio_element, total_sprout_falloff, pangio, final_active_tips);

	int nint = se->GaussPoints();
	for (int i = 0; i < nint; i++)
	{
		FEMaterialPoint* mp = se->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		FEElasticMaterialPoint* emp = mp->ExtractData<FEElasticMaterialPoint>();
		assert(mp && angio_mp);
		angio_mp->m_as.zero();
		vec3d y = mp->m_rt;
		int ntips = final_active_tips.size();
		for (int j = 0; j < ntips; j++)
		{
			Tip* tip = final_active_tips[j];
			//consider moving this out and calling it less
			vec3d x = tip->GetPosition(mesh);
			vec3d r = y - x;
			double l = r.unit();
			//consider moving this out and calling it less
			double theta = acos(tip->GetDirection(mesh) * r);

			//sprout s mag replaced with giving correct coeficients for scale
			double density_scale = this->GetDensScale(angio_element, tip, mesh, 1);
			double p = density_scale * sprout_mag * pow(cos(theta / 2.0), fan_exponential) * exp(-l / sprout_range);
			angio_mp->m_as += dyad(r) * p;
		}
	}
}

bool LoadCurveRefDenAngioStressPolicy::Init()
{
	return true;
}

void LoadCurveRefDenAngioStressPolicy::UpdateScale()
{
	//! Implementation blank
}

void LoadCurveRefDenAngioStressPolicy::AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh)
{
	std::vector<Tip*> final_active_tips;
	auto se = angio_element->_elem;
	double total_sprout_falloff = sprout_range * sprout_radius_multiplier;
	FEAngio::GetActiveFinalTipsInRadius(angio_element, total_sprout_falloff, pangio, final_active_tips);

	int nint = se->GaussPoints();
	for (int i = 0; i < nint; i++)
	{
		FEMaterialPoint* mp = se->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		FEElasticMaterialPoint* emp = mp->ExtractData<FEElasticMaterialPoint>();
		assert(mp && angio_mp);
		angio_mp->m_as.zero();
		vec3d y = mp->m_rt;
		int ntips = final_active_tips.size();

		for (int j = 0; j < ntips; j++)
		{
			Tip* tip = final_active_tips[j];
			//consider moving this out and calling it less
			vec3d x = tip->GetPosition(mesh);
			vec3d r = y - x;
			double l = r.unit();
			//consider moving this out and calling it less
			double theta = acos(tip->GetDirection(mesh) * r);

			//sprout s mag replaced with giving correct coeficients for scale
			double density_scale = this->GetDensScale(angio_element, tip, mesh, 0);
			double p = density_scale * sprout_mag * pow(cos(theta / 2.0), fan_exponential) * exp(-l / sprout_range);
			angio_mp->m_as += dyad(r) * p;
		}
	}
}

bool GrownSegmentsAngioStressPolicy::Init()
{
	return true;
}

void GrownSegmentsAngioStressPolicy::UpdateScale()
{
	//! Implementation blank
}

void GrownSegmentsAngioStressPolicy::
	AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh)
{
	std::vector<Tip*> grown_tips;
	auto se = angio_element->_elem;
	double total_sprout_falloff = sprout_range * sprout_radius_multiplier;
	FEAngio::GetGrownTipsInRadius(angio_element, total_sprout_falloff, pangio, grown_tips);
	int nint = se->GaussPoints();
	for (int i = 0; i < nint; i++)
	{
		FEMaterialPoint* mp = se->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		FEElasticMaterialPoint* emp = mp->ExtractData<FEElasticMaterialPoint>();
		assert(mp && angio_mp);
		angio_mp->m_as.zero();
		vec3d y = mp->m_rt;
		for (int j = 0; j < grown_tips.size(); j++)
		{
			Tip* tip = grown_tips[j];
			//consider moving this out and calling it less
			vec3d x = tip->GetPosition(mesh);
			vec3d r = y - x;
			double l = r.unit();
			double theta = acos(tip->GetDirection(mesh) * r);
			//consider moving this out and calling it less
			double p = sprout_mag * pow(cos(theta / 2.0), fan_exponential) * exp(-l / sprout_range);
			angio_mp->m_as += dyad(r) * p;
		}
	}
}


bool GrownSegmentsVelAngioStressPolicy::Init()
{
	return true;
}

void GrownSegmentsVelAngioStressPolicy::UpdateScale()
{
	//! Implementation blank
}

void GrownSegmentsVelAngioStressPolicy::AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh)
{
	std::vector<Tip*> grown_tips;
	auto se = angio_element->_elem;
	double total_sprout_falloff = sprout_range * sprout_radius_multiplier;
	FEAngio::GetGrownTipsInRadius(angio_element, total_sprout_falloff, pangio, grown_tips);
	
	int nint = se->GaussPoints();
	for (int i = 0; i < nint; i++)
	{
		FEMaterialPoint* mp = se->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		FEElasticMaterialPoint* emp = mp->ExtractData<FEElasticMaterialPoint>();
		assert(mp && angio_mp);
		angio_mp->m_as.zero();
		vec3d y = mp->m_rt;
		for (int j = 0; j < grown_tips.size(); j++)
		{
			Tip* tip = grown_tips[j];
			//consider moving this out and calling it less
			vec3d x = tip->GetPosition(mesh);
			vec3d r = y - x;
			double l = r.unit();
			//consider moving this out and calling it less
			double theta = acos(tip->GetDirection(mesh) * r);
			double p = tip->growth_velocity * sprout_mag * pow(cos(theta / 2.0), fan_exponential) * exp(-l / sprout_range);
			angio_mp->m_as += dyad(r) * p;
		}
	}
}