#include "AngioStressPolicy.h"
#include "AngioElement.h"
#include "FEAngio.h"
#include "Tip.h"
#include <FECore/FEDataLoadCurve.h>

void AngioStressPolicy::UpdateToLoadCurve(const char* param_name, double& value)
{
	//if load curves are used they must use step interpolation
	FEParam * m = FindParameter(ParamString(param_name));
	assert(m);
	int mlci = m->GetLoadCurve();

	FEModel * model = GetFEModel();
	if (mlci >= 0)
	{
		FEDataLoadCurve * mlc = dynamic_cast<FEDataLoadCurve*>(model->GetLoadCurve(mlci));
		assert(mlc);
		value = mlc->Value(model->GetTime().currentTime);
	}
}

double AngioStressPolicy::GetDensScale(AngioElement* angio_element, Tip* tip, FEMesh* mesh, int is_ref) {
	std::vector<double> density_at_integration_points;

	//this is really slow and it isn't apparent why.
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		FEElasticMaterialPoint* elastic_mp = gauss_point->ExtractData<FEElasticMaterialPoint>();
		// determines whether to use the referential or current density
		double m_J = pow(elastic_mp->m_J,is_ref);
		density_at_integration_points.push_back(angio_mp->ref_ecm_density * (1.0 / m_J));
	}
	//Get interpolation method
	FEModel* m_pfem = this->GetFEModel();
	PerElementVI interp = PerElementVI(m_pfem);
	double density_at_point = interp.Interpolate(angio_element->_elem, density_at_integration_points, tip->GetLocalPosition(), mesh);
	double density_scale = m_density_scale_factor.x + m_density_scale_factor.y * exp(-m_density_scale_factor.z * density_at_point);
	if (density_scale < 0) { return 0; }
	else { return density_scale; }
}

bool SigmoidAngioStressPolicy::Init()
{
	return true;
}

void SigmoidAngioStressPolicy::UpdateScale()
{
	double time = GetFEModel()->GetTime().currentTime;
	scale = y0 + a / (1 + exp(-(time - x0) / b));;
}

void SigmoidAngioStressPolicy::AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh)
{
	std::vector<Tip*> final_active_tips;
	// Get active tips within a radius
	FEAngio::GetActiveFinalTipsInRadius(angio_element, sprout_range * sprout_radius_multiplier, pangio, final_active_tips);
	// for each integration point
	for(int i=0;i < angio_element->_elem->GaussPoints();i++)
	{
		// get the material point of the gauss point
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		// get eh angio material point at the gauss point
		FEAngioMaterialPoint * angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		// get the elastic material point at the gauss point
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		// get the density scale for the integration point
		double den_scale = angio_element->_angio_mat->FindDensityScale(angio_mp);
		// assert that the material point and gauss point are in the same point
		assert(mp && angio_mp);
		// zero the angio point
		angio_mp->m_as.zero();
		// get global position of elastic material point?
		vec3d y =  emp->m_rt;
		// for each active tip
		for(int j=0; j < final_active_tips.size();j++)
		{
			Tip * tip = final_active_tips[j];
			// get the tip position in the mesh
			vec3d x = tip->GetPosition(mesh);//consider moving this out and calling it less
			// determine vector from tip to integration point
			vec3d r = y - x;
			double l = r.unit();
			// get the angle between the tip and the sprout direction
			double theta = acos(tip->GetDirection(mesh) * r);//same for GetDirection

			//sprout s mag replaced with giving correct coeficients for scale
			double p = den_scale * scale *sprout_mag * pow(cos(theta / 2), fan_exponential)* exp(-l / sprout_range);
			// make a dyad times the pressure
			angio_mp->m_as += dyad(r)*p;
		}
	}
	
}

BEGIN_PARAMETER_LIST(SigmoidAngioStressPolicy, AngioStressPolicy)
ADD_PARAMETER(sprout_mag, FE_PARAM_DOUBLE, "sprout_mag");
ADD_PARAMETER(fan_exponential, FE_PARAM_DOUBLE, "fan_exponential");
ADD_PARAMETER(sprout_range, FE_PARAM_DOUBLE, "sprout_range");
ADD_PARAMETER(sprout_radius_multiplier, FE_PARAM_DOUBLE, "sprout_radius_multiplier");
ADD_PARAMETER(a, FE_PARAM_DOUBLE, "a");
ADD_PARAMETER(b, FE_PARAM_DOUBLE, "b");
ADD_PARAMETER(x0, FE_PARAM_DOUBLE, "x0");
ADD_PARAMETER(y0, FE_PARAM_DOUBLE, "y0");

END_PARAMETER_LIST();


bool LoadCurveVelAngioStressPolicy::Init()
{
	return true;
}

void LoadCurveVelAngioStressPolicy::UpdateScale()
{
}

void LoadCurveVelAngioStressPolicy::AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh)
{
	std::vector<Tip*> final_active_tips;
	// get all tips within a radius
	FEAngio::GetActiveFinalTipsInRadius(angio_element, sprout_range * sprout_radius_multiplier, pangio, final_active_tips);
	// for each gauss point in the element
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		// get the material point
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		// get the angio material point
		FEAngioMaterialPoint * angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		// get the elastic material point
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		assert(mp && angio_mp);
		angio_mp->m_as.zero();
		// get the global position of the material point
		vec3d y = emp->m_rt;
		for (int j = 0; j < final_active_tips.size(); j++)
		{
			Tip * tip = final_active_tips[j];
			// get global position of tip
			vec3d x = tip->GetPosition(mesh);//consider moving this out and calling it less
			// get vector from tip to gauss point
			vec3d r = y - x;
			double l = r.unit();
			// get angle between tip direction and r
			double theta = acos(tip->GetDirection(mesh) * r);//same for GetDirection

			//sprout s mag replaced with giving correct coeficients for scale
			// this growth velocity may need to be a scaling term that is a function of velocity
			double p = tip->growth_velocity *sprout_mag * pow(cos(theta / 2), fan_exponential)* exp(-l / sprout_range);
			angio_mp->m_as += dyad(r)*p;
		}
	}

}


BEGIN_PARAMETER_LIST(LoadCurveVelAngioStressPolicy, AngioStressPolicy)
ADD_PARAMETER(sprout_mag, FE_PARAM_DOUBLE, "sprout_mag");
ADD_PARAMETER(fan_exponential, FE_PARAM_DOUBLE, "fan_exponential");
ADD_PARAMETER(sprout_range, FE_PARAM_DOUBLE, "sprout_range");
ADD_PARAMETER(sprout_radius_multiplier, FE_PARAM_DOUBLE, "sprout_radius_multiplier");
END_PARAMETER_LIST();

bool LoadCurveAngioStressPolicy::Init()
{
	return true;
}

void LoadCurveAngioStressPolicy::UpdateScale()
{

}

void LoadCurveAngioStressPolicy::AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh)
{
	std::vector<Tip*> final_active_tips;
	FEAngio::GetActiveFinalTipsInRadius(angio_element, sprout_range * sprout_radius_multiplier, pangio, final_active_tips);

	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint * angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		assert(mp && angio_mp);
		angio_mp->m_as.zero();
		vec3d y = emp->m_rt;
		for (int j = 0; j < final_active_tips.size(); j++)
		{
			Tip * tip = final_active_tips[j];
			vec3d x = tip->GetPosition(mesh);//consider moving this out and calling it less
			vec3d r = y - x;
			double l = r.unit();
			double theta = acos(tip->GetDirection(mesh) * r);//same for GetDirection

			//sprout s mag replaced with giving correct coeficients for scale
			double p = sprout_mag * pow(cos(theta / 2), fan_exponential)* exp(-l / sprout_range);
			angio_mp->m_as += dyad(r)*p;
		}
	}

}

BEGIN_PARAMETER_LIST(LoadCurveAngioStressPolicy, AngioStressPolicy)
ADD_PARAMETER(sprout_mag, FE_PARAM_DOUBLE, "sprout_mag");
ADD_PARAMETER(fan_exponential, FE_PARAM_DOUBLE, "fan_exponential");
ADD_PARAMETER(sprout_range, FE_PARAM_DOUBLE, "sprout_range");
ADD_PARAMETER(sprout_radius_multiplier, FE_PARAM_DOUBLE, "sprout_radius_multiplier");
END_PARAMETER_LIST();

bool LoadCurveDenAngioStressPolicy::Init()
{
	return true;
}

void LoadCurveDenAngioStressPolicy::UpdateScale()
{

}

void LoadCurveDenAngioStressPolicy::AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh)
{
	std::vector<Tip*> final_active_tips;
	FEAngio::GetActiveFinalTipsInRadius(angio_element, sprout_range * sprout_radius_multiplier, pangio, final_active_tips);

	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint * angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		assert(mp && angio_mp);
		angio_mp->m_as.zero();
		vec3d y = emp->m_rt;
		for (int j = 0; j < final_active_tips.size(); j++)
		{
			Tip * tip = final_active_tips[j]; 
			vec3d x = tip->GetPosition(mesh);//consider moving this out and calling it less
			vec3d r = y - x;
			double l = r.unit();
			double theta = acos(tip->GetDirection(mesh) * r);//same for GetDirection

			//sprout s mag replaced with giving correct coeficients for scale
			double density_scale = this->GetDensScale(angio_element, tip, mesh, 1);
			double p = density_scale * sprout_mag * pow(cos(theta / 2), fan_exponential)* exp(-l / sprout_range);
			angio_mp->m_as += dyad(r)*p;
		}
	}

}

BEGIN_PARAMETER_LIST(LoadCurveDenAngioStressPolicy, AngioStressPolicy)
ADD_PARAMETER(sprout_mag, FE_PARAM_DOUBLE, "sprout_mag");
ADD_PARAMETER(fan_exponential, FE_PARAM_DOUBLE, "fan_exponential");
ADD_PARAMETER(sprout_range, FE_PARAM_DOUBLE, "sprout_range");
ADD_PARAMETER(sprout_radius_multiplier, FE_PARAM_DOUBLE, "sprout_radius_multiplier");
END_PARAMETER_LIST();

bool LoadCurveRefDenAngioStressPolicy::Init()
{
	return true;
}

void LoadCurveRefDenAngioStressPolicy::UpdateScale()
{

}

void LoadCurveRefDenAngioStressPolicy::AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh)
{
	std::vector<Tip*> final_active_tips;
	FEAngio::GetActiveFinalTipsInRadius(angio_element, sprout_range * sprout_radius_multiplier, pangio, final_active_tips);

	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint * angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		assert(mp && angio_mp);
		angio_mp->m_as.zero();
		vec3d y = emp->m_rt;
		for (int j = 0; j < final_active_tips.size(); j++)
		{
			Tip * tip = final_active_tips[j];
			vec3d x = tip->GetPosition(mesh);//consider moving this out and calling it less
			vec3d r = y - x;
			double l = r.unit();
			double theta = acos(tip->GetDirection(mesh) * r);//same for GetDirection

			//sprout s mag replaced with giving correct coeficients for scale
			double density_scale = this->GetDensScale(angio_element, tip, mesh, 0);
			double p = density_scale * sprout_mag * pow(cos(theta / 2), fan_exponential)* exp(-l / sprout_range);
			angio_mp->m_as += dyad(r)*p;
		}
	}

}

BEGIN_PARAMETER_LIST(LoadCurveRefDenAngioStressPolicy, AngioStressPolicy)
ADD_PARAMETER(sprout_mag, FE_PARAM_DOUBLE, "sprout_mag");
ADD_PARAMETER(fan_exponential, FE_PARAM_DOUBLE, "fan_exponential");
ADD_PARAMETER(sprout_range, FE_PARAM_DOUBLE, "sprout_range");
ADD_PARAMETER(sprout_radius_multiplier, FE_PARAM_DOUBLE, "sprout_radius_multiplier");
END_PARAMETER_LIST();

bool GrownSegmentsAngioStressPolicy::Init()
{
	return true;
}

void GrownSegmentsAngioStressPolicy::UpdateScale()
{
}

void GrownSegmentsAngioStressPolicy::AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh)
{
	std::vector<Tip*> grown_tips;
	FEAngio::GetGrownTipsInRadius(angio_element, sprout_range * sprout_radius_multiplier, pangio, grown_tips);

	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint * angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		assert(mp && angio_mp);
		angio_mp->m_as.zero();
		vec3d y = emp->m_rt;
		for (int j = 0; j < grown_tips.size(); j++)
		{
			Tip * tip = grown_tips[j];
			vec3d x = tip->GetPosition(mesh);//consider moving this out and calling it less
			vec3d r = y - x;
			double l = r.unit();
			double theta = acos(tip->GetDirection(mesh) * r);//same for GetDirection

															 //sprout s mag replaced with giving correct coeficients for scale
			double p = sprout_mag * pow(cos(theta / 2), fan_exponential)* exp(-l / sprout_range);
			angio_mp->m_as += dyad(r)*p;
		}
	}

}


BEGIN_PARAMETER_LIST(GrownSegmentsAngioStressPolicy, AngioStressPolicy)
ADD_PARAMETER(sprout_mag, FE_PARAM_DOUBLE, "sprout_mag");
ADD_PARAMETER(fan_exponential, FE_PARAM_DOUBLE, "fan_exponential");
ADD_PARAMETER(sprout_range, FE_PARAM_DOUBLE, "sprout_range");
ADD_PARAMETER(sprout_radius_multiplier, FE_PARAM_DOUBLE, "sprout_radius_multiplier");
END_PARAMETER_LIST();


bool GrownSegmentsVelAngioStressPolicy::Init()
{
	return true;
}

void GrownSegmentsVelAngioStressPolicy::UpdateScale()
{
}

void GrownSegmentsVelAngioStressPolicy::AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh)
{
	std::vector<Tip*> grown_tips;
	FEAngio::GetGrownTipsInRadius(angio_element, sprout_range * sprout_radius_multiplier, pangio, grown_tips);

	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint * angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		assert(mp && angio_mp);
		angio_mp->m_as.zero();
		vec3d y = emp->m_rt;
		for (int j = 0; j < grown_tips.size(); j++)
		{
			Tip * tip = grown_tips[j];
			vec3d x = tip->GetPosition(mesh);//consider moving this out and calling it less
			vec3d r = y - x;
			double l = r.unit();
			double theta = acos(tip->GetDirection(mesh) * r);//same for GetDirection

															 //sprout s mag replaced with giving correct coeficients for scale
			double p = tip->growth_velocity *sprout_mag * pow(cos(theta / 2), fan_exponential)* exp(-l / sprout_range);
			angio_mp->m_as += dyad(r)*p;
		}
	}

}


BEGIN_PARAMETER_LIST(GrownSegmentsVelAngioStressPolicy, AngioStressPolicy)
ADD_PARAMETER(sprout_mag, FE_PARAM_DOUBLE, "sprout_mag");
ADD_PARAMETER(fan_exponential, FE_PARAM_DOUBLE, "fan_exponential");
ADD_PARAMETER(sprout_range, FE_PARAM_DOUBLE, "sprout_range");
ADD_PARAMETER(sprout_radius_multiplier, FE_PARAM_DOUBLE, "sprout_radius_multiplier");
END_PARAMETER_LIST();