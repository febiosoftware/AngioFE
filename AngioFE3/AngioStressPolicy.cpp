#include "AngioStressPolicy.h"
#include "AngioElement.h"
#include "FEAngio.h"
#include "Tip.h"
#include <FECore/FEDataLoadCurve.h>

void AngioStressPolicy::GetActiveFinalTipsInRadius(AngioElement* angio_element, double radius, FEAngio* pangio, std::vector<Tip *> & tips)
{
	std::set<AngioElement *> visited;
	std::set<AngioElement *> next;
	std::vector<vec3d> element_bounds;
	pangio->ExtremaInElement(angio_element->_elem, element_bounds);

	for(int i=0; i < angio_element->adjacency_list.size();i++)
	{
		next.insert(angio_element->adjacency_list[i]);
	}
	visited.insert(angio_element);
	while(next.size())
	{
		AngioElement * cur = *next.begin();
		next.erase(next.begin());
		visited.insert(cur);
		std::vector<vec3d> cur_element_bounds;
		pangio->ExtremaInElement(cur->_elem, cur_element_bounds);
		double cdist = FEAngio::MinDistance(element_bounds, cur_element_bounds);
		if(cdist <= radius)
		{
			//add the tips and the add all unvisited adjacent elements to next
			for(int i =0; i < cur->final_active_tips.size();i++)
			{
				tips.push_back(cur->final_active_tips[i]);
			}

			for(int i=0;i < cur->adjacency_list.size();i++)
			{
				if(!visited.count(cur->adjacency_list[i]))
				{
					next.insert(cur->adjacency_list[i]);
				}
			}
		}
	}
}
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
	GetActiveFinalTipsInRadius(angio_element, sprout_range * sprout_radius_multiplier, pangio, final_active_tips);

	for(int i=0;i < angio_element->_elem->GaussPoints();i++)
	{
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint * angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		double den_scale = angio_element->_angio_mat->FindDensityScale(angio_mp);
		assert(mp && angio_mp);
		angio_mp->m_as.zero();
		vec3d y =  emp->m_rt;
		for(int j=0; j < final_active_tips.size();j++)
		{
			Tip * tip = final_active_tips[j];
			vec3d x = tip->GetPosition(mesh);//consider moving this out and calling it less
			vec3d r = y - x;
			double l = r.unit();
			double theta = acos(tip->GetDirection(mesh) * r);//same for GetDirection

			//sprout s mag replaced with giving correct coeficients for scale
			double p = den_scale * scale *sprout_mag * pow(cos(theta / 2), sprout_width)* exp(-l / sprout_range);
			angio_mp->m_as += dyad(r)*p;
		}
	}
	
}

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
	GetActiveFinalTipsInRadius(angio_element, sprout_range * sprout_radius_multiplier, pangio, final_active_tips);

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
			double p = tip->growth_velocity *sprout_mag * pow(cos(theta / 2), sprout_width)* exp(-l / sprout_range);
			angio_mp->m_as += dyad(r)*p;
		}
	}

}


BEGIN_PARAMETER_LIST(LoadCurveVelAngioStressPolicy, AngioStressPolicy)
ADD_PARAMETER(sprout_mag, FE_PARAM_DOUBLE, "sprout_mag");
ADD_PARAMETER(sprout_width, FE_PARAM_DOUBLE, "sprout_width");
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
	GetActiveFinalTipsInRadius(angio_element, sprout_range * sprout_radius_multiplier, pangio, final_active_tips);

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
			double p = sprout_mag * pow(cos(theta / 2), sprout_width)* exp(-l / sprout_range);
			angio_mp->m_as += dyad(r)*p;
		}
	}

}

BEGIN_PARAMETER_LIST(LoadCurveAngioStressPolicy, AngioStressPolicy)
ADD_PARAMETER(sprout_mag, FE_PARAM_DOUBLE, "sprout_mag");
ADD_PARAMETER(sprout_width, FE_PARAM_DOUBLE, "sprout_width");
ADD_PARAMETER(sprout_range, FE_PARAM_DOUBLE, "sprout_range");
ADD_PARAMETER(sprout_radius_multiplier, FE_PARAM_DOUBLE, "sprout_radius_multiplier");
END_PARAMETER_LIST();