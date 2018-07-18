#include "PositionDependentDirection.h"
#include "angio3d.h"
#include "Tip.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FESolidDomain.h>
#include "FEAngioMaterialPoint.h"
#include "FEAngio.h"


BEGIN_PARAMETER_LIST(PositionDependentDirection, FEMaterial)
ADD_PARAMETER(contribution, FE_PARAM_DOUBLE, "contribution");
END_PARAMETER_LIST();

vec3d PositionDependentDirectionManager::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int buffer, double& alpha, FEMesh* mesh, FEAngio* pangio)
{
	for (int i = 0; i < pdd_modifiers.size(); i++)
	{
		prev = pdd_modifiers[i]->ApplyModifiers(prev, angio_element, local_pos, initial_fragment_id , buffer , alpha, mesh, pangio);
	}
	return prev;
}

BEGIN_PARAMETER_LIST(ECMDensityGradientPDD, PositionDependentDirection)
ADD_PARAMETER(threshold, FE_PARAM_DOUBLE, "threshold");
ADD_PARAMETER(alpha_override, FE_PARAM_BOOL, "alpha_override");
END_PARAMETER_LIST();

vec3d ECMDensityGradientPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, FEMesh* mesh, FEAngio* pangio)
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
	if (gradnorm > threshold)
	{
		vec3d currentDirectionGradientPlane = grad ^ prev;
		currentDirectionGradientPlane.unit();
		vec3d perpendicularToGradient = currentDirectionGradientPlane ^ grad;
		perpendicularToGradient.unit();
		if (alpha_override)
		{
			alpha = contribution;
		}
		return mix(prev, perpendicularToGradient, contribution);
	}
	return prev;
}

BEGIN_PARAMETER_LIST(RepulsePDD, PositionDependentDirection)
ADD_PARAMETER(threshold, FE_PARAM_DOUBLE, "threshold");
ADD_PARAMETER(alpha_override, FE_PARAM_BOOL, "alpha_override");
END_PARAMETER_LIST();

vec3d RepulsePDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, FEMesh* mesh, FEAngio* pangio)
{
	std::vector<double> repulse_at_integration_points;

	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		FEElasticMaterialPoint* elastic_mp = gauss_point->ExtractData<FEElasticMaterialPoint>();

		repulse_at_integration_points.push_back(angio_mp->repulse_value * (1.0 / elastic_mp->m_J));
	}
	if(grad_threshold)
	{
		vec3d grad = pangio->gradient(angio_element->_elem, repulse_at_integration_points, local_pos);
		double gradnorm = grad.norm();
		if (gradnorm > threshold)
		{
			grad = -grad;
			if (alpha_override)
			{
				alpha = contribution;
			}
			return mix(prev,grad, contribution);
		}
	}
	else
	{
		double nodal_repulse[FEElement::MAX_NODES];
		double H[FEElement::MAX_NODES];
		angio_element->_elem->project_to_nodes(&repulse_at_integration_points[0],nodal_repulse);
		double val = 0;
		angio_element->_elem->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
		for(int i=0; i < angio_element->_elem->Nodes();i++)
		{
			val += H[i] * nodal_repulse[i];
		}
		vec3d grad = pangio->gradient(angio_element->_elem, repulse_at_integration_points, local_pos);
		if(val > threshold)
		{
			grad = -grad;
			if (alpha_override)
			{
				alpha = contribution;
			}
			return mix(prev, grad, contribution);
		}
	}
	return prev;
	
}

BEGIN_PARAMETER_LIST(ConcentrationGradientPDD, PositionDependentDirection)
ADD_PARAMETER(threshold, FE_PARAM_DOUBLE, "threshold");
ADD_PARAMETER(alpha_override, FE_PARAM_BOOL, "alpha_override");
ADD_PARAMETER(sol_id, FE_PARAM_INT, "sol_id");
END_PARAMETER_LIST();

vec3d ConcentrationGradientPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, FEMesh* mesh, FEAngio* pangio)
{
	std::vector<double> concentration_at_integration_points;

	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);

		concentration_at_integration_points.push_back(pangio->GetConcentration(angio_element->_mat, gauss_point, sol_id));
	}
	vec3d grad = pangio->gradient(angio_element->_elem, concentration_at_integration_points, local_pos);
	double gradnorm = grad.norm();
	if (gradnorm > threshold)
	{
		if (alpha_override)
		{
			alpha = contribution;
		}
		return mix(prev, grad, contribution);
	}
	return prev;
}

//used in the anastamosis modifier
double distance2(FESolidElement * se, vec3d local_pos, Tip * tip, FEMesh* mesh)
{
	vec3d pos;
	double H[FEElement::MAX_NODES];
	se->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
	for(int i=0; i < se->Nodes();i++)
	{
		pos += (mesh->Node(se->m_node[i]).m_rt * H[i]);
	}
	return (tip->GetPosition(mesh) - pos).norm2();
}

vec3d AnastomosisPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, FEMesh* mesh, FEAngio* pangio)
{
	//find all elements within the radius 
	std::vector<Tip*> tips;
	FEAngio::GetActiveTipsInRadius(angio_element, anastomosis_radius, current_buffer, pangio, tips, initial_fragment_id);
	//if there are any valid tips pick the nearest and grow towards it
	if(tips.size())
	{
		int best_index = 0;
		double best_distance = distance2(angio_element->_elem,local_pos, tips[0], mesh);
		for(int i=1; i < tips.size();i++)
		{
			double pd = distance2(angio_element->_elem, local_pos, tips[i], mesh);
			if(best_distance > pd)
			{
				best_index = i;
				best_distance = pd;
			}
		}
		//get the direction
		vec3d pos;
		double H[FEElement::MAX_NODES];
		angio_element->_elem->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
		for (int i = 0; i < angio_element->_elem->Nodes(); i++)
		{
			pos += (mesh->Node(angio_element->_elem->m_node[i]).m_rt * H[i]);
		}
		vec3d dir = (tips[best_index]->GetPosition(mesh) - pos);
		dir.unit();
		if (alpha_override)
		{
			alpha = contribution;
		}
		return mix(prev, dir, contribution);
	}


	return prev;
}

BEGIN_PARAMETER_LIST(AnastomosisPDD, PositionDependentDirection)
ADD_PARAMETER(anastomosis_radius, FE_PARAM_DOUBLE, "anastomosis_radius");
ADD_PARAMETER(contribution, FE_PARAM_DOUBLE, "contribution");
END_PARAMETER_LIST();

void FiberPDD::Update(FEMesh * mesh, FEAngio* angio)
{

}

vec3d FiberPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, FEMesh* mesh, FEAngio* pangio)
{
	std::vector<quatd> gauss_data;
	std::vector<quatd> nodal_data;
	for (int i = 0; i< angio_element->_elem->GaussPoints(); i++)
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