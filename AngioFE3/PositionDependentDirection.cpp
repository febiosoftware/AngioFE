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

vec3d PositionDependentDirectionManager::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int buffer, bool& continue_growth, vec3d& tip_dir, double& alpha, FEMesh* mesh, FEAngio* pangio)
{
	for (int i = 0; i < pdd_modifiers.size(); i++)
	{
		prev = pdd_modifiers[i]->ApplyModifiers(prev, angio_element, local_pos, initial_fragment_id , buffer , alpha, continue_growth , tip_dir , mesh, pangio);
	}
	return prev;
}

void PositionDependentDirectionManager::Update(FEMesh * mesh, FEAngio* angio)
{
	for (int i = 0; i < pdd_modifiers.size(); i++)
	{
		pdd_modifiers[i]->Update(mesh, angio);
	}
}

BEGIN_PARAMETER_LIST(ECMDensityGradientPDD, PositionDependentDirection)
ADD_PARAMETER(threshold, FE_PARAM_DOUBLE, "threshold");
ADD_PARAMETER(alpha_override, FE_PARAM_BOOL, "alpha_override");
END_PARAMETER_LIST();

vec3d ECMDensityGradientPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continute_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio)
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
		return mix3d(prev, perpendicularToGradient, contribution);
		//return mix(prev, perpendicularToGradient, contribution);
	}
	return prev;
}

BEGIN_PARAMETER_LIST(RepulsePDD, PositionDependentDirection)
ADD_PARAMETER(threshold, FE_PARAM_DOUBLE, "threshold");
ADD_PARAMETER(alpha_override, FE_PARAM_BOOL, "alpha_override");
ADD_PARAMETER(grad_threshold, FE_PARAM_BOOL, "grad_threshold");
END_PARAMETER_LIST();

vec3d RepulsePDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continute_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio)
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
			return mix3d(prev, grad, contribution);
			// return mix(prev,grad, contribution);
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
			return mix3d(prev, grad, contribution);
			// return mix(prev, grad, contribution);
		}
	}
	return prev;
	
}

BEGIN_PARAMETER_LIST(ConcentrationGradientPDD, PositionDependentDirection)
ADD_PARAMETER(threshold, FE_PARAM_DOUBLE, "threshold");
ADD_PARAMETER(alpha_override, FE_PARAM_BOOL, "alpha_override");
ADD_PARAMETER(sol_id, FE_PARAM_INT, "sol_id");
END_PARAMETER_LIST();

vec3d ConcentrationGradientPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continute_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio)
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
		return mix3d(prev, grad, contribution);
		// return mix(prev, grad, contribution);
	}
	return prev;
}

//used in the anastamosis modifier
double AnastamosisPDD::distance2(FESolidElement * se, vec3d local_pos, Tip * tip, FEMesh* mesh)
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

vec3d AnastamosisPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continute_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio)
{
	vec3d tip_pos;
	double H[FEElement::MAX_NODES];
	angio_element->_elem->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
	for (int i = 0; i < angio_element->_elem->Nodes(); i++)
	{
		tip_pos += (mesh->Node(angio_element->_elem->m_node[i]).m_rt * H[i]);
	}
	Tip* best_tip = FuseWith(angio_element, pangio, mesh, tip_pos, tip_dir, initial_fragment_id, anastamosis_radius);
	if(best_tip)
	{
		vec3d dir = (best_tip->GetPosition(mesh) - tip_pos);
		double len = dir.unit();
		if (alpha_override)
		{
			alpha = contribution;
		}
		if (len < fuse_radius && fuse_angle  > dir * tip_dir)
		{
			angio_element->anastamoses++;
			continute_growth = false;
		}
		return mix3d(prev, dir, contribution);
		// return mix(prev, dir, contribution);
	}
	return prev;
}

Tip* AnastamosisPDD::FuseWith(class AngioElement* angio_element, FEAngio* pangio, class FEMesh* mesh, class vec3d tip_pos, class vec3d tip_dir, int exclude, double radius)
{
	double best_possible_distance;
	Tip * best = BestInElement(angio_element, pangio, mesh, tip_pos,tip_dir, exclude, best_possible_distance);
	double best_distance = std::numeric_limits<double>::max();
	if (best)
	{
		best_distance = best_possible_distance;
	}

	std::unordered_set<AngioElement *> visited;
	std::set<AngioElement *> next;
	std::vector<vec3d> element_bounds;
	visited.reserve(1000);
	pangio->ExtremaInElement(angio_element->_elem, element_bounds);



	for (int i = 0; i < angio_element->adjacency_list.size(); i++)
	{
		next.insert(angio_element->adjacency_list[i]);
	}
	visited.insert(angio_element);
	while (next.size())
	{
		AngioElement * cur = *next.begin();
		next.erase(next.begin());
		visited.insert(cur);
		std::vector<vec3d> cur_element_bounds;
		pangio->ExtremaInElement(cur->_elem, cur_element_bounds);
		double cdist = FEAngio::MinDistance(element_bounds, cur_element_bounds);
		if (cdist <= radius && cdist <= best_distance)
		{
			Tip * pos_best = BestInElement(cur, pangio, mesh, tip_pos, tip_dir, exclude, best_possible_distance);
			if (pos_best && best_possible_distance < best_distance)
			{
				best = pos_best;
				best_distance = best_possible_distance;
			}


			for (int i = 0; i < cur->adjacency_list.size(); i++)
			{
				if (!visited.count(cur->adjacency_list[i]))
				{
					next.insert(cur->adjacency_list[i]);
				}
			}
		}
		//otherwise let uneeded elements die
	}
	return best;
}

bool AnastamosisPDD::ValidTip(Tip* tip, vec3d tip_dir, FEMesh* mesh)
{
	return fuse_angle > tip->GetDirection(mesh) * tip_dir;
}

Tip * AnastamosisPDD::BestInElement(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh, vec3d tip_origin, vec3d tip_dir, int exclude, double& best_distance)
{
	best_distance = std::numeric_limits<double>::max();
	int best_index = -1;
	for(int i=0; i < angio_element->grown_segments.size();i++)
	{
		Tip * tip = angio_element->grown_segments[i]->front;
		if(tip->initial_fragment_id != exclude)
		{
			vec3d pos = tip->GetPosition(mesh);
			double dist = (tip_origin - pos).norm2();
			if (dist < best_distance)
			{
				best_distance = dist;
				best_index = i;
			}
		}
	}
	if(best_index != -1)
	{
		best_distance = sqrt(best_distance);
		return angio_element->grown_segments[best_index]->front;
	}
	return nullptr;
}


BEGIN_PARAMETER_LIST(AnastamosisPDD, PositionDependentDirection)
ADD_PARAMETER(anastamosis_radius, FE_PARAM_DOUBLE, "anastamosis_radius");
ADD_PARAMETER(contribution, FE_PARAM_DOUBLE, "contribution");
ADD_PARAMETER(fuse_radius, FE_PARAM_DOUBLE, "fuse_radius");
ADD_PARAMETER(fuse_angle, FE_PARAM_DOUBLE, "fuse_angle");
END_PARAMETER_LIST();

void FiberPDD::Update(FEMesh * mesh, FEAngio* angio)
{

}

vec3d FiberPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continute_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio)
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
	return mix3d(prev, fiber_direction, contribution);
	// return mix(prev, fiber_direction, contribution);
}