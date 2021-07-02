#include "PositionDependentDirection.h"
#include "angio3d.h"
#include "Tip.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FESolidDomain.h>
#include "FEAngioMaterialPoint.h"
#include "FEAngio.h"
#include "FEAngioMaterial.h"
#include "FECore/FEDomainMap.h"
#include <iostream>
#include <algorithm>
#include <FECore/mathalg.h>
#include <FECore/FEDomainMap.h>

BEGIN_FECORE_CLASS(PositionDependentDirection, FEMaterial)
ADD_PARAMETER(contribution, "contribution");
END_FECORE_CLASS();

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
	mat3ds eye;
}

/*bool sortinrev(const pair<double, int> &a, const pair<double, int> &b)
{
	return (a.first > b.first);
}*/

BEGIN_FECORE_CLASS(LaGrangePStrainPDD, PositionDependentDirection)
ADD_PARAMETER(beta, "beta");
END_FECORE_CLASS();

struct EigComp {
	bool operator()(const std::pair<double, vec3d>& x , std::pair<double, vec3d>& y) const{
		if (std::abs(x.first) > std::abs(y.first)) return true;
		else return false;
	}
};


// PDD that determines the principal strain of least magnitude and mixes it with the previous direction
// TODO: Determine contribution and direction based on ratio of each direction's magnitude.
vec3d LaGrangePStrainPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio) 
{
	// for each integration point
	//std::vector<double> FirstPSMag;
	//std::vector<double> SecondPSMag;
	//std::vector<double> ThirdPSMag;
	//std::vector<vec3d> FirstPSDir;
	//std::vector<vec3d> SecondPSDir;
	mat3ds E;
	//vec3d MagVecs[3];

	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++) {
		FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		FEElasticMaterialPoint* elastic_mp = gauss_point->ExtractData<FEElasticMaterialPoint>();
		E += elastic_mp->Strain();

	}
	E = E / (angio_element->_elem->GaussPoints());
	// get eigenvalues and eigenvectors
	double d[3]; vec3d r[3];
	E.eigen(d,r);
	// Create vector of pairs of eigenvalues and eigenvectors
	std::vector<std::pair<double, vec3d>> v = { {d[0],r[0]}, {d[1],r[1]}, {d[2],r[2]} };
	// Sort based on absolute value of eigenvalues
	std::sort(v.begin(), v.end(), EigComp());
	// multiply magnitudes by directions and then determine the scaling based on this value
	//MagVecs[0] = v[0].second*v[0].first; MagVecs[1] = v[1].second*v[1].first; MagVecs[2] = v[2].second*v[2].first;
	vec3d least_strain = v[2].second;
	least_strain.unit();
	return angio_element->_angio_mat->mix_method->ApplyMixAxis(prev, least_strain, beta);
}

BEGIN_FECORE_CLASS(ECMDensityGradientPDD, PositionDependentDirection)
ADD_PARAMETER(threshold, "threshold");
ADD_PARAMETER(alpha_override, "alpha_override");
END_FECORE_CLASS();

vec3d ECMDensityGradientPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio)
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
		vec3d new_dir = angio_element->_angio_mat->mix_method->ApplyMix(prev, perpendicularToGradient, contribution);
		if (prev* perpendicularToGradient < 0) { new_dir = -new_dir; }
		return new_dir;
		//return mix3d(prev, perpendicularToGradient, contribution);
	}
	return prev;
}

BEGIN_FECORE_CLASS(RepulsePDD, PositionDependentDirection)
ADD_PARAMETER(threshold, "threshold");
ADD_PARAMETER(alpha_override, "alpha_override");
ADD_PARAMETER(grad_threshold, "grad_threshold");
END_FECORE_CLASS();

vec3d RepulsePDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio)
{
	std::vector<double> repulse_at_integration_points;

	// find which angio elements have repulse values
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		FEElasticMaterialPoint* elastic_mp = gauss_point->ExtractData<FEElasticMaterialPoint>();

		// Get the replse value and scale it by any deformation
		repulse_at_integration_points.push_back(angio_mp->repulse_value * (1.0 / elastic_mp->m_J));
	}
	// if the threshold method should be used
	if(grad_threshold)
	{
		// get the gradient for the element from the repulse value at the integration points from the current position
		vec3d grad = pangio->gradient(angio_element->_elem, repulse_at_integration_points, local_pos);
		// get the magnitude and compare to the threshold
		double gradnorm = grad.norm();
		//std::cout << "gradnorm is " << gradnorm << endl;
		if (gradnorm > threshold)
		{	
			grad = -grad;
			if (alpha_override)
			{
				alpha = contribution;
			}
			// use the negative of the gradient of the repulse values to tell vessels which direction to grow away from
			vec3d new_dir = angio_element->_angio_mat->mix_method->ApplyMix(prev, grad, alpha);
			if (prev* grad < 0) { new_dir = -new_dir; }
			return new_dir;
			//return angio_element->_angio_mat->mix_method->ApplyMix(prev, grad, contribution);
			//return mix(prev,grad, contribution);
		}	
	}
	// if the gradient method is not being used
	else
	{
		double nodal_repulse[FESolidElement::MAX_NODES];
		double H[FESolidElement::MAX_NODES];
		// project the repulse value from integration points to the nodes
		angio_element->_elem->project_to_nodes(&repulse_at_integration_points[0],nodal_repulse);
		double val = 0;
		// determine shape function value for the local position
		angio_element->_elem->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
		// for each node in an angio element determine the contributions from the shape functions
		for(int i=0; i < angio_element->_elem->Nodes();i++)
		{
			val += H[i] * nodal_repulse[i];
		}
		//SL: moving this into the loop below so it is only calculated when necessary 
		// vec3d grad = pangio->gradient(angio_element->_elem, repulse_at_integration_points, local_pos);
		// if the value from the shape functions is greater than the threshold then the vessel will be repulsed
		if(val > threshold)
		{
			// determine the gradient in each angio element of the repulse value at a given local position
			vec3d grad = pangio->gradient(angio_element->_elem, repulse_at_integration_points, local_pos);
			// new direction will be in the opposite of the other direction.
			grad = -grad;
			// if alpha_override is set to 1 then set alpha to 1. Seems like an unneccessary step
			if (alpha_override)
			{
				alpha = contribution;
			}
			// new vector direction is mix of previous and deflected direction.
			vec3d new_dir = angio_element->_angio_mat->mix_method->ApplyMix(prev, grad, 1);
			// if previous direction dotted with grad is -ve i.e. the angle between them is greater than 90 degrees. not sure if this makes sense. 
			// if (prev* grad < 0) { new_dir = -new_dir; }
			return new_dir;
			//return angio_element->_angio_mat->mix_method->ApplyMix(prev, grad, contribution);
			//return mix(prev, grad, contribution);
		}
	}
	return prev;
	
}

BEGIN_FECORE_CLASS(ConcentrationGradientPDD, PositionDependentDirection)
ADD_PARAMETER(threshold, "threshold");
ADD_PARAMETER(alpha_override, "alpha_override");
ADD_PARAMETER(sol_id, "sol_id");
END_FECORE_CLASS();

vec3d ConcentrationGradientPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio)
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
		vec3d new_dir = angio_element->_angio_mat->mix_method->ApplyMix(prev, grad, contribution);
		if (prev* grad< 0) { new_dir = -new_dir; }
		return new_dir;
		//return angio_element->_angio_mat->mix_method->ApplyMix(prev, grad, contribution);
		//return mix(prev, grad, contribution);
	}
	return prev;
}

//used in the anastamosis modifier
double AnastamosisPDD::distance2(FESolidElement * se, vec3d local_pos, Tip * tip, FEMesh* mesh)
{
	vec3d pos;
	double H[FESolidElement::MAX_NODES];
	se->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
	for(int i=0; i < se->Nodes();i++)
	{
		pos += (mesh->Node(se->m_node[i]).m_rt * H[i]);
	}
	return (tip->GetPosition(mesh) - pos).norm2();
}

vec3d AnastamosisPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio)
{
	vec3d tip_pos;
	double H[FESolidElement::MAX_NODES];
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
			continue_growth = false;
		}
		vec3d new_dir = angio_element->_angio_mat->mix_method->ApplyMix(prev, dir, contribution);
		if (prev* dir< 0) { new_dir = -new_dir; }
		return new_dir;
		//return angio_element->_angio_mat->mix_method->ApplyMix(prev, dir, contribution);
		//return mix(prev, dir, contribution);
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


BEGIN_FECORE_CLASS(AnastamosisPDD, PositionDependentDirection)
ADD_PARAMETER(anastamosis_radius, "anastamosis_radius");
ADD_PARAMETER(contribution, "contribution");
ADD_PARAMETER(fuse_radius, "fuse_radius");
ADD_PARAMETER(fuse_angle, "fuse_angle");
END_FECORE_CLASS();

void FiberPDD::Update(FEMesh * mesh, FEAngio* angio)
{

}

vec3d FiberPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio)
{
	std::vector<quatd> gauss_data;
	std::vector<quatd> nodal_data;
	//FEElasticMaterial* pmm = dynamic_cast<FEElasticMaterial*> (angio_element->_mat->GetElasticMaterial);
	//std::cout << "material is " << angio_element->_angio_mat->GetMatrixMaterial()->GetElasticMaterial()->FindProperty("m_Base") << endl;
	for (int i = 0; i< angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		vec3d axis(1, 0, 0);

		//Ask Steve about this
		//get the FE domain
		FEDomain* Dom = dynamic_cast<FEDomain*>(angio_element->_elem->GetMeshPartition());
		//
		FEElementSet* elset = mesh->FindElementSet(Dom->GetName());
		int local_index = elset->GetLocalIndex(*angio_element->_elem);

		FEAngioMaterial* Mat_ang = Dom->GetMaterial()->ExtractProperty<FEAngioMaterial>();
		FEMaterial * Mat_a = Mat_ang->GetMatrixMaterial();
		// assumes that materials mat_axis is already mapped which we'll need to do somewhere else.
		FEParam* matax = Mat_a->FindParameter("mat_axis");
		FEParamMat3d& p = matax->value<FEParamMat3d>();
		FEMappedValueMat3d* val = dynamic_cast<FEMappedValueMat3d*>(p.valuator());
		FEDomainMap* map = dynamic_cast<FEDomainMap*>(val->dataMap());
		//mat3d m_Q = Mat_a->GetLocalCS(*mp);

		//mat3d m_Q = map->valueMat3d(emp);

		//axis = emp->m_F * m_Q * axis;
		gauss_data.push_back({ axis });

		/*axis = emp->m_F * emp->m_Q * axis;
		gauss_data.push_back({ axis });*/
	}

	quatd rv = interpolation_prop->Interpolate(angio_element->_elem, gauss_data, local_pos, mesh);
	vec3d fiber_direction = rv.GetVector();
	return angio_element->_angio_mat->mix_method->ApplyMixAxis(tip_dir, fiber_direction, contribution);
	// return mix(prev, fiber_direction, contribution);
}

void FractionalAnisotropyPDD::Update(FEMesh * mesh, FEAngio* angio)
{

}

vec3d FractionalAnisotropyPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio)
{
	angio_element->UpdateSPD();
	angio_element->UpdateAngioFractionalAnisotropy();

	mat3d ax;
	ax.setCol(0, vec3d(angio_element->angioSPD.xx(), angio_element->angioSPD.xy(), angio_element->angioSPD.xz()));
	ax.setCol(1, vec3d(angio_element->angioSPD.xy(), angio_element->angioSPD.yy(), angio_element->angioSPD.yz()));
	ax.setCol(2, vec3d(angio_element->angioSPD.xz(), angio_element->angioSPD.yz(), angio_element->angioSPD.zz()));
	
	// get the vectors of the principal directions and sort in descending order
	std::vector<pair<double, int>> v;
	v.push_back(pair<double, int>(ax.col(0).norm(), 0));
	v.push_back(pair<double, int>(ax.col(1).norm(), 1));
	v.push_back(pair<double, int>(ax.col(2).norm(), 2));
	sort(v.begin(), v.end(), sortinrev);

	// store the indices
	int i = v[0].second;
	int j = v[1].second;
	int k = v[2].second;

	vec3d axis_0 = ax.col(i); axis_0.unit();
	vec3d axis_1 = ax.col(j); axis_1.unit();
	vec3d axis_2 = ax.col(k); axis_2.unit();
	double r0 = (ax.col(i).norm());
	double r1 = (ax.col(j).norm());
	double r2 = (ax.col(k).norm());

	double theta_12 = angio_element->GetEllipseAngle(r0, r1,-PI/2,PI,180);
	double theta_13 = angio_element->GetEllipseAngle(r0, r2,-PI/2,PI,180);

	// rotate the primary direction by theta_12 about the normal between them
	vec3d axis = mix3d_t(axis_0, axis_1, theta_12); axis.unit();
	vec3d fiber_dir = mix3d_t(axis, axis_2, theta_13); fiber_dir.unit();
	// determine the fiber directon contribution from the fractional anisotropy
	/*if (alpha_override)
	{
		alpha = contribution;
	}*/
	//double FA_cont = std::min(angio_element->angioFA, 0.36);
	//std::cout << "contribution is " << contribution << endl;
	alpha = std::min((contribution*(((0.75 - 0.36) / (0.6 - 0))*angio_element->angioFA + 0.36)),0.75);
	/*if (angio_element->angioFA > 0.5) 
	{ 
		alpha = contribution*(((1 - 0.36) / (0.75 - 0))*FA_cont + 0.36); 
	}*/
	return angio_element->_angio_mat->mix_method->ApplyMixAxis(tip_dir, fiber_dir, alpha);
}

BEGIN_FECORE_CLASS(FractionalAnisotropyPDD,PositionDependentDirection)
ADD_PARAMETER(alpha_override, "alpha_override");
END_FECORE_CLASS();

void FractionalAnisotropyMatPointPDD::Update(FEMesh * mesh, FEAngio* angio)
{

}

vec3d FractionalAnisotropyMatPointPDD::ApplyModifiers(vec3d prev, AngioElement* angio_element, vec3d local_pos, int initial_fragment_id, int current_buffer, double& alpha, bool& continue_growth, vec3d& tip_dir, FEMesh* mesh, FEAngio* pangio)
{
	// vector containing the SPD for each gauss point in the element
	std::vector<mat3ds> SPDs_gausspts;

	// get each gauss point's SPD
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		// get the angio point
		FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		//std::cout << angio_mp->angioSPD.xx() << ", " << angio_mp->angioSPD.yy() << ", " << angio_mp->angioSPD.zz() << ", " << angio_mp->angioSPD.xy() << ", " << angio_mp->angioSPD.yz() << ", " << angio_mp->angioSPD.xz() << endl;
		// Get the SPD
		//angio_mp->nhit = 1;
		angio_mp->UpdateSPD();
		//std::cout << angio_mp->angioSPD.xx() << ", " << angio_mp->angioSPD.yy() << ", " << angio_mp->angioSPD.zz() << ", " << angio_mp->angioSPD.xy() << ", " << angio_mp->angioSPD.yz() << ", " << angio_mp->angioSPD.xz() << endl;
		SPDs_gausspts.push_back(angio_mp->angioSPD);
	}
	// array containing the SPD for each node in the element
	mat3ds SPDs_nodes[FESolidElement::MAX_NODES];
	// array for the shape function values
	double H[FESolidElement::MAX_NODES];
	// project the spds from integration points to the nodes
	angio_element->_elem->project_to_nodes(&SPDs_gausspts[0], SPDs_nodes);
	// determine shape function value for the local position
	angio_element->_elem->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
	//angio_element->_elem->shape_fnc(H, 0, 0, 0);
	// Get the interpolated SPD from the shape function-weighted Average Structure Tensor
	mat3ds SPD_int = weightedAverageStructureTensor(SPDs_nodes, H, angio_element->_elem->Nodes());
	// get the vectors of the principal directions and sort in descending order
	std::vector<pair<double, int>> v;
	mat3d ax;
	ax.setCol(0, vec3d(SPD_int.xx(), SPD_int.xy(), SPD_int.xz()));
	ax.setCol(1, vec3d(SPD_int.xy(), SPD_int.yy(), SPD_int.yz()));
	ax.setCol(2, vec3d(SPD_int.xz(), SPD_int.yz(), SPD_int.zz()));
	v.push_back(pair<double, int>(ax.col(0).norm(), 0));
	v.push_back(pair<double, int>(ax.col(1).norm(), 1));
	v.push_back(pair<double, int>(ax.col(2).norm(), 2));
	sort(v.begin(), v.end(), sortinrev);

	// store the indices
	int i = v[0].second;
	int j = v[1].second;
	int k = v[2].second;

	vec3d axis_0 = ax.col(i); axis_0.unit();
	vec3d axis_1 = ax.col(j); axis_1.unit();
	vec3d axis_2 = ax.col(k); axis_2.unit();
	double r0 = (ax.col(i).norm());
	double r1 = (ax.col(j).norm());
	double r2 = (ax.col(k).norm());
	double theta_12 = angio_element->GetEllipseAngle(r0, r1, -PI / 2, PI, 180);
	double theta_13 = angio_element->GetEllipseAngle(r0, r2, -PI / 2, PI, 180);

	// rotate the primary direction by theta_12 about the normal between them
	vec3d axis = mix3d_t(axis_0, axis_1, theta_12); axis.unit();
	vec3d fiber_dir = mix3d_t(axis, axis_2, theta_13); fiber_dir.unit();
	// determine the fiber directon contribution from the fractional anisotropy
	/*if (alpha_override)
	{
	alpha = contribution;
	}*/
	//double FA_cont = std::min(angio_element->angioFA, 0.36);
	//std::cout << "contribution is " << contribution << endl;
	
	// calculate the fractional anisotropy
	//double angioFA_int = sqrt(0.5)*(sqrt(pow(r0 - r1, 2) + pow(r1 - r2, 2) + pow(r2 - r0, 2)) / (sqrt(pow(r0, 2) + pow(r1, 2) + pow(r2, 2))));
	//double angioFA_int = 1.0 - (r1 / r0);
	//alpha = std::min(angioFA_int, 0.75);
	//alpha = std::min((contribution*(((0.64 - 0.36) / (0.64 - 0))*angioFA_int + 0.36)), 0.64);
	//alpha = 1.0;
	alpha = contribution;
	return angio_element->_angio_mat->mix_method->ApplyMixAxis(tip_dir, fiber_dir, alpha);
}

BEGIN_FECORE_CLASS(FractionalAnisotropyMatPointPDD, PositionDependentDirection)
ADD_PARAMETER(alpha_override, "alpha_override");
END_FECORE_CLASS();