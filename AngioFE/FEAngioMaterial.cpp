#include "StdAfx.h"
#include "FEAngioMaterial.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include "FEBioMech/FEElasticMixture.h"
#include "FEBioMech/FEFiberMaterialPoint.h"
#include "FECore/FEElementTraits.h"
#include "FEBioMech/FEViscoElasticMaterial.h"
#include "FECore/FESPRProjection.h"
#include <iostream>
#include <algorithm>
#include "angio3d.h"
#include "AngioElement.h"
#include "Segment.h"
#include "Tip.h"
#include "FECell.h"
#include "VariableInterpolation.h"


//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEAngioMaterial, FEElasticMaterial)
//ADD_PARAMETER(initial_segment_velocity, "initial_segment_velocity");
ADD_PARAMETER(vessel_radius, "vessel_radius");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEAngioMaterial::FEAngioMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	AddClassProperty(this, &common_properties, "common_properties");
	AddClassProperty(this, &matrix_material, "matrix");
	AddClassProperty(this, &angio_stress_policy, "angio_stress_policy");

	AddClassProperty(this, &proto_pdd_manager, "proto_pdd_manager");
	AddClassProperty(this, &pdd_manager, "pdd_manager");
	AddClassProperty(this, &proto_psc_manager, "proto_psc_manager");
	AddClassProperty(this, &psc_manager, "psc_manager");
	AddClassProperty(this, &proto_cm_manager, "proto_cm_manager");
	AddClassProperty(this, &cm_manager, "cm_manager");
	AddClassProperty(this, &mix_method, "mix_method");
	AddClassProperty(this, &velocity_manager, "velocity_manager");
	AddClassProperty(this, &im_manager, "im_manager", FEProperty::Optional);
	AddClassProperty(this, &nodedata_interpolation_manager, "nodedata_interpolation_manager", FEProperty::Optional);
	AddClassProperty(this, &branch_policy, "branch_policy", FEProperty::Optional);
	AddClassProperty(this, &proto_branch_policy, "proto_branch_policy", FEProperty::Optional);
	AddClassProperty(this, &cell_species_manager, "cell_species_manager", FEProperty::Optional);
	AddClassProperty(this, &cell_reaction_manager, "cell_reaction_manager", FEProperty::Optional);
}

FEAngioMaterial::~FEAngioMaterial()
{
}

//-----------------------------------------------------------------------------
bool FEAngioMaterial::Init()
{
	// Make sure the angio class was allocated
	if (m_pangio == nullptr) return false;

	if(!matrix_material->Init()) return false;

	if (FEElasticMaterial::Init() == false) return false;

	if(!common_properties->vessel_material->Init()) return false;

	
	if(proto_branch_policy && branch_policy)
	{
		if(proto_branch_policy->GetSuperClassID() != branch_policy->GetSuperClassID())
		{
			//the type of the proto branch policy must match the branch policy
			return false;
		}
	}
	if(proto_branch_policy && !branch_policy)
	{
		//a branch policy must be defineded if a proto branch policy is defined
		return false;
	}

	//culture must be initialized here  so pangio is defined
	assert(m_pangio);



	// add the user sprouts
	std::vector<int> matls;
	
	FEMesh& mesh = GetFEModel()->GetMesh();
	
	FECoreBase * base = GetParent();
	if(base)
	{
		matls.emplace_back(base->GetID());
	}
	else
	{
		matls.emplace_back(this->GetID_ang());
	}



	return true;
}

//-----------------------------------------------------------------------------
mat3ds FEAngioMaterial::AngioStress(FEAngioMaterialPoint& angioPt)
{
	mat3ds val;
	val.zero();
	return val;
}

double FEAngioMaterial::FindDensityScale(FEAngioMaterialPoint * mp)
{
	//!SL check this.
	return 1.0;
}

mat3ds FEAngioMaterial::Stress(FEMaterialPoint& mp)
{
	FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
	FEElasticMaterialPoint& elastic_pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3ds s;
	s.zero();
	//should always be true but we should check
	assert(angioPt);
	if(angioPt)
	{
		FEElasticMaterialPoint& vessel_elastic = *angioPt->vessPt->ExtractData<FEElasticMaterialPoint>();
		FEElasticMaterialPoint& matrix_elastic = *angioPt->matPt->ExtractData<FEElasticMaterialPoint>();

		vessel_elastic.m_rt = elastic_pt.m_rt; //spatial position
		vessel_elastic.m_r0 = elastic_pt.m_r0; //material position
		vessel_elastic.m_F = elastic_pt.m_F; //deformation gradient
		vessel_elastic.m_J = elastic_pt.m_J; //determinant
		
		matrix_elastic.m_rt = elastic_pt.m_rt;
		matrix_elastic.m_r0 = elastic_pt.m_r0;
		matrix_elastic.m_F = elastic_pt.m_F;
		matrix_elastic.m_J = elastic_pt.m_J;

		vessel_elastic.m_s = common_properties->vessel_material->Stress(*(angioPt->vessPt));
		matrix_elastic.m_s = matrix_material->Stress(*(angioPt->matPt));

		s = angioPt->m_as + angioPt->vessel_weight*vessel_elastic.m_s + angioPt->matrix_weight*matrix_elastic.m_s;
	}
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEAngioMaterial::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& elastic_pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
	tens4ds s(0.0);
	if(angioPt)
	{
		FEElasticMaterialPoint& vessel_elastic = *angioPt->vessPt->ExtractData<FEElasticMaterialPoint>();
		vessel_elastic.m_rt = elastic_pt.m_rt;
		vessel_elastic.m_r0 = elastic_pt.m_r0;
		vessel_elastic.m_F = elastic_pt.m_F;
		vessel_elastic.m_J = elastic_pt.m_J;
		FEElasticMaterialPoint& matrix_elastic = *angioPt->matPt->ExtractData<FEElasticMaterialPoint>();
		matrix_elastic.m_rt = elastic_pt.m_rt;
		matrix_elastic.m_r0 = elastic_pt.m_r0;
		matrix_elastic.m_F = elastic_pt.m_F;
		matrix_elastic.m_J = elastic_pt.m_J;
		s = angioPt->vessel_weight*common_properties->vessel_material->Tangent(*angioPt->vessPt) + angioPt->matrix_weight*matrix_material->Tangent(*angioPt->matPt);
	}
	return s;
}

double FEAngioMaterial::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& elastic_pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
    
	// calculate strain energy density
	double sed = 0.0;
	if(angioPt)
	{
		FEElasticMaterialPoint& vessel_elastic = *angioPt->vessPt->ExtractData<FEElasticMaterialPoint>();
		vessel_elastic.m_rt = elastic_pt.m_rt;
		vessel_elastic.m_r0 = elastic_pt.m_r0;
		vessel_elastic.m_F = elastic_pt.m_F;
		vessel_elastic.m_J = elastic_pt.m_J;
		FEElasticMaterialPoint& matrix_elastic = *angioPt->matPt->ExtractData<FEElasticMaterialPoint>();
		matrix_elastic.m_rt = elastic_pt.m_rt;
		matrix_elastic.m_r0 = elastic_pt.m_r0;
		matrix_elastic.m_F = elastic_pt.m_F;
		matrix_elastic.m_J = elastic_pt.m_J;
		sed = angioPt->vessel_weight*common_properties->vessel_material->ExtractProperty<FEElasticMaterial>()->StrainEnergyDensity(*angioPt->vessPt) + angioPt->matrix_weight*matrix_material->ExtractProperty<FEElasticMaterial>()->StrainEnergyDensity(*angioPt->matPt);
	}
	return sed;
}

vec3d FEAngioMaterial::CheckFaceProximity(vec3d pos, vec3d dir) {
	double val = 0.99; 
	if ((pos.x > val) & (dir.x >= 0)) { pos.x = val; }
	else if ((pos.x < -val) & (dir.x <= 0)) { pos.x = -val; }
	if ((pos.y > val) & (dir.y >= 0)) { pos.y = val; }
	else if ((pos.y < -val) & (dir.y <= 0)){ pos.y = -val; }
	if ((pos.z > val) & (dir.z >= 0)) { pos.z = val; }
	else if ((pos.z < -val) & (dir.z <= 0)) { pos.z = -val; }
	return pos;
}

std::vector<double> FEAngioMaterial::CornerPossibleValues(vec3d local_pos, vec3d nat_dir) {
	//! need to update this in the future for other element types.
	std::vector<double> possible_values;
	vec3d upper_bounds = vec3d(1, 1, 1);
	vec3d lower_bounds = vec3d(-1, -1, -1);

	if (nat_dir.x > 0)
	{
		double temp = (upper_bounds.x - local_pos.x) / nat_dir.x;
		if (temp >= 0)
			possible_values.push_back(temp);
	}
	else if (nat_dir.x < 0)
	{
		double temp = (lower_bounds.x - local_pos.x) / nat_dir.x;
		if (temp >= 0)
			possible_values.push_back(temp);
	}
	if (nat_dir.y > 0)
	{
		double temp = (upper_bounds.y - local_pos.y) / nat_dir.y;
		if (temp >= 0)
			possible_values.push_back(temp);
	}
	else if (nat_dir.y < 0)
	{
		double temp = (lower_bounds.y - local_pos.y) / nat_dir.y;
		if (temp >= 0)
			possible_values.push_back(temp);
	}
	if (nat_dir.z > 0)
	{
		double temp = (upper_bounds.z - local_pos.z) / nat_dir.z;
		if (temp >= 0)
			possible_values.push_back(temp);
	}
	else if (nat_dir.z < 0)
	{
		double temp = (lower_bounds.z - local_pos.z) / nat_dir.z;
		if (temp >= 0)
			possible_values.push_back(temp);
	}
	return possible_values;
}

void FEAngioMaterial::SetSeeds(AngioElement* angio_elem)
{
	int offset = angio_elem->_elem->GetID();
	static int seed = int (m_pangio->m_fem->GetGlobalConstant("seed"));
	std::seed_seq sq{ offset + seed };
	angio_elem->_rengine.seed(sq);
	//offset++;
}

double FEAngioMaterial::GetSegmentVelocity(AngioElement * angio_element, vec3d local_pos, double time_shift, FEMesh* mesh)
{
	double const vel_init = 1;
	return velocity_manager->ApplyModifiers(vel_init, local_pos, angio_element, time_shift, mesh);
}

double FEAngioMaterial::GetMin_dt(AngioElement* angio_elem, FEMesh* mesh)
{
	assert(angio_elem);
	assert(angio_elem->_elem);

	// Calculate the maximum growth velocity based on the potential velocities calculated at each of the Gauss points.
	double max_grow_velocity = 0.0;
	for (int i = 0; i < angio_elem->_elem->GaussPoints();i++)
	{
		vec3d pos(angio_elem->_elem->gr(i), angio_elem->_elem->gs(i), angio_elem->_elem->gt(i));
		double vel = angio_elem->_angio_mat->GetSegmentVelocity(angio_elem, pos, 0.0, mesh);
		max_grow_velocity = std::max(max_grow_velocity, vel);
	}

	double min_side_length_so_far = std::numeric_limits<double>::max();

	for (int i = 0; i < angio_elem->_elem->GaussPoints(); i++)
	{
		double Gr[FESolidElement::MAX_NODES];
		double Gs[FESolidElement::MAX_NODES];
		double Gt[FESolidElement::MAX_NODES];
		angio_elem->_elem->shape_deriv(Gr, Gs, Gt, angio_elem->_elem->gr(i), angio_elem->_elem->gs(i), angio_elem->_elem->gt(i));

		// Calculate the length of each side of the element.
		vec3d er, es, et; // Basis vectors of the natural coordinates
		for (int j = 0; j < angio_elem->_elem->Nodes(); j++)
		{
			er += mesh->Node(angio_elem->_elem->m_node[j]).m_rt * Gr[j];
		}
		for (int j = 0; j < angio_elem->_elem->Nodes(); j++)
		{
			es += mesh->Node(angio_elem->_elem->m_node[j]).m_rt * Gs[j];
		}
		for (int j = 0; j < angio_elem->_elem->Nodes(); j++)
		{
			et += mesh->Node(angio_elem->_elem->m_node[j]).m_rt * Gt[j];
		}

		min_side_length_so_far = std::min(min_side_length_so_far, er.norm2()); // norm2() is magnitude without sqrt.
		min_side_length_so_far = std::min(min_side_length_so_far, es.norm2());
		min_side_length_so_far = std::min(min_side_length_so_far, et.norm2());
	}

	min_side_length_so_far = sqrt(min_side_length_so_far);

	// Basis vectors already represented correctly for HEX8 elements, need to take half of this value for TET4's.
	switch(angio_elem->_elem->Type())
	{
		case FE_Element_Type::FE_TET4G1:
		case FE_Element_Type::FE_TET4G4:
		case FE_Element_Type::FE_TET10G1:
		case FE_Element_Type::FE_TET10G4:
		case FE_Element_Type::FE_TET10G4RI1:
		case FE_Element_Type::FE_TET10G8:
		case FE_Element_Type::FE_TET10G8RI4:
		case FE_Element_Type::FE_TET10GL11:
			min_side_length_so_far /= 2;
		default:
			break;
	}

	// Return dt = d / v.
	return (min_side_length_so_far * 0.5 * dt_safety_multiplier) / max_grow_velocity;
}

void FEAngioMaterial::GrowSegments(AngioElement * angio_elem, double end_time, int buffer_index, double min_scale_factor, double bounds_tolerance)
{
	assert(angio_elem);

	// Grow tips that originate in this element and have not reached an adjacent face.
	// create a vector for the active tips in the element
	std::vector<Tip*> & tips = angio_elem->active_tips[buffer_index].at(angio_elem);
	// for each tip ensure that the tip belongs to the element being analyzed. Next call growth in the element.
	for(int j=0; j < tips.size();j++)
	{
		assert(tips[j]->angio_element == angio_elem);
		// call GrowthInElement providing the end time, the active tip, -1 to indicate the tip originated within this element, the buffer index, etc.
		GrowthInElement(end_time, tips[j], -1, buffer_index, min_scale_factor,bounds_tolerance);
	}
}

void FEAngioMaterial::ProtoGrowSegments(AngioElement * angio_elem, double end_time, int buffer_index, double min_scale_factor, double bounds_tolerance)
{
	assert(angio_elem);

	// copy the vector of active tips in the current angio element
	std::vector<Tip*> & tips = angio_elem->active_tips[buffer_index].at(angio_elem);\
	// for each tip in the current element grow the tips.
	for (int j = 0; j < tips.size(); j++)
	{
		assert(tips[j]->angio_element == angio_elem);
		ProtoGrowthInElement(end_time, tips[j], -1, buffer_index, min_scale_factor, bounds_tolerance);
	}
}

//Grow a provided active tip based on the dt, originating element, etc.
void FEAngioMaterial::GrowthInElement(double end_time, Tip * active_tip, int source_index, int buffer_index, double min_scale_factor, double bounds_tolerance)
{
	////////// Initialization //////////

	// get the current angio element and make sure it's active
	auto angio_element = active_tip->angio_element;
	
	assert(active_tip);
	// get the mesh, set the min segment length, and set the next buffer (change to write buffer from read).
	auto mesh = m_pangio->GetMesh();
	int next_buffer_index = (buffer_index + 1) % 2;
	// calculate the change in time
	double dt = end_time - active_tip->time;
	if (dt < 0) { return; }
	// determine the length the segment grows.
	double grow_vel = this->GetSegmentVelocity(angio_element, active_tip->GetLocalPosition(), active_tip->seg_time, mesh);
	//If the velocity is negative wrt the tip direction then end the tip/
	if(grow_vel < 0)
	{
		return;
	}
	double grow_len = grow_vel*dt;

	////////// Determine the growth length and new position //////////

	// get the local position of the tip
	double Gr[FESolidElement::MAX_NODES];
	double Gs[FESolidElement::MAX_NODES];
	double Gt[FESolidElement::MAX_NODES];
	vec3d local_pos = active_tip->GetLocalPosition();
	// get the shape function derivative values for the element type
	angio_element->_elem->shape_deriv(Gr, Gs, Gt, local_pos.x, local_pos.y, local_pos.z);
	vec3d er, es, et; // Basis vectors of the natural coordinates

	// for each node determine the contribution to the conversion between the natural and global coordinate systems.
	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		er += mesh->Node(angio_element->_elem->m_node[j]).m_rt* Gr[j];
	}
	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		es += mesh->Node(angio_element->_elem->m_node[j]).m_rt* Gs[j];
	}
	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		et += mesh->Node(angio_element->_elem->m_node[j]).m_rt* Gt[j];
	}

	// convert natural coords to global
	mat3d natc_to_global(er, es, et);
	mat3d global_to_natc = natc_to_global.inverse();
	bool continue_growth = true;

	// get the growth direction
	double alpha = cm_manager->ApplyModifiers(dt, active_tip->angio_element, active_tip->GetLocalPosition(), mesh);
	// determine contributions from previous segments
	vec3d psc_dir = psc_manager->ApplyModifiers(vec3d(1, 0, 0), active_tip->angio_element, active_tip->GetLocalPosition(), active_tip->GetDirection(mesh), mesh);
	// determine contributions from local stimuli
	vec3d pdd_dir = pdd_manager->ApplyModifiers(vec3d(1, 0, 0), active_tip->angio_element, active_tip->GetLocalPosition() , active_tip->initial_fragment_id , buffer_index , continue_growth, psc_dir, alpha, mesh, m_pangio);
	if(!continue_growth)
	{
		return;
	}
	// get the new direction in global coordinates
	vec3d global_dir = mix_method->ApplyMix(psc_dir, pdd_dir, alpha);
	global_dir.unit();
	vec3d nat_dir = global_to_natc* global_dir;
	
	////////// Calculate variables to determine which element to grow into //////////
	
	double factor;         
	// Determine the scale factor needed for the vessel to reach the face of an element
	bool proj_success = m_pangio->ScaleFactorToProjectToNaturalCoordinates(active_tip->angio_element->_elem, nat_dir, local_pos, factor, min_scale_factor); 
	// calculate the position the vessel will grow to in the natural coords
	vec3d possible_natc = local_pos + (nat_dir*factor);
	// determine the hypothetical length to grow
	double possible_grow_length = m_pangio->InElementLength(active_tip->angio_element->_elem, local_pos, possible_natc);
	
	////////// Determine which element to grow into //////////

	// if the vessel is going to stay in this element:
	if ((possible_grow_length > grow_len) && proj_success)
	{
		#ifndef NDEBUG
			#pragma omp critical
			std::cout << "case 0: growth within element" << endl;
		#endif
		//determine the new natural coordinate position in this element. Not sure that this calculation needs to be this convoluted.
		vec3d real_natc = local_pos + (nat_dir * (factor * (grow_len / possible_grow_length)));
		
		//create the new tip
		Tip * next = new Tip(active_tip, mesh);
		next->time = end_time;
		Segment * seg = new Segment();
		next->SetLocalPosition(real_natc, mesh);
		next->TipCell->time = next->time;
		next->TipCell->UpdateSpecies(mesh);
		next->parent = seg;
		assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G4 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);
		assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G1 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);

		seg->back = active_tip;
		seg->front = next;
		seg->is_branch_base = seg->back->is_branch;
		if (active_tip->parent)
		{
			seg->parent = active_tip->parent;
		}
		if(angio_element->_angio_mat->branch_policy && angio_element->vessel_weight < angio_element->_angio_mat->thresh_vess_weight){
			if (seg->is_branch_base) { seg->LengthSinceBranch = 0.0; seg->time_emerged = seg->back->time; }
			else {
				seg->LengthSinceBranch = seg->parent->LengthSinceBranch + seg->RefLength(mesh);
			}
			double l2b = dynamic_cast<DelayBranchInfo*>(angio_element->branch_info)->length_to_branch;
			if (seg->LengthSinceBranch >= l2b) {
				double rel = (l2b - seg->parent->LengthSinceBranch) / (seg->LengthSinceBranch - seg->parent->LengthSinceBranch);
				double branch_time = active_tip->time + rel * dt;
				vec3d real_natc_b = local_pos + (nat_dir * (factor * rel * (grow_len / possible_grow_length)));
				double final_time = branch_time;
				angio_element->_angio_mat->branch_policy->AddBranchTipEFD(angio_element, real_natc_b, seg, branch_time, next->initial_fragment_id, buffer_index, mesh);
				seg->LengthSinceBranch = seg->LengthSinceBranch - l2b;
			}
		}

		// cleanup and add recent segments and final tips.
		angio_element->recent_segments.push_back(seg);
		angio_element->reference_frame_segment_length += seg->RefLength(mesh);
		angio_element->next_tips.at(angio_element).push_back(next);
		angio_element->final_active_tips.push_back(next);
		assert(next->angio_element);
	}
	// if the vessel is going to grow into another element
	else if (proj_success)
	{
		// first grow the segment onto the face
		#ifndef NDEBUG
		#pragma omp critical
			std::cout << "case 1: face encountered" << endl;
		#endif
		//the segment only grows for a portion of dt. This portion is the amount needed to hit the face.
		vec3d next_natc = local_pos + (nat_dir * factor);
		Tip* next = new Tip(active_tip, mesh);		
		double tip_time_start = active_tip->time + (possible_grow_length / grow_len) * dt;
		assert(tip_time_start <= end_time);
		next->time = tip_time_start;
		Segment * seg = new Segment();
		next->SetLocalPosition(next_natc, mesh);
		next->TipCell->UpdateSpecies(mesh);
		next->TipCell->time = next->time;
		next->parent = seg;
		seg->back = active_tip;
		seg->front = next;
		seg->is_branch_base = seg->back->is_branch;
		if (active_tip->parent)
		{
			seg->parent = active_tip->parent;
		}
		if (angio_element->_angio_mat->branch_policy && angio_element->vessel_weight < angio_element->_angio_mat->thresh_vess_weight) {
			if (seg->is_branch_base) { seg->LengthSinceBranch = 0.0; seg->time_emerged = seg->back->time;}
			else {
				seg->LengthSinceBranch = seg->parent->LengthSinceBranch + seg->RefLength(mesh);
			}
			double l2b = dynamic_cast<DelayBranchInfo*>(angio_element->branch_info)->length_to_branch;
			if (seg->LengthSinceBranch >= l2b) {
				double rel = (l2b - seg->parent->LengthSinceBranch) / (seg->LengthSinceBranch - seg->parent->LengthSinceBranch);
				double branch_time = active_tip->time + rel * dt;
				vec3d real_natc_b = local_pos + (nat_dir * (rel * (possible_grow_length / grow_len)));
				double final_time = branch_time;
				angio_element->_angio_mat->branch_policy->AddBranchTipEFD(angio_element, real_natc_b, seg, branch_time, next->initial_fragment_id, buffer_index, mesh);
				seg->LengthSinceBranch = seg->LengthSinceBranch - l2b;
			}
		}

		// prep for the intermediary segment/tip
		vec3d pos = next->GetPosition(mesh);
		local_pos = next->GetLocalPosition();

		// find the next location we can go to
		std::vector<AngioElement*> possible_locations;
		std::vector<vec3d> possible_local_coordinates;
		for (int i = 0; i < angio_element->adjacency_list.size(); i++)
		{
			AngioElement* ang_elem = angio_element->adjacency_list[i];
			if (ang_elem)
			{
				FESolidDomain* sd = dynamic_cast<FESolidDomain*>(ang_elem->_elem->GetMeshPartition());
				if (sd)
				{
					double r[3];
					bool natc_proj = m_pangio->ProjectToElement(*ang_elem->_elem, pos, mesh, r);
					if (natc_proj && m_pangio->IsInBounds(ang_elem->_elem, r, bounds_tolerance))
					{
						possible_locations.push_back(angio_element->adjacency_list[i]);
						possible_local_coordinates.push_back(FEAngio::clamp_natc(angio_element->adjacency_list[i]->_elem->Type(), vec3d(r[0], r[1], r[2])));
					}
				}
			}
		}
		// if there is at least one location that we can grow to
		if (possible_locations.size())
		{
			//choose the next location based on which location would require the least change in direction
			int index = SelectNextTip(possible_locations, possible_local_coordinates, next, dt, buffer_index, mesh, min_scale_factor);
			if (index != -1)
			{
				#ifndef NDEBUG
					#pragma omp critical
					std::cout << "case 1a: adjacent element found" << endl;
				#endif
				vec3d real_natc_a = possible_local_coordinates[index];
				next->angio_element = possible_locations[index];
				next->TipCell->angio_element = possible_locations[index];
				next->SetLocalPosition(real_natc_a, mesh);
				next->TipCell->UpdateSpecies(mesh);
				next->face = next->angio_element;
				local_pos = next->GetLocalPosition();
				assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G4 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);
				assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G1 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);

				// cleanup and add recent segments and final tips.
				next->angio_element->recent_segments.push_back(seg);
				//SL: This line crashes occasionally...
				next->angio_element->reference_frame_segment_length += seg->RefLength(mesh);
				next->angio_element->active_tips[next_buffer_index].at(next->angio_element).push_back(next);
				}
			}
		// If the tip is near a geometric boundary:
		else
		{
			#ifndef NDEBUG
				#pragma omp critical
				std::cout << "case 1b: no adjacent element" << endl;
			#endif
			FESolidElement* se = next->angio_element->_elem;
			//! Todo: In future need to adapt this to work for any element type
			vec3d upper_bounds = vec3d(1, 1, 1);
			vec3d lower_bounds = vec3d(-1, -1, -1);
			std::vector<double> possible_values = CornerPossibleValues(local_pos, nat_dir);

			int minElemIndex = std::min_element(possible_values.begin(), possible_values.end()) - possible_values.begin();
			double d[3]; d[0] = 0; d[1] = 0; d[2] = 0;
			if (possible_values.at(minElemIndex) > 0.0) { d[minElemIndex] = 1.0; }
			else { d[minElemIndex] = -1.0; }
			vec3d face_norm = vec3d(d[0], d[1], d[2]);
			// If bounce, then bounce off face. Else grow along the face.
			vec3d plane_projection;
			if (m_pangio->bounce) {
				plane_projection = nat_dir - (face_norm * ((nat_dir * face_norm) / face_norm.norm2()) * 2.0);
			} 
			else {
				plane_projection = nat_dir - (face_norm * ((nat_dir * face_norm) / face_norm.norm2()) * 1.1);
			}
			plane_projection.unit();
			next->direction = plane_projection;
			next->use_direction = true;
			next->angio_element->active_tips[next_buffer_index].at(next->angio_element).push_back(next);
			next->angio_element->recent_segments.push_back(seg);
			next->angio_element->reference_frame_segment_length += seg->RefLength(mesh);
			}
	}
	// generally not hit unless there is large matrix deformation

	else{
		#ifndef NDEBUG
			#pragma omp critical
			std::cout << "case 2: growth into next element, projected failure" << endl;
		#endif
		angio_element->final_active_tips.push_back(active_tip);
		angio_element->next_tips.at(angio_element).push_back(active_tip);
	}
}

void FEAngioMaterial::ProtoGrowthInElement(double end_time, Tip * active_tip, int source_index, int buffer_index, double min_scale_factor, double bounds_tolerance)
{
	////////// Initialization //////////

	// get the current angio element and make sure it's active
	auto angio_element = active_tip->angio_element;

	assert(active_tip);
	// get the mesh, set the min segment length, and set the next buffer (change to write buffer from read).
	auto mesh = m_pangio->GetMesh();
	int next_buffer_index = (buffer_index + 1) % 2;
	// calculate the change in time
	double dt = end_time - active_tip->time;
	if (dt < 0) { return; }
	// determine the length the segment grows.
	double grow_vel = active_tip->GetProtoGrowthLength();
	//double grow_vel = angio_element->_angio_mat->GetInitialVelocity(angio_element);
	if (grow_vel < 0)
	{
		return;
	}
	double grow_len = grow_vel*dt;

	////////// Determine the growth length and new position //////////

	// get the local position of the tip
	double Gr[FESolidElement::MAX_NODES];
	double Gs[FESolidElement::MAX_NODES];
	double Gt[FESolidElement::MAX_NODES];
	vec3d local_pos = active_tip->GetLocalPosition();
	// get the shape function derivative values for the element type
	angio_element->_elem->shape_deriv(Gr, Gs, Gt, local_pos.x, local_pos.y, local_pos.z);
	vec3d er, es, et; // Basis vectors of the natural coordinates

	// for each node determine the contribution to the conversion between the natural and global coordinate systems.
	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		er += mesh->Node(angio_element->_elem->m_node[j]).m_rt* Gr[j];
	}
	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		es += mesh->Node(angio_element->_elem->m_node[j]).m_rt* Gs[j];
	}
	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		et += mesh->Node(angio_element->_elem->m_node[j]).m_rt* Gt[j];
	}
	
	// convert natural coords to global
	mat3d natc_to_global(er, es, et);
	mat3d global_to_natc = natc_to_global.inverse();
	bool continue_growth = true;

	// get the growth direction
	double alpha = proto_cm_manager->ApplyModifiers(dt, active_tip->angio_element, active_tip->GetLocalPosition(), mesh);
	// determine contributions from previous segments
	vec3d proto_psc_dir = proto_psc_manager->ApplyModifiers(vec3d(1, 0, 0), active_tip->angio_element, active_tip->GetLocalPosition(), active_tip->GetDirection(mesh), mesh);
	// determine contributions from local stimuli
	vec3d proto_pdd_dir = proto_pdd_manager->ApplyModifiers(vec3d(1, 0, 0), active_tip->angio_element, active_tip->GetLocalPosition(), active_tip->initial_fragment_id, buffer_index, continue_growth, proto_psc_dir, alpha, mesh, m_pangio);
	// get the new direction in global coordinates
	vec3d global_dir = mix_method->ApplyMix(proto_psc_dir, proto_pdd_dir, alpha);
	global_dir.unit();

	vec3d nat_dir = global_to_natc * global_dir;

	////////// Calculate variables to determine which element to grow into //////////

	double factor;
	// Determine if the ray cast run along the face of an element.
	bool proj_success = m_pangio->ScaleFactorToProjectToNaturalCoordinates(active_tip->angio_element->_elem, nat_dir, local_pos, factor, min_scale_factor); // TODO: HERE
	// Calculate the position on the face that the tip will grow to.
	vec3d possible_natc = local_pos + (nat_dir*factor);
	// determine the hypothetical length to grow
	double possible_grow_length = m_pangio->InElementLength(active_tip->angio_element->_elem, local_pos, possible_natc);
	
	////////// Determine which element to grow into //////////
	
	// if the vessel is going to stay in this element (or hit the face) and a face to grow to was found:
	if ((possible_grow_length >= grow_len) && proj_success)
	{
		#ifndef NDEBUG
			#pragma omp critical
			std::cout << "case 0: growth within element" << endl;
		#endif
		//determine the new natural coordinate position in this element. Not sure that this calculation needs to be this convoluted.
		vec3d real_natc = local_pos + (nat_dir * (factor * (grow_len / possible_grow_length)));

		//create the new tip
		Tip * next = new Tip(active_tip, mesh);
		next->time = end_time;
		Segment* seg = new Segment();
		next->SetLocalPosition(real_natc, mesh);
		next->TipCell->time = next->time;
		next->TipCell->ProtoUpdateSpecies(mesh);
		next->parent = seg;
		next->SetProtoGrowthLength(active_tip);
		assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G4 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);
		assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G1 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);
		
		seg->back = active_tip;
		seg->front = next;
		if (active_tip->parent)
		{
			seg->parent = active_tip->parent;
		}

		// cleanup and add recent segments and final tips.
		angio_element->recent_segments.push_back(seg);
		angio_element->reference_frame_segment_length += seg->RefLength(mesh);
		angio_element->next_tips.at(angio_element).push_back(next);
		angio_element->final_active_tips.push_back(next);
		assert(next->angio_element);
	}
	// if the vessel is going to grow into another element
	else if (proj_success)
	{
		// first grow the segment onto the face
		#ifndef NDEBUG
			#pragma omp critical
			std::cout << "case 1: face encountered" << endl;
		#endif
		//the segment only grows for a portion of dt. This portion is the amount needed to hit the face.
		vec3d next_natc = local_pos + (nat_dir * factor);
		Tip* next = new Tip(active_tip, mesh);
		double tip_time_start = active_tip->time + (possible_grow_length / grow_len) * dt;
		assert(tip_time_start < end_time);
		next->time = tip_time_start;
		Segment* seg = new Segment();
		next->SetLocalPosition(next_natc, mesh);
		next->TipCell->ProtoUpdateSpecies(mesh);
		next->TipCell->time = next->time;
		next->parent = seg;
		next->SetProtoGrowthLength(active_tip);
		seg->back = active_tip;
		seg->front = next;
		if (active_tip->parent)
		{
			seg->parent = active_tip->parent;
		}

		// prep for the intermediary segment/tip
		vec3d pos = next->GetPosition(mesh);
		local_pos = next->GetLocalPosition();
		
		// find the next location we can go to
		std::vector<AngioElement*> possible_locations;
		std::vector<vec3d> possible_local_coordinates;
		int current_mat = this->GetID_ang();
		for (int i = 0; i < angio_element->adjacency_list.size(); i++)
		{
			AngioElement * ang_elem = angio_element->adjacency_list[i];
			// if it is an angio material and it is from the same material as the current one
			int adj_mat = ang_elem->_angio_mat->GetID_ang();
			bool proto_cross = common_properties->fseeder->proto_mat_cross;
			bool same_mat = (adj_mat == current_mat);
			if (ang_elem && (proto_cross || same_mat))
			{
				FESolidDomain * sd = dynamic_cast<FESolidDomain*>(ang_elem->_elem->GetMeshPartition());
				if (sd)
				{
					double r[3];//new natural coordinates
					bool natc_proj = m_pangio->ProjectToElement(*ang_elem->_elem, pos, mesh, r);
					if (natc_proj && m_pangio->IsInBounds(ang_elem->_elem, r, bounds_tolerance))
					{
						possible_locations.push_back(angio_element->adjacency_list[i]);
						possible_local_coordinates.push_back(FEAngio::clamp_natc(angio_element->adjacency_list[i]->_elem->Type(), vec3d(r[0], r[1], r[2])));
					}
				}
			}
		}
		// if there is at least one location that we can grow to
		if (possible_locations.size())
		{
			//choose the next location based on which location would require the least change in direction
			int index = SelectNextTip(possible_locations, possible_local_coordinates, next, dt, buffer_index , mesh, min_scale_factor);
			if (index != -1)
			{
				#ifndef NDEBUG
					#pragma omp critical
					std::cout << "case 1a: adjacent element found" << endl;
				#endif
				vec3d real_natc_a = possible_local_coordinates[index];
				next->angio_element = possible_locations[index];
				next->TipCell->angio_element = possible_locations[index];
				next->SetLocalPosition(real_natc_a, mesh);
				next->TipCell->ProtoUpdateSpecies(mesh);
				next->face = next->angio_element;
				local_pos = next->GetLocalPosition();
				assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G4 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);
				assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G1 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);

				// cleanup and add recent segments and final tips.
				next->angio_element->recent_segments.push_back(seg);
				next->angio_element->reference_frame_segment_length += seg->RefLength(mesh);
				next->angio_element->active_tips[next_buffer_index].at(next->angio_element).push_back(next);
			}
		}
		// if there are no adjacent elements to grow into then 
		// currently this just keeps the tip active.
		else
		{
			#ifndef NDEBUG
				#pragma omp critical
				std::cout << "case 1b: no adjacent element" << endl;
			#endif
			//! Find which face we are encountering
			//! Todo: In future need to adapt this to work for any element type
			FESolidElement* se = active_tip->angio_element->_elem;
			vec3d upper_bounds = vec3d(1, 1, 1);
			vec3d lower_bounds = vec3d(-1, -1, -1);
			std::vector<double> possible_values = CornerPossibleValues(local_pos, nat_dir);

			int minElemIndex = std::min_element(possible_values.begin(), possible_values.end()) - possible_values.begin();
			double d[3]; d[0] = 0; d[1] = 0; d[2] = 0;
			if (possible_values.at(minElemIndex) > 0.0) { d[minElemIndex] = 1.0; }
			else { d[minElemIndex] = -1.0; }
			vec3d face_norm = vec3d(d[0], d[1], d[2]);
			// If bounce, then bounce off face. Else grow along the face.
			vec3d plane_projection;
			if (m_pangio->bounce) {
				plane_projection = nat_dir - (face_norm * ((nat_dir * face_norm) / face_norm.norm2()) * 2.0);
			}
			else {
				plane_projection = nat_dir - (face_norm * ((nat_dir * face_norm) / face_norm.norm2()) * 1.0);
			}
			plane_projection.unit();
			next->direction = plane_projection;
			next->use_direction = true;
			next->angio_element->active_tips[next_buffer_index].at(next->angio_element).push_back(next);
			next->angio_element->recent_segments.push_back(seg);
			next->angio_element->reference_frame_segment_length += seg->RefLength(mesh);
		}
	}
	// generally not hit unless there is large matrix deformation
	else {
		#ifndef NDEBUG
			#pragma omp critical
			std::cout << "case 2: growth into next element, projected failure" << endl;
		#endif
			
		angio_element->final_active_tips.push_back(active_tip);
		angio_element->next_tips.at(angio_element).push_back(active_tip);
	}
}

//needs to be selected on some criteria 
//SL: selects direction based on closest and least change in angle from current trajectory.

int FEAngioMaterial::SelectNextTip(std::vector<AngioElement*> & possible_locations, std::vector<vec3d> & possible_local_coordinates, Tip* tip, double dt, int buffer, FEMesh* mesh, double min_scale_factor)
{
#ifndef NDEBUG
#pragma omp critical
	std::cout << "Selecting new tip" << endl;
#endif

	assert(possible_locations.size() == possible_local_coordinates.size());
	if (possible_locations.size() == 1)
		return 0;
	auto dir = tip->GetDirection(mesh); dir.unit();
	std::vector<double> angles;
	std::vector<double> dists;
	bool continue_growth = true;
	for(int i=0; i < possible_locations.size();i++)
	{
		if(!continue_growth)
		{
			return i;
		}
		vec3d nat_dir = tip->GetLocalPosition() - possible_local_coordinates[i]; nat_dir.unit();
		//vec3d possible_dir = possible_locations[i]->_angio_mat->pdd_manager->ApplyModifiers({ 1,0,0 }, possible_locations[i], possible_local_coordinates[i], tip->initial_fragment_id, buffer, continue_growth, psc_dir, alpha, mesh, m_pangio);
		double Gr[FESolidElement::MAX_NODES]; double Gs[FESolidElement::MAX_NODES]; double Gt[FESolidElement::MAX_NODES];
		vec3d local_pos = possible_local_coordinates[i];
		possible_locations[i]->_elem->shape_deriv(Gr, Gs, Gt, local_pos.x, local_pos.y, local_pos.z);
		vec3d er, es, et; // Basis vectors of the natural coordinates
		for (int j = 0; j < possible_locations[i]->_elem->Nodes(); j++)
		{
			er += mesh->Node(possible_locations[i]->_elem->m_node[j]).m_rt* Gr[j];
		}
		for (int j = 0; j < possible_locations[i]->_elem->Nodes(); j++)
		{
			es += mesh->Node(possible_locations[i]->_elem->m_node[j]).m_rt* Gs[j];
		}
		for (int j = 0; j < possible_locations[i]->_elem->Nodes(); j++)
		{
			et += mesh->Node(possible_locations[i]->_elem->m_node[j]).m_rt* Gt[j];
		}
		mat3d natc_to_global(er, es, et);
		nat_dir = natc_to_global * nat_dir; nat_dir.unit();
		vec3d global_pos = this->m_pangio->Position(possible_locations[i]->_elem, local_pos);
		dists.push_back((global_pos - tip->GetPosition(mesh)).norm());
		angles.push_back(dir * nat_dir);
	}
	double min = dists[0];
	double mina = 0;
	int index = 0;
	for(int i=0;i < dists.size();i++)
	{
		if(dists[i] < min && angles[i] > mina)
		{
			min = dists[i];
			mina = angles[i];
			index = i;
		}
	}
	return index;
}

void FEAngioMaterial::PostGrowthUpdate(AngioElement* angio_elem, double end_time, double final_time, double min_scale_factor, int buffer_index, FEMesh* mesh, FEAngio* feangio)
{	

	//do the buffer management
	for(auto iter = angio_elem->active_tips[buffer_index].begin(); iter != angio_elem->active_tips[buffer_index].end(); ++iter)
	{
		iter->second.clear();
	}
}

void FEAngioMaterial::ProtoPostGrowthUpdate(AngioElement* angio_elem, double end_time, double min_scale_factor, int buffer_index, FEMesh* mesh, FEAngio* feangio)
{

	//do the buffer management
	for (auto iter = angio_elem->active_tips[buffer_index].begin(); iter != angio_elem->active_tips[buffer_index].end(); ++iter)
	{
		iter->second.clear();
	}
}

void FEAngioMaterial::Cleanup(AngioElement* angio_elem, double end_time, int buffer_index)
{
	for(int i=0;i < angio_elem->recent_segments.size();i++)
	{
		angio_elem->grown_segments.push_back(angio_elem->recent_segments[i]);
	}
	angio_elem->recent_segments.clear();
	angio_elem->processed_recent_segments = 0;
}

void FEAngioMaterial::PrepBuffers(AngioElement* angio_elem, double end_time, int buffer_index)
{
	angio_elem->next_tips.swap(angio_elem->active_tips[buffer_index]);
	for(auto iter=angio_elem->next_tips.begin(); iter != angio_elem->next_tips.end(); ++iter)
	{
		iter->second.clear();
	}

	std::vector<Tip*> temp;
	temp.reserve(angio_elem->final_active_tips.size());
	for(int i=0; i < angio_elem->final_active_tips.size();i++)
	{
		//recycle tips that have not been hit temporally
		if(angio_elem->final_active_tips[i]->time > end_time)
		{
			temp.push_back(angio_elem->final_active_tips[i]);
		}
	}
	angio_elem->final_active_tips = temp;
}

bool FEAngioMaterial::SeedFragments(std::vector<AngioElement *>& angio_elements, FEMesh* mesh)
{
	return common_properties->fseeder->SeedFragments(angio_elements, mesh, this, 0);
}



//-----------------------------------------------------------------------------
FEMaterialPoint* FEAngioMaterial::CreateMaterialPointData()
{
	return new FEAngioMaterialPoint((FEElasticMaterial::CreateMaterialPointData()), common_properties->vessel_material->CreateMaterialPointData(), matrix_material->CreateMaterialPointData());
}