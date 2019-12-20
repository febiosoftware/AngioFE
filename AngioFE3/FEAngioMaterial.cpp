#include "StdAfx.h"
#include "FEAngioMaterial.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include "FEBioMech/FEElasticMixture.h"
#include "FEBioMech/FEFiberMaterialPoint.h"
#include "FECore/FEElementTraits.h"
#include "FEBioMech/FEViscoElasticMaterial.h"
#include "FEBioMech/FESPRProjection.h"
#include <iostream>
#include "angio3d.h"
#include "AngioElement.h"
#include "Segment.h"
#include "Tip.h"
#include "VariableInterpolation.h"


//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEAngioMaterial, FEElasticMaterial)
	
	ADD_PARAMETER(initial_segment_velocity, FE_PARAM_DOUBLE, "initial_segment_velocity");
	ADD_PARAMETER(vessel_radius, FE_PARAM_DOUBLE, "vessel_radius");

END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEAngioMaterial::FEAngioMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	AddProperty(&common_properties, "common_properties");
	AddProperty(&matrix_material, "matrix");
	AddProperty(&angio_stress_policy, "angio_stress_policy");

	AddProperty(&pdd_manager, "pdd_manager");
	AddProperty(&psc_manager, "psc_manager");
	AddProperty(&cm_manager, "cm_manager");
	AddProperty(&mix_method, "mix_method");
	AddProperty(&velocity_manager, "velocity_manager");
	AddProperty(&im_manager, "im_manager"); im_manager.m_brequired = false;
	AddProperty(&nodedata_interpolation_manager, "nodedata_interpolation_manager"); nodedata_interpolation_manager.m_brequired = false;
	AddProperty(&branch_policy, "branch_policy"); branch_policy.m_brequired = false;
	AddProperty(&proto_branch_policy, "proto_branch_policy"); proto_branch_policy.m_brequired = false;
	AddProperty(&tip_species_manager, "tip_species_manager"); tip_species_manager.m_brequired = false;
}

FEAngioMaterial::~FEAngioMaterial()
{
}

//-----------------------------------------------------------------------------
bool FEAngioMaterial::Init()
{
	// Create symmetry vectors


	if(!matrix_material->Init()) return false;

	if(!common_properties->vessel_material->Init()) return false;

	if (FEElasticMaterial::Init() == false) return false;
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


void FEAngioMaterial::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
	// get the material's coordinate system (if defined)
	FECoordSysMap* pmap = GetCoordinateSystemMap();
	//this allows the local coordinates to work correctly
	if (pmap)
	{
		FEElasticMaterial::SetLocalCoordinateSystem(el, n, mp);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);

		FEElasticMaterialPoint& vessel_elastic = *angioPt->vessPt->ExtractData<FEElasticMaterialPoint>();
		FEElasticMaterialPoint& matrix_elastic = *angioPt->matPt->ExtractData<FEElasticMaterialPoint>();

		// compound the local map with the global material axes
		//mat3d Qlocal = pmap->LocalElementCoord(el, n);
		//pt.m_Q = pt.m_Q * Qlocal;

		vessel_elastic.m_Q = pt.m_Q;
		matrix_elastic.m_Q = pt.m_Q;

		FEElasticMaterial* vess_elastic = common_properties->vessel_material->GetElasticMaterial();
		FEElasticMaterial* mat_elastic = matrix_material->GetElasticMaterial();

		vess_elastic->SetLocalCoordinateSystem(el, n, *angioPt->vessPt);
		mat_elastic->SetLocalCoordinateSystem(el, n, *angioPt->matPt);
	}
	
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

		vessel_elastic.m_rt = elastic_pt.m_rt;//spatial position
		vessel_elastic.m_r0 = elastic_pt.m_r0;//material position
		vessel_elastic.m_F = elastic_pt.m_F;//deformation gradient
		vessel_elastic.m_J = elastic_pt.m_J;//determinate
		
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
		sed = angioPt->vessel_weight*common_properties->vessel_material->GetElasticMaterial()->StrainEnergyDensity(*angioPt->vessPt) + angioPt->matrix_weight*matrix_material->GetElasticMaterial()->StrainEnergyDensity(*angioPt->matPt);
	}
	return sed;
}

void FEAngioMaterial::SetSeeds(AngioElement* angio_elem)
{
	static int offset = 0;
	static int seed = int (m_pangio->m_fem->GetGlobalConstant("seed"));
	std::seed_seq sq{ offset + seed };
	angio_elem->_rengine.seed(sq);
	offset++;
}

double FEAngioMaterial::GetSegmentVelocity(AngioElement * angio_element, vec3d local_pos, FEMesh* mesh)
{
	double const vel_init = 1;
	return velocity_manager->ApplyModifiers(vel_init, local_pos, angio_element, mesh);
}

double FEAngioMaterial::GetInitialVelocity() const
{
	return initial_segment_velocity;
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
		double vel = angio_elem->_angio_mat->GetSegmentVelocity(angio_elem, pos, mesh);
		max_grow_velocity = std::max(max_grow_velocity, vel);
	}

	double min_side_length_so_far = std::numeric_limits<double>::max();

	for (int i = 0; i < angio_elem->_elem->GaussPoints(); i++)
	{
		double Gr[FEElement::MAX_NODES];
		double Gs[FEElement::MAX_NODES];
		double Gt[FEElement::MAX_NODES];
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

void FEAngioMaterial::GrowSegments(AngioElement * angio_elem, double end_time, int buffer_index, double min_scale_factor, double bounds_tolerance, double min_angle)
{
	assert(angio_elem);
	for(int i=0; i < angio_elem->adjacency_list.size();i++)
	{
		if(angio_elem->adjacency_list.at(i))
		{
			std::vector<Tip*> & tips = angio_elem->adjacency_list.at(i)->active_tips[buffer_index].at(angio_elem);
			for(int j=0; j < tips.size();j++)
			{
				assert(tips.at(j)->angio_element == angio_elem);
				GrowthInElement(end_time, tips.at(j), i, buffer_index, min_scale_factor,bounds_tolerance, min_angle);
			}
			
		}
		
	}

	std::vector<Tip*> & tips = angio_elem->active_tips[buffer_index].at(angio_elem);
	for(int j=0; j < tips.size();j++)
	{
		assert(tips[j]->angio_element == angio_elem);
		GrowthInElement(end_time, tips[j], -1, buffer_index, min_scale_factor,bounds_tolerance, min_angle);
	}
}

void FEAngioMaterial::ProtoGrowSegments(AngioElement * angio_elem, double end_time, int buffer_index, double min_scale_factor, double bounds_tolerance, double min_angle)
{
	assert(angio_elem);
	for (int i = 0; i < angio_elem->adjacency_list.size(); i++)
	{
		if (angio_elem->adjacency_list.at(i))
		{
			std::vector<Tip*> & tips = angio_elem->adjacency_list.at(i)->active_tips[buffer_index].at(angio_elem);
			for (int j = 0; j < tips.size(); j++)
			{
				assert(tips.at(j)->angio_element == angio_elem);
				ProtoGrowthInElement(end_time, tips.at(j), i, buffer_index, min_scale_factor, bounds_tolerance, min_angle);
			}
		}
	}

	std::vector<Tip*> & tips = angio_elem->active_tips[buffer_index].at(angio_elem);
	for (int j = 0; j < tips.size(); j++)
	{
		assert(tips[j]->angio_element == angio_elem);
		ProtoGrowthInElement(end_time, tips[j], -1, buffer_index, min_scale_factor, bounds_tolerance,min_angle);
	}
}

//Grow a given active tip
void FEAngioMaterial::GrowthInElement(double end_time, Tip * active_tip, int source_index, int buffer_index, double min_scale_factor, double bounds_tolerance, double min_angle)
{
	// get the current angio element and make sure it's active
	auto angio_element = active_tip->angio_element;
	
	assert(active_tip);
	// get the mesh, set the min segment length, and determine the next buffer index.
	auto mesh = m_pangio->GetMesh();
	const double min_segm_len = 0.1;
	int next_buffer_index = (buffer_index + 1) % 2;
	// calculate the change in time
	double dt = end_time - active_tip->time;
	assert(dt > 0.0);
	// determine the length the segment grows.
	double grow_vel = this->GetSegmentVelocity(angio_element, active_tip->GetLocalPosition(), mesh);
	//If the velocity is negative wrt the tip direction then end the tip/
	if(grow_vel < 0)
	{
		return;
	}
	assert(grow_vel > 0.0);
	double grow_len = grow_vel*dt;

	double Gr[FEElement::MAX_NODES];
	double Gs[FEElement::MAX_NODES];
	double Gt[FEElement::MAX_NODES];
	// get the local position of the tip
	vec3d local_pos = active_tip->GetLocalPosition();
	// get the shape function derivative values for the element type
	angio_element->_elem->shape_deriv(Gr, Gs, Gt, local_pos.x, local_pos.y, local_pos.z);
	vec3d er, es, et; // Basis vectors of the natural coordinates
	// for each node determine something about the position
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
	// determine the inverse operation to go from global to natural
	mat3d global_to_natc = natc_to_global.inverse();
	bool continue_growth = true;
	// determine the alpha (mix value)
	double alpha = cm_manager->ApplyModifiers(dt, active_tip->angio_element, active_tip->GetLocalPosition(), mesh);
	// determine contributions from previous segments
	vec3d psc_dir = psc_manager->ApplyModifiers(vec3d(1, 0, 0), active_tip->angio_element, active_tip->GetLocalPosition(), active_tip->GetDirection(mesh), mesh);
	// determine contributions from local stimuli
	vec3d pdd_dir = pdd_manager->ApplyModifiers(vec3d(1, 0, 0), active_tip->angio_element, active_tip->GetLocalPosition() , active_tip->initial_fragment_id , buffer_index , continue_growth, psc_dir, alpha, mesh, m_pangio);
		alpha = std::max(0.0,std::min(1.0, alpha));
	if(!continue_growth)
	{
		return;
	}
	// get the new direction in global coordinates
	vec3d global_dir = mix_method->ApplyMix(psc_dir, pdd_dir, alpha);
	global_dir.unit();
	// get new position in global coordinates
	vec3d global_pos;
	// transform global to local
	double H[FEElement::MAX_NODES];
	angio_element->_elem->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		global_pos += mesh->Node(angio_element->_elem->m_node[j]).m_rt* H[j];
	}
	vec3d nat_dir = global_to_natc * global_dir;
	double factor;         
	// Determine if the ray cast run along the face of an element. If you see vessels bouncing off unexposed faces this is likely the issue.
	//This also updates "factor" which scales the length allowed in the element.
	bool proj_success = m_pangio->ScaleFactorToProjectToNaturalCoordinates(active_tip->angio_element->_elem, nat_dir, local_pos, factor, min_scale_factor); // TODO: HERE
	// Calculate new possible natural coordinate location.
	vec3d possible_natc = local_pos + (nat_dir*factor);
	// determine the hypothetical length to grow
	double possible_grow_length = m_pangio->InElementLength(active_tip->angio_element->_elem, local_pos, possible_natc);
	// if the vessel is going to grow into another element and and success is projected:
	if ((possible_grow_length >= grow_len) && proj_success)
	{
		//determine the length to grow in this element
		vec3d real_natc = local_pos + (nat_dir * (factor * (grow_len / possible_grow_length)));
		//add the segment and add the new tip to the current angio element
		Tip * next = new Tip(active_tip, mesh);
		Segment * seg = new Segment();
		next->time = end_time;
		// the new tip is assigned to the current element.
		next->angio_element = angio_element;
		// assign the local position to the new tip.
		next->SetLocalPosition(real_natc, mesh);

		// Update SBM boundary condition if necessary.
		//
		//next->TipSBM->UpdatePos(natc_to_global * real_natc);
		//

		// assign the new tip as the parent of the segment.
		next->parent = seg;
		// assign the current element as location where growth began.
		next->face = angio_element;
		// assign the velocity that the tip grew at
		next->growth_velocity = grow_vel;
		// copy the parent fragment id from the original tip.
		next->initial_fragment_id = active_tip->initial_fragment_id;
		// ensure that if the element is a tet that it does not exceed the element boundaries.
		assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G4 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);
		assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G1 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);
		// add the new tip to the element.
		angio_element->next_tips.at(angio_element).push_back(next);

		//move next to use the rest of the remaining dt

		seg->back = active_tip;
		seg->front = next;
		if (active_tip->parent)
		{
			seg->parent = active_tip->parent;
		}

		angio_element->refernce_frame_segment_length += seg->RefLength(mesh);
		angio_element->recent_segments.push_back(seg);
		angio_element->final_active_tips.push_back(next);
		assert(next->angio_element);
	}
	else if(proj_success)
	{
		//the segment only grows for a portion of dt
		double tip_time_start = active_tip->time + (possible_grow_length / grow_len) * dt;
		assert(tip_time_start < end_time);
		//grow len should always be nonzero so this division should be okay
		//this segment is growing into a new element do the work to do this

		//next is just the end of the segment not the begining of the next segment
		Tip * next = new Tip();
		Segment * seg = new Segment();
		next->time = tip_time_start;

		next->angio_element = angio_element;

		next->Species = active_tip->Species;
		active_tip->Species.clear();
		
		//still need to update the local position of the tip
		next->SetLocalPosition(local_pos + (nat_dir * possible_grow_length), mesh);

		next->parent = seg;
		next->face = angio_element;
		next->initial_fragment_id = active_tip->initial_fragment_id;
		next->growth_velocity = grow_vel;

		seg->back = active_tip;
		seg->front = next;
		if (active_tip->parent)
		{
			seg->parent = active_tip->parent;
		}

		vec3d pos = next->GetPosition(mesh);
		angio_element->recent_segments.push_back(seg);
		angio_element->refernce_frame_segment_length += seg->RefLength(mesh);

		std::vector<AngioElement*> possible_locations;
		std::vector<vec3d> possible_local_coordinates;
		// check to see if the element is a hex how many exposed faces there are. 0 -> 26, 1 -> 17, 2 -> 11, 3 -> 7
		/*if (angio_element->adjacency_list.size() != 26) {
			std::cout << "exposed, size is " << angio_element->adjacency_list.size() << endl;
		}*/
		for(int i=0; i < angio_element->adjacency_list.size();i++) // TODO: HERE
		{
			AngioElement * ang_elem = angio_element->adjacency_list[i];
			if(ang_elem)
			{
				FESolidDomain * sd = dynamic_cast<FESolidDomain*>(ang_elem->_elem->GetDomain());
				if(sd)
				{
					double r[3];//new natural coordinates
					bool natc_proj = m_pangio->ProjectToElement(*ang_elem->_elem,pos,mesh, r);
					if(natc_proj &&  m_pangio->IsInBounds(ang_elem->_elem, r, bounds_tolerance))
					{
						possible_locations.push_back(angio_element->adjacency_list[i]);
						possible_local_coordinates.push_back(FEAngio::clamp_natc(angio_element->adjacency_list[i]->_elem->Type(),vec3d(r[0], r[1], r[2])));
					}
				}
			}
		}
		// if there is at least one location that we can grow to
		if(possible_locations.size())
		{
			//need some way to choose the correct element to continue the tip in
			int index = SelectNextTip(possible_locations, possible_local_coordinates, next, dt , buffer_index, mesh,min_scale_factor, min_angle);
			if(index != -1)
			{
				Tip * adj = new Tip(next, mesh);
				adj->angio_element = possible_locations[index];
				adj->face = angio_element;
				adj->SetLocalPosition(possible_local_coordinates[index], mesh);
				adj->use_direction = true;
				adj->growth_velocity = grow_vel;
				angio_element->active_tips[next_buffer_index].at(possible_locations[index]).push_back(adj);
			}
		}
		// If the tip is near a geometric boundary:
		else
		{
			// in this element assign this tip for evaluation on the next step. We will assign it at the current element by pushing back the active tip.
			angio_element->next_tips.at(angio_element).push_back(active_tip);
		}
	}
	// this is not being hit rn because proj_success always returns true.
	else
	{
	}
}

void FEAngioMaterial::ProtoGrowthInElement(double end_time, Tip * active_tip, int source_index, int buffer_index, double min_scale_factor, double bounds_tolerance, double min_angle)
{
	auto angio_element = active_tip->angio_element;

	assert(active_tip);
	auto mesh = m_pangio->GetMesh();
	const double min_segm_len = 0.1;
	int next_buffer_index = (buffer_index + 1) % 2;
	double dt = end_time - active_tip->time;
	double grow_vel = angio_element->_angio_mat->GetInitialVelocity();
	double grow_len = grow_vel *dt;

	double Gr[FEElement::MAX_NODES];
	double Gs[FEElement::MAX_NODES];
	double Gt[FEElement::MAX_NODES];
	vec3d local_pos = active_tip->GetLocalPosition();
	angio_element->_elem->shape_deriv(Gr, Gs, Gt, local_pos.x, local_pos.y, local_pos.z);
	vec3d er, es, et; // Basis vectors of the natural coordinates
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
	mat3d natc_to_global(er, es, et);
	mat3d global_to_natc = natc_to_global.inverse();
	vec3d psc_dir = psc_manager->ApplyModifiers(vec3d(1, 0, 0), active_tip->angio_element, active_tip->GetLocalPosition(), active_tip->GetDirection(mesh), mesh);
	vec3d global_dir = psc_dir;
	global_dir.unit();

	vec3d global_pos;
	double H[FEElement::MAX_NODES];
	angio_element->_elem->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		global_pos += mesh->Node(angio_element->_elem->m_node[j]).m_rt* H[j];
	}

	// Just do the transformations
	vec3d nat_dir = global_to_natc * global_dir;
	double factor;
	bool proj_success = m_pangio->ScaleFactorToProjectToNaturalCoordinates(active_tip->angio_element->_elem, nat_dir, local_pos, factor, min_scale_factor); // TODO: HERE
	vec3d possible_natc = local_pos + (nat_dir*factor);
	double possible_grow_length = m_pangio->InElementLength(active_tip->angio_element->_elem, local_pos, possible_natc);
	if ((possible_grow_length >= grow_len) && proj_success)
	{
		//this is still in the same element
		vec3d real_natc = local_pos + (nat_dir * (factor * (grow_len / possible_grow_length)));
		//add the segment and add the new tip to the current angio element
		Tip * next = new Tip(active_tip, mesh);
		Segment * seg = new Segment();
		next->time = end_time;
		next->angio_element = angio_element;
		next->SetLocalPosition(real_natc, mesh);
		next->parent = seg;
		next->face = angio_element;
		//growth velocity must be zero to work with grown segments stress policy
		//next->growth_velocity = grow_vel;
		next->growth_velocity = 0;
		next->initial_fragment_id = active_tip->initial_fragment_id;
		assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G4 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);
		assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G1 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);
		angio_element->next_tips.at(angio_element).push_back(next);

		//move next to use the rest of the remaining dt

		seg->back = active_tip;
		seg->front = next;
		if (active_tip->parent)
		{
			seg->parent = active_tip->parent;
		}


		angio_element->recent_segments.push_back(seg);
		angio_element->refernce_frame_segment_length += seg->RefLength(mesh);
		angio_element->final_active_tips.push_back(next);
		assert(next->angio_element);
	}
	else if (proj_success)
	{
		//the segment only grows for a portion of dt
		double tip_time_start = active_tip->time + (possible_grow_length / grow_len) * dt;
		assert(tip_time_start < end_time);
		//grow len should always be nonzero so this division should be okay
		//this segment is growing into a new element do the work to do this

		//next is just the end of the segment not the begining of the next segment
		Tip * next = new Tip();
		Segment * seg = new Segment();
		next->time = tip_time_start;

		next->angio_element = angio_element;
		next->Species = active_tip->Species;
		active_tip->Species.clear();

		//still need to update the local position of the tip
		next->SetLocalPosition(local_pos + (nat_dir * possible_grow_length), mesh);

		next->parent = seg;
		next->face = angio_element;
		next->initial_fragment_id = active_tip->initial_fragment_id;
		//growth velocity must be zero to work with grown segments stress policy
		next->growth_velocity = 0;

		seg->back = active_tip;
		seg->front = next;
		if (active_tip->parent)
		{
			seg->parent = active_tip->parent;
		}

		vec3d pos = next->GetPosition(mesh);
		angio_element->refernce_frame_segment_length += seg->RefLength(mesh);
		angio_element->recent_segments.push_back(seg);

		std::vector<AngioElement*> possible_locations;
		std::vector<vec3d> possible_local_coordinates;
		for (int i = 0; i < angio_element->adjacency_list.size(); i++) // TODO: HERE
		{
			AngioElement * ang_elem = angio_element->adjacency_list[i];
			if (ang_elem)
			{
				FESolidDomain * sd = dynamic_cast<FESolidDomain*>(ang_elem->_elem->GetDomain());
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
		if (possible_locations.size())
		{
			//need some way to choose the correct element to continue the tip in
			int index = SelectNextTip(possible_locations, possible_local_coordinates, next, dt, buffer_index , mesh, min_scale_factor, min_angle);
			if (index != -1)
			{
				Tip * adj = new Tip(next, mesh);
				adj->angio_element = possible_locations[index];
				adj->face = angio_element;
				adj->SetLocalPosition(possible_local_coordinates[index], mesh);
				adj->use_direction = true;
				//growth velocity must be zero to work with grown segments stress policy
				adj->growth_velocity = 0;
				angio_element->active_tips[next_buffer_index].at(possible_locations[index]).push_back(adj);
			}
		}
		else
		{
			// in this element assign this tip for evaluation on the next step. We will assign it at the current element by pushing back the active tip.
			angio_element->next_tips.at(angio_element).push_back(active_tip);
		}
	}
}

//needs to be selected on some criteria 
//possibilites include: longest possible growth length,
//most similar in direction to the tip's direction and above a 

int FEAngioMaterial::SelectNextTip(std::vector<AngioElement*> & possible_locations, std::vector<vec3d> & possible_local_coordinates, Tip* tip, double dt, int buffer, FEMesh* mesh, double min_scale_factor, double min_angle)
{
	assert(possible_locations.size() == possible_local_coordinates.size());
	if (possible_locations.size() == 1)
		return 0;
	auto dir = tip->GetDirection(mesh);
	std::vector<double> angles;
	bool continue_growth = true;
	for(int i=0; i < possible_locations.size();i++)
	{
		double alpha = cm_manager->ApplyModifiers(dt, possible_locations[i], possible_local_coordinates[i], mesh);
		vec3d psc_dir = psc_manager->ApplyModifiers(vec3d(1, 0, 0), possible_locations[i], possible_local_coordinates[i], tip->GetDirection(mesh), mesh);
		vec3d pdd_dir = pdd_manager->ApplyModifiers(vec3d(1, 0, 0), possible_locations[i], possible_local_coordinates[i],tip->initial_fragment_id , buffer , continue_growth , psc_dir , alpha, mesh, m_pangio);
		if(!continue_growth)
		{
			return i;
		}
		vec3d global_dir = mix_method->ApplyMix(psc_dir, pdd_dir, (alpha));
		//vec3d global_dir = mix(psc_dir, pdd_dir, (alpha));
		global_dir.unit();
		vec3d possible_dir = possible_locations[i]->_angio_mat->pdd_manager->ApplyModifiers({ 1,0,0 }, possible_locations[i], possible_local_coordinates[i], tip->initial_fragment_id, buffer, continue_growth, psc_dir, alpha, mesh, m_pangio);
		double Gr[FEElement::MAX_NODES];
		double Gs[FEElement::MAX_NODES];
		double Gt[FEElement::MAX_NODES];
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
		mat3d global_to_natc = natc_to_global.inverse();
		vec3d nat_dir = global_to_natc * possible_dir;
		double factor;
		bool proj_success = m_pangio->ScaleFactorToProjectToNaturalCoordinates(possible_locations[i]->_elem, nat_dir, local_pos, factor, min_scale_factor);
		if(proj_success && factor > min_scale_factor && possible_dir*dir >= min_angle)
		{
			angles.push_back(factor);
		}
		else
		{
			angles.push_back(-1);
		}
	}
	double min = -1;
	int index = -1;
	for(int i=0;i < angles.size();i++)
	{
		if(angles[i] > min)
		{
			min = angles[i];
			index = i;
		}
	}

	

	return index;
}

void FEAngioMaterial::PostGrowthUpdate(AngioElement* angio_elem, double end_time, double final_time, double min_scale_factor, int buffer_index, FEMesh* mesh, FEAngio* feangio)
{
	//do the branching
	if(angio_elem->_angio_mat->branch_policy)
	{
		angio_elem->_angio_mat->branch_policy->AddBranches(angio_elem, buffer_index, end_time , final_time, min_scale_factor, mesh, feangio);
	}
	

	//do the buffer management
	for(auto iter = angio_elem->active_tips[buffer_index].begin(); iter != angio_elem->active_tips[buffer_index].end(); ++iter)
	{
		iter->second.clear();
	}
}

void FEAngioMaterial::ProtoPostGrowthUpdate(AngioElement* angio_elem, double end_time, double min_scale_factor, int buffer_index, FEMesh* mesh, FEAngio* feangio)
{
	//do the branching
	if (angio_elem->_angio_mat->proto_branch_policy)
	{
		angio_elem->_angio_mat->proto_branch_policy->AddBranches(angio_elem, buffer_index, end_time, std::numeric_limits<int>::max(), min_scale_factor , mesh, feangio);
	}


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
	return new FEAngioMaterialPoint(new FEFiberMaterialPoint(FEElasticMaterial::CreateMaterialPointData()), common_properties->vessel_material->CreateMaterialPointData(), matrix_material->CreateMaterialPointData());
}