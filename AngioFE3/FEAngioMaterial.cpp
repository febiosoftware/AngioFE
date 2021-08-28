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
#include "angio3d.h"
#include "AngioElement.h"
#include "Segment.h"
#include "Tip.h"
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
	AddClassProperty(this, &tip_species_manager, "tip_species_manager", FEProperty::Optional);
	
}

FEAngioMaterial::~FEAngioMaterial()
{
}

//-----------------------------------------------------------------------------
bool FEAngioMaterial::Init()
{
	// Create symmetry vectors


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

void FEAngioMaterial::GrowSegments(AngioElement * angio_elem, double end_time, int buffer_index, double min_scale_factor, double bounds_tolerance, double min_angle)
{
	assert(angio_elem);

	// Grow tips on shared faces that are growing into/out of? the element. This should handle passing of tips between adjacent elements.
	// for each adjacent element
	for(int i=0; i < angio_elem->adjacency_list.size();i++)
	{
		// if there are adjacent elements
		if(angio_elem->adjacency_list.at(i))
		{
			// create a vector of tips on the shared face?
			std::vector<Tip*> & tips = angio_elem->adjacency_list.at(i)->active_tips[buffer_index].at(angio_elem);
			// for each tip in the vector ensure that ... then call growth in element.
			for(int j=0; j < tips.size();j++)
			{
				assert(tips.at(j)->angio_element == angio_elem);
				// call GrowthInElement providing the end time of the step, the active tip, the source index of the adjacent element, the buffer index, and other 3 things.
				GrowthInElement(end_time, tips.at(j), i, buffer_index, min_scale_factor,bounds_tolerance, min_angle);
			}
			
		}
		
	}

	// Grow tips that originate in this element and have not reached an adjacent face.
	// create a vector for the active tips in the element
	std::vector<Tip*> & tips = angio_elem->active_tips[buffer_index].at(angio_elem);
	// for each tip ensure that the tip belongs to the element being analyzed. Next call growth in the element.
	for(int j=0; j < tips.size();j++)
	{
		assert(tips[j]->angio_element == angio_elem);
		// call GrowthInElement providing the end time, the active tip, -1 to indicate the tip originated within this element, the buffer index, etc.
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
			//copy the vector of active tips in the adjacent angio element
			std::vector<Tip*> & tips = angio_elem->adjacency_list.at(i)->active_tips[buffer_index].at(angio_elem);
			// for each tip in the adjacent element grow the tips
			for (int j = 0; j < tips.size(); j++)
			{
				assert(tips.at(j)->angio_element == angio_elem);
				ProtoGrowthInElement(end_time, tips.at(j), i, buffer_index, min_scale_factor, bounds_tolerance, min_angle);
			}
		}
	}

	// copy the vector of active tips in the current angio element
	std::vector<Tip*> & tips = angio_elem->active_tips[buffer_index].at(angio_elem);\
	// for each tip in the current element grow the tips.
	for (int j = 0; j < tips.size(); j++)
	{
		assert(tips[j]->angio_element == angio_elem);
		ProtoGrowthInElement(end_time, tips[j], -1, buffer_index, min_scale_factor, bounds_tolerance,min_angle);
	}
}

//Grow a provided active tip based on the dt, originating element, etc.
void FEAngioMaterial::GrowthInElement(double end_time, Tip * active_tip, int source_index, int buffer_index, double min_scale_factor, double bounds_tolerance, double min_angle)
{
	////////// Initialization //////////

	// get the current angio element and make sure it's active
	auto angio_element = active_tip->angio_element;
	
	assert(active_tip);
	// get the mesh, set the min segment length, and set the next buffer (change to write buffer from read).
	auto mesh = m_pangio->GetMesh();
	const double min_segm_len = 5;
	
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
	// determine the inverse operation to go from global to natural
	mat3d global_to_natc = natc_to_global.inverse();
	bool continue_growth = true;

	// get the growth direction
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
	double H[FESolidElement::MAX_NODES];
	angio_element->_elem->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		global_pos += mesh->Node(angio_element->_elem->m_node[j]).m_rt* H[j];
	}
	vec3d nat_dir = global_to_natc * global_dir;
	
	////////// Calculate variables to determine which element to grow into //////////
	
	double factor;         
	// Determine if the ray cast run along the face of an element. If you see vessels bouncing off unexposed faces this is likely the issue.
	//This also updates "factor" which scales the length allowed in the element. It is the shortest of the 3 directions that an element can grow before exceeding the bounds while also staying above the min_scale_factor
	//If it fails that means that all possible factors are less than the min_scale_factor.
	bool proj_success = m_pangio->ScaleFactorToProjectToNaturalCoordinates(active_tip->angio_element->_elem, nat_dir, local_pos, factor, min_scale_factor); // TODO: HERE
	//bool proj_success = true;
																																							// Calculate the position on the face that the tip will grow to.
	vec3d possible_natc = local_pos + (nat_dir*factor);
	// determine the hypothetical length to grow
	double possible_grow_length = m_pangio->InElementLength(active_tip->angio_element->_elem, local_pos, possible_natc);
		
	////////// Determine which element to grow into //////////

	// if the vessel is going to stay in this element (or hit the face) and a face to grow to was found:
	if ((possible_grow_length >= grow_len) && proj_success)
	{
		//determine the new natural coordinate position in this element. Not sure that this calculation needs to be this convoluted.
		vec3d real_natc = local_pos + (nat_dir * (factor * (grow_len / possible_grow_length)));
		//add the segment and add the new tip to the current angio element
		Tip * next = new Tip(active_tip, mesh);
		Segment * seg = new Segment();
		next->time = end_time;
		// the new tip is assigned to the current element.
		next->angio_element = angio_element;
		// assign the new position to the new tip.
		next->SetLocalPosition(real_natc, mesh);
		// assign the new tip as the parent of the segment.
		next->parent = seg;
		// assign the current element as location where growth began.
		next->face = angio_element;
		// assign the velocity that the tip grew at
		next->growth_velocity = grow_vel;
		// copy the parent fragment id from the original tip.
		next->initial_fragment_id = active_tip->initial_fragment_id;
		// ensure that if the element is a tet that it does not exceed the element boundaries. This should already have been done in ScaleFactor calculation.
		assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G4 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);
		assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G1 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);
		// add the new tip to the element so the element manages the tip.
		angio_element->next_tips.at(angio_element).push_back(next);
		// make the active tip the back of this segment and make next the front.
		seg->back = active_tip;
		seg->front = next;
		// assign the segment's parent seg as the same as the active tip's parent seg
		if (active_tip->parent)
		{
			seg->parent = active_tip->parent;
		}

		// add the segment to the recent segments vector
		angio_element->recent_segments.push_back(seg);
		// calculate the segment length of seg
		angio_element->refernce_frame_segment_length += seg->RefLength(mesh);
		// add next to the final active tips
		angio_element->final_active_tips.push_back(next);
		assert(next->angio_element);
	}
	else if (proj_success)
	{
		//the segment only grows for a portion of dt. This portion is the amount needed to hit the face.
		//grow len should always be nonzero so this division should be okay. This will need to be accounted for in the future if we explicitly model sprouting.
		double tip_time_start = active_tip->time + (possible_grow_length / grow_len) * dt;
		assert(tip_time_start < end_time);
		double ct = m_pangio->GetFEModel()->GetTime().currentTime;
		//std::cout << "current time is " << ct << endl;
		// check that the tip isn't stuck against an exterior face.
		if (tip_time_start > ct) {
			

			// create a tip on the element face
			Tip * next = new Tip();
			Segment * seg = new Segment();
			// the tip will appear at the time it takes to reach the face.
			next->time = tip_time_start;
			// assign the tip to this element.
			next->angio_element = angio_element;
			// update species.
			next->Species = active_tip->Species;
			active_tip->Species.clear();

			//still need to update the local position of the tip. 
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
			for (int i = 0; i < angio_element->adjacency_list.size(); i++)
			{
				AngioElement * ang_elem = angio_element->adjacency_list[i];
				if (ang_elem)
				{
					FESolidDomain * sd = dynamic_cast<FESolidDomain*>(ang_elem->_elem->GetMeshPartition());
					if (sd)
					{
						double r[3];//new natural coordinates
						bool natc_proj = m_pangio->ProjectToElement(*ang_elem->_elem, pos, mesh, r);
						if (natc_proj &&  m_pangio->IsInBounds(ang_elem->_elem, r, bounds_tolerance))
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
				//need some way to choose the correct element to continue the tip in
				int index = SelectNextTip(possible_locations, possible_local_coordinates, next, dt, buffer_index, mesh, min_scale_factor, min_angle);
				if (index != -1)
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
				angio_element->final_active_tips.push_back(next);
				angio_element->next_tips.at(angio_element).push_back(next);
				//next->growth_velocity = 0;
			}
		}
		// if the tips weren't growing re-add them to the active list for now. In the future may want to have this be the case to handle model boundaries and other stuff.
		else
		{ 
			// add the tip that hit an external face to the active tips
			angio_element->final_active_tips.push_back(active_tip);
			angio_element->next_tips.at(angio_element).push_back(active_tip); 
		}
	}
	// generally not hit unless there is large matrix deformation

	else{
		
		angio_element->final_active_tips.push_back(active_tip);
		angio_element->next_tips.at(angio_element).push_back(active_tip);
	}
}

void FEAngioMaterial::ProtoGrowthInElement(double end_time, Tip * active_tip, int source_index, int buffer_index, double min_scale_factor, double bounds_tolerance, double min_angle)
{
	////////// Initialization //////////

	// get the current angio element and make sure it's active
	auto angio_element = active_tip->angio_element;

	assert(active_tip);
	// get the mesh, set the min segment length, and set the next buffer (change to write buffer from read).
	auto mesh = m_pangio->GetMesh();
	const double min_segm_len = 5;
	int next_buffer_index = (buffer_index + 1) % 2;
	// calculate the change in time
	double dt = end_time - active_tip->time;
	// determine the length the segment grows.
	double grow_vel = active_tip->GetProtoGrowthLength();
	//double grow_vel = angio_element->_angio_mat->GetInitialVelocity(angio_element);
	if (grow_vel < 0)
	{
		return;
	}
	assert(grow_vel > 0.0);
	double grow_len = grow_vel*dt;

	////////// Determine the growth length and new position //////////

	// get the local position of the tip
	double Gr[FESolidElement::MAX_NODES];
	double Gs[FESolidElement::MAX_NODES];
	double Gt[FESolidElement::MAX_NODES];
	// get the local position of the tip
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
	// determine the inverse operation to go from global to natural
	mat3d global_to_natc = natc_to_global.inverse();
	bool continue_growth = true;

	// get the growth direction
	// determine the alpha (mix value)
	double alpha = proto_cm_manager->ApplyModifiers(dt, active_tip->angio_element, active_tip->GetLocalPosition(), mesh);
	vec3d proto_psc_dir = proto_psc_manager->ApplyModifiers(vec3d(1, 0, 0), active_tip->angio_element, active_tip->GetLocalPosition(), active_tip->GetDirection(mesh), mesh);
	vec3d proto_pdd_dir = proto_pdd_manager->ApplyModifiers(vec3d(1, 0, 0), active_tip->angio_element, active_tip->GetLocalPosition(), active_tip->initial_fragment_id, buffer_index, continue_growth, proto_psc_dir, alpha, mesh, m_pangio);
	
	// assign the global direction from the persistence component
	vec3d global_dir = mix_method->ApplyMix(proto_psc_dir, proto_pdd_dir, alpha);
	global_dir.unit();

	// get the new position in global coordinates
	vec3d global_pos;
	double H[FESolidElement::MAX_NODES];
	angio_element->_elem->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		global_pos += mesh->Node(angio_element->_elem->m_node[j]).m_rt* H[j];
	}

	vec3d nat_dir = global_to_natc * global_dir;

	////////// Calculate variables to determine which element to grow into //////////

	double factor;
	// Determine if the ray cast run along the face of an element. If you see vessels bouncing off unexposed faces this is likely the issue.
	//This also updates "factor" which scales the length allowed in the element. It is the shortest of the 3 directions that an element can grow before exceeding the bounds while also staying above the min_scale_factor
	//If it fails that means that all possible factors are less than the min_scale_factor.
	bool proj_success = m_pangio->ScaleFactorToProjectToNaturalCoordinates(active_tip->angio_element->_elem, nat_dir, local_pos, factor, min_scale_factor); // TODO: HERE
	// Calculate the position on the face that the tip will grow to.
	vec3d possible_natc = local_pos + (nat_dir*factor);
	// determine the hypothetical length to grow
	double possible_grow_length = m_pangio->InElementLength(active_tip->angio_element->_elem, local_pos, possible_natc);
	
	////////// Determine which element to grow into //////////
	
	// if the vessel is going to stay in this element (or hit the face) and a face to grow to was found:
	if ((possible_grow_length >= grow_len) && proj_success)
	{
		//determine the new natural coordinate position in this element. Not sure that this calculation needs to be this convoluted.
		vec3d real_natc = local_pos + (nat_dir * (factor * (grow_len / possible_grow_length)));
		//add the segment and add the new tip to the current angio element
		Tip * next = new Tip(active_tip, mesh);
		Segment * seg = new Segment();
		next->time = end_time;
		// the new tip is assigned to the current element.
		next->angio_element = angio_element;
		// assign the new position to the new tip.
		next->SetLocalPosition(real_natc, mesh);
		// assign the new tip as the parent of the segment.
		next->parent = seg;
		// assign the current element as location where growth began.
		next->face = angio_element;
		// assign the velocity that the tip grew at
		//SL!
		//growth velocity must be zero to work with grown segments stress policy. Otherwise it would be grow_vel. TODO: Maybe make it use different methods based on the policy?
		next->growth_velocity = 0;
		next->SetProtoGrowthLength(active_tip);
		// copy the parent fragment id from the original tip.
		next->initial_fragment_id = active_tip->initial_fragment_id;
		// ensure that if the element is a tet that it does not exceed the element boundaries. This should already have been done in ScaleFactor calculation.
		assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G4 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);
		assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G1 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);
		// add the new tip to the element so the element manages the tip.
		angio_element->next_tips.at(angio_element).push_back(next);
		// make the active tip the back of this segment and make next the front.

		seg->back = active_tip;
		seg->front = next;
		// assign the segment's parent seg as the same as the active tip's parent seg
		if (active_tip->parent)
		{
			seg->parent = active_tip->parent;
		}

		// add the segment to the recent segments vector
		angio_element->recent_segments.push_back(seg);
		// calculate the segment length of seg
		angio_element->refernce_frame_segment_length += seg->RefLength(mesh);
		// add next to the final active tips
		angio_element->final_active_tips.push_back(next);
		assert(next->angio_element);
	}
	else if (proj_success)
	{
		//the segment only grows for a portion of dt. This portion is the amount needed to hit the face.
		//grow len should always be nonzero so this division should be okay. This will need to be accounted for in the future if we explicitly model sprouting.
		double tip_time_start = active_tip->time + (possible_grow_length / grow_len) * dt;
		assert(tip_time_start < end_time);
		//grow len should always be nonzero so this division should be okay
		//this segment is growing into a new element do the work to do this

		//next is just the end of the segment not the begining of the next segment
	
		// create a tip on the element face
		Tip * next = new Tip();
		Segment * seg = new Segment();
		// the tip will appear at the time it takes to reach the face.
		next->time = tip_time_start;
		// assign the tip to this element.
		next->angio_element = angio_element;
		// update species.
		next->Species = active_tip->Species;
		active_tip->Species.clear();

		//still need to update the local position of the tip
		next->SetLocalPosition(local_pos + (nat_dir * possible_grow_length), mesh);

		next->parent = seg;
		next->face = angio_element;
		next->initial_fragment_id = active_tip->initial_fragment_id;
		//growth velocity must be zero to work with grown segments stress policy
		next->growth_velocity = 0;
		next->SetProtoGrowthLength(active_tip);

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
		for (int i = 0; i < angio_element->adjacency_list.size(); i++) // TODO: HERE //!SL
		{
			AngioElement * ang_elem = angio_element->adjacency_list[i];
			if (ang_elem)
			{
				FESolidDomain * sd = dynamic_cast<FESolidDomain*>(ang_elem->_elem->GetMeshPartition());
				if (sd)
				{
					double r[3];//new natural coordinates
					bool natc_proj = m_pangio->ProjectToElement(*ang_elem->_elem, pos, mesh, r);
					bool in_bounds = m_pangio->IsInBounds(ang_elem->_elem, r, bounds_tolerance);
					if (natc_proj && in_bounds)
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
			int index = ProtoSelectNextTip(possible_locations, possible_local_coordinates, next, dt, buffer_index , mesh, min_scale_factor, min_angle);
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
		// if there are no adjacent elements to grow into then 
		// currently this just keeps the tip active.
		else
		{
			// in this element assign this tip for evaluation on the next step. We will assign it at the current element by pushing back the active tip.
			// let's grow it somewhere
			angio_element->final_active_tips.push_back(next);
			angio_element->next_tips.at(angio_element).push_back(active_tip);
			//next->growth_velocity = 0;
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
		global_dir.unit();
		vec3d possible_dir = possible_locations[i]->_angio_mat->pdd_manager->ApplyModifiers({ 1,0,0 }, possible_locations[i], possible_local_coordinates[i], tip->initial_fragment_id, buffer, continue_growth, psc_dir, alpha, mesh, m_pangio);
		double Gr[FESolidElement::MAX_NODES];
		double Gs[FESolidElement::MAX_NODES];
		double Gt[FESolidElement::MAX_NODES];
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

int FEAngioMaterial::ProtoSelectNextTip(std::vector<AngioElement*> & possible_locations, std::vector<vec3d> & possible_local_coordinates, Tip* tip, double dt, int buffer, FEMesh* mesh, double min_scale_factor, double min_angle)
{
	assert(possible_locations.size() == possible_local_coordinates.size());
	if (possible_locations.size() == 1)
		return 0;
	auto dir = tip->GetDirection(mesh);
	std::vector<double> angles;
	bool continue_growth = true;
	for (int i = 0; i < possible_locations.size(); i++)
	{
		double alpha = proto_cm_manager->ApplyModifiers(dt, possible_locations[i], possible_local_coordinates[i], mesh);
		vec3d proto_psc_dir = proto_psc_manager->ApplyModifiers(vec3d(1, 0, 0), possible_locations[i], possible_local_coordinates[i], tip->GetDirection(mesh), mesh);
		vec3d proto_pdd_dir = proto_pdd_manager->ApplyModifiers(vec3d(1, 0, 0), possible_locations[i], possible_local_coordinates[i], tip->initial_fragment_id, buffer, continue_growth, proto_psc_dir, alpha, mesh, m_pangio);
		if (!continue_growth)
		{
			return i;
		}
		vec3d global_dir = mix_method->ApplyMix(proto_psc_dir, proto_pdd_dir, (alpha));
		global_dir.unit();
		vec3d possible_dir = possible_locations[i]->_angio_mat->proto_pdd_manager->ApplyModifiers({ 1,0,0 }, possible_locations[i], possible_local_coordinates[i], tip->initial_fragment_id, buffer, continue_growth, proto_psc_dir, alpha, mesh, m_pangio);
		double Gr[FESolidElement::MAX_NODES];
		double Gs[FESolidElement::MAX_NODES];
		double Gt[FESolidElement::MAX_NODES];
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
		if (proj_success && factor > min_scale_factor && possible_dir*dir >= min_angle)
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
	for (int i = 0; i < angles.size(); i++)
	{
		if (angles[i] > min)
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
	return new FEAngioMaterialPoint((FEElasticMaterial::CreateMaterialPointData()), common_properties->vessel_material->CreateMaterialPointData(), matrix_material->CreateMaterialPointData());
}